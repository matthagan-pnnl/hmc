//! Background simulation worker.
//!
//! Owns a thread that loops on incoming `SimParams` requests, runs the
//! experiment, and ships progress + the final `ExperimentResult` back
//! over channels. The UI only ever calls `submit` and `try_drain`.

use std::sync::mpsc::{self, Receiver, Sender, TryRecvError};
use std::thread;
use std::time::Instant;

use crate::run::{run_experiment, ExperimentResult, SimParams, Stage};

#[derive(Debug, Clone)]
pub enum Progress {
    Started,
    Stage(Stage),
}

pub enum WorkerMessage {
    Progress(Progress),
    Finished { result: Box<ExperimentResult>, elapsed_ms: u128 },
}

pub struct SimWorker {
    request_tx: Sender<SimParams>,
    inbox: Receiver<WorkerMessage>,
}

impl SimWorker {
    pub fn spawn() -> Self {
        let (request_tx, request_rx) = mpsc::channel::<SimParams>();
        let (msg_tx, inbox) = mpsc::channel::<WorkerMessage>();

        thread::Builder::new()
            .name("hmc-worker".into())
            .spawn(move || worker_loop(request_rx, msg_tx))
            .expect("spawn worker thread");

        Self { request_tx, inbox }
    }

    pub fn submit(&self, params: SimParams) {
        // The receiver lives for the lifetime of the worker thread, which
        // we never join — so send failures are non-recoverable bugs.
        let _ = self.request_tx.send(params);
    }

    /// Drain all messages currently in the inbox.
    pub fn try_drain(&self) -> Vec<WorkerMessage> {
        let mut out = Vec::new();
        loop {
            match self.inbox.try_recv() {
                Ok(msg) => out.push(msg),
                Err(TryRecvError::Empty) | Err(TryRecvError::Disconnected) => break,
            }
        }
        out
    }
}

fn worker_loop(request_rx: Receiver<SimParams>, msg_tx: Sender<WorkerMessage>) {
    while let Ok(mut params) = request_rx.recv() {
        // Coalesce: if the UI rapidly pushed several requests, only
        // honor the most recent one.
        while let Ok(newer) = request_rx.try_recv() {
            params = newer;
        }

        let _ = msg_tx.send(WorkerMessage::Progress(Progress::Started));
        let started = Instant::now();
        let tx_for_stages = msg_tx.clone();
        let result = run_experiment(&params, |stage| {
            let _ = tx_for_stages.send(WorkerMessage::Progress(Progress::Stage(stage)));
        });
        let elapsed_ms = started.elapsed().as_millis();
        let _ = msg_tx.send(WorkerMessage::Finished {
            result: Box::new(result),
            elapsed_ms,
        });
    }
}
