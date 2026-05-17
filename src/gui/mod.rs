//! egui-based interactive frontend for the HMC + MH experiments.
//!
//! The app spawns a worker thread that runs `run_experiment` on demand;
//! the UI itself never blocks on the simulation. Plots are drawn live
//! with `egui_plot`. An "Export PDFs" button hands the latest result to
//! the existing `kuva`-based pipeline in `crate::export`.

pub mod app;
pub mod worker;

pub use app::HmcApp;
