//! `HmcApp`: the egui application. Owns the parameter state, the
//! background worker, and the most-recent experiment result for plotting.

use std::time::Duration;

use eframe::egui::{self, CollapsingHeader, RichText, ScrollArea, Slider};
use egui_plot::{Corner, Legend, Line, LineStyle, Plot, PlotPoints, Points};

use crate::export::{write_pdfs, ExportPaths};
use crate::potentials::Potential;
use crate::run::{ExperimentResult, SimParams, Stage};

use super::worker::{Progress, SimWorker, WorkerMessage};

/// Status string for the bottom bar.
enum JobStatus {
    Idle,
    Running(String),
    Done { elapsed_ms: u128 },
    Exported(String),
    Error(String),
}

pub struct HmcApp {
    params: SimParams,
    last_submitted: Option<SimParams>,
    worker: SimWorker,
    in_flight: bool,
    auto_rerun: bool,
    result: Option<ExperimentResult>,
    last_elapsed_ms: Option<u128>,
    status: JobStatus,
    enable_accept_temp: bool,
    accept_temp_value: f64,
    // Potential-specific scratch space (currently only DoubleWell).
    x_left: f64,
    x_right: f64,
    /// Which mass index in `result.momentum_scan_hmc` the momentum
    /// histogram is currently showing. Bounded against the result each frame.
    momentum_mass_index: usize,
}

impl Default for HmcApp {
    fn default() -> Self {
        let params = SimParams::default();
        let (x_left, x_right) = match params.potential {
            Potential::DoubleWell { x_left, x_right } => (x_left, x_right),
            _ => (-1.0, 1.0),
        };
        Self {
            params,
            last_submitted: None,
            worker: SimWorker::spawn(),
            in_flight: false,
            auto_rerun: false,
            result: None,
            last_elapsed_ms: None,
            status: JobStatus::Idle,
            enable_accept_temp: false,
            accept_temp_value: 1.0,
            x_left,
            x_right,
            momentum_mass_index: 0,
        }
    }
}

impl HmcApp {
    pub fn new(_cc: &eframe::CreationContext<'_>) -> Self {
        Self::default()
    }

    fn submit_current(&mut self) {
        // Sync UI-only state back into params.
        self.params.acceptance_temp = if self.enable_accept_temp {
            Some(self.accept_temp_value)
        } else {
            None
        };
        self.params.potential = Potential::DoubleWell {
            x_left: self.x_left,
            x_right: self.x_right,
        };
        self.in_flight = true;
        self.status = JobStatus::Running("queued".into());
        self.last_submitted = Some(self.params.clone());
        self.worker.submit(self.params.clone());
    }

    fn ingest_worker(&mut self) {
        for msg in self.worker.try_drain() {
            match msg {
                WorkerMessage::Progress(Progress::Started) => {
                    self.status = JobStatus::Running("starting".into());
                }
                WorkerMessage::Progress(Progress::Stage(stage)) => {
                    self.status = JobStatus::Running(stage_label(&stage));
                }
                WorkerMessage::Finished { result, elapsed_ms } => {
                    self.result = Some(*result);
                    self.last_elapsed_ms = Some(elapsed_ms);
                    self.in_flight = false;
                    self.status = JobStatus::Done { elapsed_ms };
                }
            }
        }
    }

    fn parameters_changed(&self) -> bool {
        match &self.last_submitted {
            None => true,
            Some(prev) => {
                let mut current = self.params.clone();
                current.acceptance_temp = if self.enable_accept_temp {
                    Some(self.accept_temp_value)
                } else {
                    None
                };
                current.potential = Potential::DoubleWell {
                    x_left: self.x_left,
                    x_right: self.x_right,
                };
                &current != prev
            }
        }
    }
}

fn stage_label(stage: &Stage) -> String {
    match stage {
        Stage::Hmc { chain, total_chains, done, total } => format!(
            "HMC chain {}/{}: {}/{}",
            chain + 1,
            total_chains,
            done,
            total
        ),
        Stage::Mh { chain, total_chains, done, total } => format!(
            "MH chain {}/{}: {}/{}",
            chain + 1,
            total_chains,
            done,
            total
        ),
        Stage::MomentumScanHmc => "scanning masses (HMC)".into(),
        Stage::MomentumScanMh => "scanning masses (MH)".into(),
        Stage::Histogram => "building histograms".into(),
        Stage::Done => "rendering".into(),
    }
}

impl eframe::App for HmcApp {
    fn ui(&mut self, ui: &mut egui::Ui, _frame: &mut eframe::Frame) {
        self.ingest_worker();

        if self.in_flight {
            ui.ctx().request_repaint_after(Duration::from_millis(50));
        }

        egui::Panel::left("controls")
            .resizable(true)
            .default_size(340.0)
            .show_inside(ui, |ui| {
                ScrollArea::vertical().show(ui, |ui| {
                    ui.heading("Parameters");

                    CollapsingHeader::new("Sampler")
                        .default_open(true)
                        .show(ui, |ui| {
                            slider_f64(ui, "mass", &mut self.params.mass, 0.05..=10.0);
                            slider_f64(ui, "beta", &mut self.params.beta, 0.05..=10.0);
                            ui.add(
                                Slider::new(&mut self.params.step_size, 1e-4..=1.0)
                                    .text("step_size (ε)")
                                    .logarithmic(true),
                            );
                            slider_f64(
                                ui,
                                "avg_sim_time",
                                &mut self.params.avg_sim_time,
                                0.01..=20.0,
                            );
                            slider_usize(ui, "dimensions", &mut self.params.dimensions, 1..=8);
                            ui.horizontal(|ui| {
                                ui.checkbox(&mut self.enable_accept_temp, "acceptance_temp");
                                if self.enable_accept_temp {
                                    ui.add(
                                        Slider::new(&mut self.accept_temp_value, 0.05..=5.0)
                                            .text("T_a"),
                                    );
                                }
                            });
                        });

                    CollapsingHeader::new("Chain")
                        .default_open(true)
                        .show(ui, |ui| {
                            slider_usize(
                                ui,
                                "iterations",
                                &mut self.params.iterations,
                                10..=20_000,
                            );
                            slider_usize(ui, "num_chains", &mut self.params.num_chains, 1..=16);
                            ui.label("Sample window for momentum check (% of chain):");
                            ui.horizontal(|ui| {
                                ui.add(
                                    Slider::new(
                                        &mut self.params.burn_in_lo_percent,
                                        0.0..=100.0,
                                    )
                                    .text("lo"),
                                );
                            });
                            ui.horizontal(|ui| {
                                ui.add(
                                    Slider::new(
                                        &mut self.params.burn_in_hi_percent,
                                        0.0..=100.0,
                                    )
                                    .text("hi"),
                                );
                            });
                            if self.params.burn_in_lo_percent
                                > self.params.burn_in_hi_percent
                            {
                                std::mem::swap(
                                    &mut self.params.burn_in_lo_percent,
                                    &mut self.params.burn_in_hi_percent,
                                );
                            }
                            ui.label("Initial-position box (broadcast over dims):");
                            ui.horizontal(|ui| {
                                ui.add(
                                    egui::DragValue::new(&mut self.params.init_box.0)
                                        .speed(0.05)
                                        .prefix("lo: "),
                                );
                                ui.add(
                                    egui::DragValue::new(&mut self.params.init_box.1)
                                        .speed(0.05)
                                        .prefix("hi: "),
                                );
                            });
                            if self.params.init_box.0 > self.params.init_box.1 {
                                std::mem::swap(
                                    &mut self.params.init_box.0,
                                    &mut self.params.init_box.1,
                                );
                            }
                        });

                    CollapsingHeader::new("MH baseline")
                        .default_open(false)
                        .show(ui, |ui| {
                            slider_f64(
                                ui,
                                "proposal_std",
                                &mut self.params.proposal_std,
                                0.01..=5.0,
                            );
                        });

                    CollapsingHeader::new("Mass scan")
                        .default_open(false)
                        .show(ui, |ui| {
                            slider_f64(
                                ui,
                                "mass_upper_bound",
                                &mut self.params.mass_upper_bound,
                                0.1..=10.0,
                            );
                            slider_usize(
                                ui,
                                "num_masses",
                                &mut self.params.num_masses,
                                2..=100,
                            );
                            slider_usize(
                                ui,
                                "check_iterations",
                                &mut self.params.check_iterations,
                                10..=5_000,
                            );
                        });

                    CollapsingHeader::new("Potential (double well)")
                        .default_open(false)
                        .show(ui, |ui| {
                            slider_f64(ui, "x_left", &mut self.x_left, -5.0..=5.0);
                            slider_f64(ui, "x_right", &mut self.x_right, -5.0..=5.0);
                        });

                    ui.separator();

                    ui.horizontal(|ui| {
                        let run_label = if self.in_flight { "Running…" } else { "Run" };
                        let run = ui
                            .add_enabled(!self.in_flight, egui::Button::new(run_label))
                            .clicked();
                        if run {
                            self.submit_current();
                        }
                        ui.checkbox(&mut self.auto_rerun, "Auto-rerun on change");
                    });

                    if self.auto_rerun && !self.in_flight && self.parameters_changed() {
                        self.submit_current();
                    }

                    ui.horizontal(|ui| {
                        let can_export = self.result.is_some();
                        if ui
                            .add_enabled(can_export, egui::Button::new("Export PDFs"))
                            .clicked()
                        {
                            if let Some(result) = &self.result {
                                let mut paths = ExportPaths::default();
                                paths.momentum_mass_index = Some(self.momentum_mass_index);
                                match write_pdfs(result, &paths) {
                                    Ok(_) => {
                                        self.status = JobStatus::Exported(format!(
                                            "wrote {}, {}, {}, {}",
                                            paths.position_histogram,
                                            paths.mh_position_histogram,
                                            paths.momentum_histogram,
                                            paths.momentum_scan,
                                        ));
                                    }
                                    Err(e) => {
                                        self.status =
                                            JobStatus::Error(format!("export failed: {e}"));
                                    }
                                }
                            }
                        }
                    });
                });
            });

        egui::Panel::bottom("status").show_inside(ui, |ui| match &self.status {
            JobStatus::Idle => {
                ui.label("Idle. Adjust parameters and press Run.");
            }
            JobStatus::Running(label) => {
                ui.horizontal(|ui| {
                    ui.spinner();
                    ui.label(format!("Running: {label}"));
                });
            }
            JobStatus::Done { elapsed_ms } => {
                ui.label(format!(
                    "Done in {:.2}s. {} HMC samples, {} MH samples.",
                    *elapsed_ms as f64 / 1000.0,
                    self.result.as_ref().map(|r| r.hmc_samples.len()).unwrap_or(0),
                    self.result.as_ref().map(|r| r.mh_samples.len()).unwrap_or(0),
                ));
            }
            JobStatus::Exported(msg) => {
                ui.label(RichText::new(msg).strong());
            }
            JobStatus::Error(msg) => {
                ui.label(RichText::new(msg).color(egui::Color32::RED));
            }
        });

        egui::CentralPanel::default().show_inside(ui, |ui| {
            if let Some(result) = &self.result {
                // Clamp the momentum-mass selector against the current result.
                let max_idx = result.momentum_scan_hmc.len().saturating_sub(1);
                if self.momentum_mass_index > max_idx {
                    self.momentum_mass_index = max_idx;
                }
                draw_plots(ui, result, &mut self.momentum_mass_index);
            } else {
                ui.vertical_centered(|ui| {
                    ui.add_space(60.0);
                    ui.heading("No simulation run yet");
                    ui.label("Press Run on the left panel.");
                });
            }
        });
    }
}

fn draw_plots(ui: &mut egui::Ui, result: &ExperimentResult, momentum_mass_index: &mut usize) {
    let avail = ui.available_size();
    let plot_height = (avail.y / 4.0 - 36.0).max(110.0);

    CollapsingHeader::new("HMC position histogram + V(x)")
        .default_open(true)
        .show(ui, |ui| {
            Plot::new("hmc_position_histogram")
                .height(plot_height)
                .legend(Legend::default())
                .show(ui, |plot_ui| {
                    plot_ui.line(
                        Line::new(
                            "HMC hist (scaled)",
                            to_plot_points(&result.histogram_bins),
                        )
                        .color(egui::Color32::from_rgb(70, 130, 180)),
                    );
                    plot_ui.line(
                        Line::new("V(x)", to_plot_points(&result.potential_curve))
                            .color(egui::Color32::from_rgb(255, 127, 80)),
                    );
                });
        });

    CollapsingHeader::new("MH position histogram + V(x)")
        .default_open(true)
        .show(ui, |ui| {
            Plot::new("mh_position_histogram")
                .height(plot_height)
                .legend(Legend::default())
                .show(ui, |plot_ui| {
                    plot_ui.line(
                        Line::new(
                            "MH hist (scaled)",
                            to_plot_points(&result.mh_histogram_bins),
                        )
                        .color(egui::Color32::from_rgb(60, 179, 113)),
                    );
                    plot_ui.line(
                        Line::new("V(x)", to_plot_points(&result.potential_curve))
                            .color(egui::Color32::from_rgb(255, 127, 80)),
                    );
                });
        });

    CollapsingHeader::new("Momentum histogram (per mass)")
        .default_open(true)
        .show(ui, |ui| {
            let max_idx = result.momentum_scan_hmc.len().saturating_sub(1);
            ui.horizontal(|ui| {
                ui.label("Mass index");
                ui.add(egui::Slider::new(momentum_mass_index, 0..=max_idx).integer());
                if let Some((mass, var)) = result.momentum_scan_hmc.get(*momentum_mass_index) {
                    let expected_var = mass / result.params.beta;
                    ui.label(format!(
                        "m = {:.3}  |  var = {:.3}  |  expected m/β = {:.3}",
                        mass, var, expected_var
                    ));
                }
            });

            Plot::new("momentum_hist")
                .height(plot_height)
                .legend(Legend::default())
                .show(ui, |plot_ui| {
                    let hmc_momenta = result.hmc_momenta_per_mass.get(*momentum_mass_index);
                    let mh_momenta = result.mh_momenta_per_mass.get(*momentum_mass_index);
                    let Some(&(mass, _)) = result.momentum_scan_hmc.get(*momentum_mass_index)
                    else {
                        return;
                    };

                    // Shared x-range so both histograms (and the expected
                    // Gaussian) sit on the same axis without auto-zoom jumps.
                    let mut combined: Vec<f64> = Vec::new();
                    if let Some(h) = hmc_momenta {
                        combined.extend_from_slice(h);
                    }
                    if let Some(m) = mh_momenta {
                        combined.extend_from_slice(m);
                    }
                    if combined.is_empty() {
                        return;
                    }

                    let hmc_marks = hmc_momenta
                        .map(|m| crate::export::histogram_centers(m, 50))
                        .unwrap_or_default();
                    let mh_marks = mh_momenta
                        .map(|m| crate::export::histogram_centers(m, 50))
                        .unwrap_or_default();

                    let expected = crate::export::expected_gaussian_curve(
                        mass,
                        result.params.beta,
                        &combined,
                    );
                    let bin_peak = hmc_marks
                        .iter()
                        .chain(mh_marks.iter())
                        .map(|(_, y)| *y)
                        .fold(0.0_f64, f64::max);
                    let expected_peak = expected
                        .iter()
                        .map(|(_, y)| *y)
                        .fold(0.0_f64, f64::max)
                        .max(1e-12);
                    let scale = bin_peak / expected_peak;
                    let scaled_expected: Vec<(f64, f64)> = expected
                        .into_iter()
                        .map(|(x, y)| (x, y * scale))
                        .collect();

                    if !hmc_marks.is_empty() {
                        plot_ui.points(
                            Points::new("HMC momenta", to_plot_points(&hmc_marks))
                                .radius(3.0)
                                .color(egui::Color32::from_rgb(70, 130, 180)),
                        );
                    }
                    if !mh_marks.is_empty() {
                        plot_ui.points(
                            Points::new("MH momenta", to_plot_points(&mh_marks))
                                .radius(3.0)
                                .color(egui::Color32::from_rgb(60, 179, 113)),
                        );
                    }
                    plot_ui.line(
                        Line::new("Expected N(0, √(m/β))", to_plot_points(&scaled_expected))
                            .color(egui::Color32::from_rgb(255, 127, 80)),
                    );

                    // Gaussian fits to the empirical momenta (mean 0, σ from
                    // the population std). Scaled to the same peak as the markers.
                    if let Some(m) = hmc_momenta {
                        if let Some(fit) = crate::export::gaussian_fit_curve(m, &combined, scale) {
                            plot_ui.line(
                                Line::new("HMC fit", to_plot_points(&fit))
                                    .style(LineStyle::dashed_loose())
                                    .color(egui::Color32::from_rgb(70, 130, 180)),
                            );
                        }
                    }
                    if let Some(m) = mh_momenta {
                        if let Some(fit) = crate::export::gaussian_fit_curve(m, &combined, scale) {
                            plot_ui.line(
                                Line::new("MH fit", to_plot_points(&fit))
                                    .style(LineStyle::dashed_loose())
                                    .color(egui::Color32::from_rgb(60, 179, 113)),
                            );
                        }
                    }
                });
        });

    CollapsingHeader::new("Momentum variance vs mass (expected slope = 1/β)")
        .default_open(true)
        .show(ui, |ui| {
            Plot::new("momentum_scan")
                .height(plot_height)
                .legend(Legend::default().position(Corner::RightBottom))
                .show(ui, |plot_ui| {
                    plot_ui.line(
                        Line::new("HMC", to_plot_points(&result.momentum_scan_hmc))
                            .color(egui::Color32::from_rgb(70, 130, 180)),
                    );
                    plot_ui.line(
                        Line::new("MH", to_plot_points(&result.momentum_scan_mh))
                            .color(egui::Color32::from_rgb(60, 179, 113)),
                    );
                    plot_ui.line(
                        Line::new(
                            "Expected m/β",
                            to_plot_points(&result.momentum_scan_expected),
                        )
                        .color(egui::Color32::from_rgb(255, 127, 80)),
                    );
                });
        });
}

fn to_plot_points(data: &[(f64, f64)]) -> PlotPoints<'static> {
    PlotPoints::from_iter(data.iter().map(|(x, y)| [*x, *y]))
}


fn slider_f64(
    ui: &mut egui::Ui,
    label: &str,
    value: &mut f64,
    range: std::ops::RangeInclusive<f64>,
) -> egui::Response {
    ui.add(Slider::new(value, range).text(label))
}

fn slider_usize(
    ui: &mut egui::Ui,
    label: &str,
    value: &mut usize,
    range: std::ops::RangeInclusive<usize>,
) -> egui::Response {
    ui.add(Slider::new(value, range).text(label))
}
