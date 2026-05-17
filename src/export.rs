//! PDF export of an `ExperimentResult` via the `kuva` plotting crate.
//! Writes three PDFs: position histogram, per-mass momentum histogram,
//! and the momentum-std-vs-mass scan.

use kuva::prelude::*;

use crate::run::ExperimentResult;

pub struct ExportPaths {
    pub position_histogram: String,
    pub mh_position_histogram: String,
    pub momentum_histogram: String,
    pub momentum_scan: String,
    /// Which mass index to histogram momenta for. Defaults to the index
    /// whose mass is closest to `result.params.mass`.
    pub momentum_mass_index: Option<usize>,
}

impl Default for ExportPaths {
    fn default() -> Self {
        Self {
            position_histogram: "hmc_histogram.pdf".into(),
            mh_position_histogram: "mh_histogram.pdf".into(),
            momentum_histogram: "momentum_histogram.pdf".into(),
            momentum_scan: "momentum_std_vs_mass.pdf".into(),
            momentum_mass_index: None,
        }
    }
}

pub fn write_pdfs(result: &ExperimentResult, paths: &ExportPaths) -> std::io::Result<()> {
    // Position histogram with potential overlay.
    let position_hist_line = LinePlot::new()
        .with_data(result.histogram_bins.clone())
        .with_color("steelblue");
    let potential_line = LinePlot::new()
        .with_data(result.potential_curve.clone())
        .with_color("coral");
    let pos_plots: Vec<Plot> = vec![position_hist_line.into(), potential_line.into()];
    let pos_layout = Layout::auto_from_plots(&pos_plots)
        .with_title("Histogram of Chain Samples with Potential Function")
        .with_x_label("Position")
        .with_y_label("Count / Potential Energy");
    let pdf = render_to_pdf(pos_plots, pos_layout).expect("render position histogram");
    std::fs::write(&paths.position_histogram, pdf)?;

    // MH position histogram with potential overlay.
    let mh_hist_line = LinePlot::new()
        .with_data(result.mh_histogram_bins.clone())
        .with_color("mediumseagreen");
    let mh_pot_line = LinePlot::new()
        .with_data(result.potential_curve.clone())
        .with_color("coral");
    let mh_plots: Vec<Plot> = vec![mh_hist_line.into(), mh_pot_line.into()];
    let mh_layout = Layout::auto_from_plots(&mh_plots)
        .with_title("MH Histogram of Chain Samples with Potential Function")
        .with_x_label("Position")
        .with_y_label("Count / Potential Energy");
    let pdf = render_to_pdf(mh_plots, mh_layout).expect("render MH histogram");
    std::fs::write(&paths.mh_position_histogram, pdf)?;

    // Momentum histogram at the requested mass.
    let mass_idx = paths
        .momentum_mass_index
        .unwrap_or_else(|| nearest_mass_index(result, result.params.mass));
    if let Some(&(mass, _var)) = result.momentum_scan_hmc.get(mass_idx) {
        let hmc_m = result.hmc_momenta_per_mass.get(mass_idx);
        let mh_m = result.mh_momenta_per_mass.get(mass_idx);
        let mut combined: Vec<f64> = Vec::new();
        if let Some(h) = hmc_m {
            combined.extend_from_slice(h);
        }
        if let Some(m) = mh_m {
            combined.extend_from_slice(m);
        }
        if !combined.is_empty() {
            let hmc_marks = hmc_m.map(|m| histogram_centers(m, 50)).unwrap_or_default();
            let mh_marks = mh_m.map(|m| histogram_centers(m, 50)).unwrap_or_default();
            let expected = expected_gaussian_curve(mass, result.params.beta, &combined);
            let bin_peak = peak_height(&hmc_marks).max(peak_height(&mh_marks));
            let scale = bin_peak / peak_height(&expected).max(1e-12);
            let scaled_expected: Vec<(f64, f64)> =
                expected.into_iter().map(|(x, y)| (x, y * scale)).collect();
            let mut mom_plots: Vec<Plot> = Vec::new();
            if !hmc_marks.is_empty() {
                mom_plots.push(
                    ScatterPlot::new()
                        .with_data(hmc_marks)
                        .with_color("steelblue")
                        .with_size(4.0)
                        .with_legend("HMC momenta")
                        .into(),
                );
            }
            if !mh_marks.is_empty() {
                mom_plots.push(
                    ScatterPlot::new()
                        .with_data(mh_marks)
                        .with_color("mediumseagreen")
                        .with_size(4.0)
                        .with_legend("MH momenta")
                        .into(),
                );
            }
            mom_plots.push(
                LinePlot::new()
                    .with_data(scaled_expected)
                    .with_color("coral")
                    .with_legend(format!("N(0, √(m/β)) at m={:.3}", mass).as_str())
                    .into(),
            );

            // Empirical Gaussian fits, dashed, same color as the markers.
            if let Some(samples) = hmc_m {
                if let Some(fit) = gaussian_fit_curve(samples, &combined, scale) {
                    mom_plots.push(
                        LinePlot::new()
                            .with_data(fit)
                            .with_color("steelblue")
                            .with_dashed()
                            .with_legend("HMC fit")
                            .into(),
                    );
                }
            }
            if let Some(samples) = mh_m {
                if let Some(fit) = gaussian_fit_curve(samples, &combined, scale) {
                    mom_plots.push(
                        LinePlot::new()
                            .with_data(fit)
                            .with_color("mediumseagreen")
                            .with_dashed()
                            .with_legend("MH fit")
                            .into(),
                    );
                }
            }
            let mom_layout = Layout::auto_from_plots(&mom_plots)
                .with_title(&format!("Momentum Histogram (mass = {:.3})", mass))
                .with_x_label("Momentum")
                .with_y_label("Count");
            let pdf = render_to_pdf(mom_plots, mom_layout).expect("render momentum histogram");
            std::fs::write(&paths.momentum_histogram, pdf)?;
        }
    }

    // Momentum-std-vs-mass scan.
    let hmc_line = LinePlot::new()
        .with_data(result.momentum_scan_hmc.clone())
        .with_color("steelblue")
        .with_legend("HMC");
    let mh_line = LinePlot::new()
        .with_data(result.momentum_scan_mh.clone())
        .with_color("mediumseagreen")
        .with_legend("MH");
    let expected_line = LinePlot::new()
        .with_data(result.momentum_scan_expected.clone())
        .with_color("coral")
        .with_legend("Expected");
    let scan_plots: Vec<Plot> = vec![hmc_line.into(), mh_line.into(), expected_line.into()];
    let scan_layout = Layout::auto_from_plots(&scan_plots)
        .with_title("Momentum Variance vs Mass (expected slope = 1/β)")
        .with_x_label("Mass")
        .with_y_label("Variance");
    let pdf = render_to_pdf(scan_plots, scan_layout).expect("render momentum scan");
    std::fs::write(&paths.momentum_scan, pdf)?;

    Ok(())
}

pub fn nearest_mass_index(result: &ExperimentResult, target_mass: f64) -> usize {
    result
        .momentum_scan_hmc
        .iter()
        .enumerate()
        .min_by(|(_, (a, _)), (_, (b, _))| {
            (a - target_mass)
                .abs()
                .partial_cmp(&(b - target_mass).abs())
                .unwrap()
        })
        .map(|(i, _)| i)
        .unwrap_or(0)
}

/// Bins `samples` into `num_bins` and returns `(bin_center, count)` pairs.
/// Useful when you'd rather plot one marker per bin than a stair-step.
pub fn histogram_centers(samples: &[f64], num_bins: usize) -> Vec<(f64, f64)> {
    if samples.is_empty() || num_bins == 0 {
        return Vec::new();
    }
    let lo = samples.iter().cloned().fold(f64::INFINITY, f64::min);
    let hi = samples.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if !(hi > lo) {
        return Vec::new();
    }
    let bin_w = (hi - lo) / num_bins as f64;
    let mut counts = vec![0usize; num_bins];
    for s in samples {
        let idx = (((s - lo) / bin_w).floor() as isize).clamp(0, num_bins as isize - 1) as usize;
        counts[idx] += 1;
    }
    counts
        .into_iter()
        .enumerate()
        .map(|(i, c)| (lo + (i as f64 + 0.5) * bin_w, c as f64))
        .collect()
}

/// Builds a 50-bin histogram in stair-step (x, y) form for line plotting.
pub fn build_histogram(samples: &[f64], num_bins: usize) -> Vec<(f64, f64)> {
    if samples.is_empty() || num_bins == 0 {
        return Vec::new();
    }
    let lo = samples.iter().cloned().fold(f64::INFINITY, f64::min);
    let hi = samples.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if !(hi > lo) {
        return Vec::new();
    }
    let bin_w = (hi - lo) / num_bins as f64;
    let mut counts = vec![0usize; num_bins];
    for s in samples {
        let idx = (((s - lo) / bin_w).floor() as isize).clamp(0, num_bins as isize - 1) as usize;
        counts[idx] += 1;
    }
    let mut out = Vec::with_capacity(num_bins * 4);
    for (i, &c) in counts.iter().enumerate() {
        let xs = lo + i as f64 * bin_w;
        let xe = lo + (i as f64 + 1.0) * bin_w;
        let y = c as f64;
        out.push((xs, 0.0));
        out.push((xs, y));
        out.push((xe, y));
        out.push((xe, 0.0));
    }
    out
}

/// Sample the analytic momentum distribution Normal(0, sqrt(m/β)) across
/// the same range the empirical samples occupy.
pub fn expected_gaussian_curve(mass: f64, beta: f64, samples: &[f64]) -> Vec<(f64, f64)> {
    if samples.is_empty() {
        return Vec::new();
    }
    let lo = samples.iter().cloned().fold(f64::INFINITY, f64::min);
    let hi = samples.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let pad = 0.05 * (hi - lo).abs().max(1.0);
    let lo = lo - pad;
    let hi = hi + pad;
    let sigma = (mass / beta).sqrt();
    let inv_norm = 1.0 / (sigma * (2.0 * std::f64::consts::PI).sqrt());
    (0..200)
        .map(|i| {
            let x = lo + (hi - lo) * i as f64 / 199.0;
            let z = x / sigma;
            (x, inv_norm * (-0.5 * z * z).exp())
        })
        .collect()
}

fn peak_height(data: &[(f64, f64)]) -> f64 {
    data.iter().map(|(_, y)| *y).fold(0.0_f64, f64::max)
}

/// Gaussian PDF (mean 0, σ from the population std of `samples`), sampled
/// over the x-range of `range_samples` and scaled by `scale`.
pub fn gaussian_fit_curve(
    samples: &[f64],
    range_samples: &[f64],
    scale: f64,
) -> Option<Vec<(f64, f64)>> {
    if samples.len() < 2 || range_samples.is_empty() {
        return None;
    }
    let (_mean, sigma) = crate::avg_and_std(samples);
    if !(sigma > 0.0) {
        return None;
    }
    let lo = range_samples.iter().cloned().fold(f64::INFINITY, f64::min);
    let hi = range_samples
        .iter()
        .cloned()
        .fold(f64::NEG_INFINITY, f64::max);
    let pad = 0.05 * (hi - lo).abs().max(1.0);
    let lo = lo - pad;
    let hi = hi + pad;
    let inv_norm = 1.0 / (sigma * (2.0 * std::f64::consts::PI).sqrt());
    Some(
        (0..200)
            .map(|i| {
                let x = lo + (hi - lo) * i as f64 / 199.0;
                let z = x / sigma;
                (x, inv_norm * (-0.5 * z * z).exp() * scale)
            })
            .collect(),
    )
}
