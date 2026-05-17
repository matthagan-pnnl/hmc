//! Drives a full HMC + MH + momentum-scan experiment from a single
//! `SimParams` value, returning everything the CLI and GUI need to plot.

use crate::metropolis_hastings::MetropolisHastings;
use crate::momentum_checker::MomentumChecker;
use crate::potential_sampler::PoissonHMC;
use crate::potentials::Potential;

/// Snapshot of every knob the experiment exposes. Mirrors the CLI struct
/// in `src/main.rs` but is Clone/PartialEq so the GUI can compare-and-rerun.
#[derive(Debug, Clone, PartialEq)]
pub struct SimParams {
    pub mass: f64,
    pub beta: f64,
    pub dimensions: usize,
    pub avg_sim_time: f64,
    pub step_size: f64,
    pub iterations: usize,
    pub num_chains: usize,
    pub acceptance_temp: Option<f64>,
    pub mass_upper_bound: f64,
    pub num_masses: usize,
    pub check_iterations: usize,
    pub proposal_std: f64,
    /// Lower edge of the sample window (percent of chain) fed to the
    /// momentum checker. Samples before this point are discarded.
    pub burn_in_lo_percent: f64,
    /// Upper edge of the sample window (percent of chain) fed to the
    /// momentum checker. Samples after this point are discarded.
    pub burn_in_hi_percent: f64,
    pub potential: Potential,
    /// Bounding box for the first position sample of every chain (HMC and MH).
    /// Broadcast to all `dimensions`.
    pub init_box: (f64, f64),
}

impl Default for SimParams {
    fn default() -> Self {
        Self {
            mass: 1.0,
            beta: 1.0,
            dimensions: 1,
            avg_sim_time: 1.0,
            step_size: 0.1,
            iterations: 1000,
            num_chains: 1,
            acceptance_temp: None,
            mass_upper_bound: 2.0,
            num_masses: 20,
            check_iterations: 500,
            proposal_std: 0.5,
            burn_in_lo_percent: 10.0,
            burn_in_hi_percent: 100.0,
            potential: Potential::default(),
            init_box: (-1.0, 1.0),
        }
    }
}

/// Coarse progress signal emitted while a job is running.
#[derive(Debug, Clone)]
pub enum Stage {
    Hmc { chain: usize, total_chains: usize, done: usize, total: usize },
    Mh { chain: usize, total_chains: usize, done: usize, total: usize },
    MomentumScanHmc,
    MomentumScanMh,
    Histogram,
    Done,
}

/// Bundle of everything the plotting layer needs.
pub struct ExperimentResult {
    pub params: SimParams,
    pub hmc_samples: Vec<Vec<f64>>,
    pub mh_samples: Vec<Vec<f64>>,
    pub samples_1d: Vec<f64>,
    pub momentum_scan_hmc: Vec<(f64, f64)>,
    pub momentum_scan_mh: Vec<(f64, f64)>,
    pub momentum_scan_expected: Vec<(f64, f64)>,
    /// Per-mass raw momentum samples from the HMC-positions momentum check,
    /// indexed identically to `momentum_scan_hmc`. Used for the
    /// per-mass momentum histogram in the GUI.
    pub hmc_momenta_per_mass: Vec<Vec<f64>>,
    /// Per-mass raw momentum samples from the MH-positions momentum check.
    pub mh_momenta_per_mass: Vec<Vec<f64>>,
    pub potential_curve: Vec<(f64, f64)>,
    pub histogram_bins: Vec<(f64, f64)>,
    pub mh_histogram_bins: Vec<(f64, f64)>,
    pub max_potential: f64,
}

/// Run the full experiment. `on_stage` is invoked at each major step
/// transition so callers can render progress.
pub fn run_experiment<P>(params: &SimParams, mut on_stage: P) -> ExperimentResult
where
    P: FnMut(Stage),
{
    let potential = params.potential;
    let potential_fn = move |x: &[f64]| potential.eval(x);
    let init_bounds: Vec<(f64, f64)> = (0..params.dimensions).map(|_| params.init_box).collect();

    let hmc = PoissonHMC::new(
        potential_fn,
        params.mass,
        params.beta,
        params.dimensions,
        params.avg_sim_time,
        params.step_size,
        params.acceptance_temp,
        Some(init_bounds.clone()),
    );

    let mut hmc_samples: Vec<Vec<f64>> = Vec::new();
    let mut samples_1d: Vec<f64> = Vec::new();
    for chain_id in 0..params.num_chains {
        let total_chains = params.num_chains;
        let iters = params.iterations;
        let chain = hmc.run_chain_with_progress(iters, |done, total| {
            on_stage(Stage::Hmc { chain: chain_id, total_chains, done, total });
        });
        samples_1d.extend(chain.iter().map(|x| x[0]));
        hmc_samples.extend(chain);
    }

    let mh = MetropolisHastings::new(
        move |x: &[f64]| potential.eval(x),
        params.beta,
        params.dimensions,
        params.proposal_std,
        Some(init_bounds),
    );
    let mut mh_samples: Vec<Vec<f64>> = Vec::new();
    for chain_id in 0..params.num_chains {
        let total_chains = params.num_chains;
        let iters = params.iterations;
        let chain = mh.run_chain_with_progress(iters, |done, total| {
            on_stage(Stage::Mh { chain: chain_id, total_chains, done, total });
        });
        mh_samples.extend(chain);
    }

    // Restrict to the [lo%, hi%] window of each chain before feeding into
    // the momentum checker. Bounds are clamped and ordered defensively.
    let hmc_post_burn = window_slice(
        &hmc_samples,
        params.burn_in_lo_percent,
        params.burn_in_hi_percent,
    );
    let mh_post_burn = window_slice(
        &mh_samples,
        params.burn_in_lo_percent,
        params.burn_in_hi_percent,
    );

    let checker = MomentumChecker::new(
        move |x: &[f64]| potential.eval(x),
        params.beta,
        params.mass,
        params.step_size,
        params.avg_sim_time,
        params.dimensions,
    );

    on_stage(Stage::MomentumScanHmc);
    let hmc_scan_full = checker.scan_masses_with_samples(
        params.mass_upper_bound,
        params.num_masses,
        hmc_post_burn,
        params.check_iterations,
    );
    let mut momentum_scan_hmc: Vec<(f64, f64)> = Vec::with_capacity(hmc_scan_full.len());
    let mut hmc_momenta_per_mass: Vec<Vec<f64>> = Vec::with_capacity(hmc_scan_full.len());
    for (m, std, momenta) in hmc_scan_full {
        momentum_scan_hmc.push((m, std));
        hmc_momenta_per_mass.push(momenta);
    }

    on_stage(Stage::MomentumScanMh);
    let mh_scan_full = checker.scan_masses_with_samples(
        params.mass_upper_bound,
        params.num_masses,
        mh_post_burn,
        params.check_iterations,
    );
    let mut momentum_scan_mh: Vec<(f64, f64)> = Vec::with_capacity(mh_scan_full.len());
    let mut mh_momenta_per_mass: Vec<Vec<f64>> = Vec::with_capacity(mh_scan_full.len());
    for (m, var, momenta) in mh_scan_full {
        momentum_scan_mh.push((m, var));
        mh_momenta_per_mass.push(momenta);
    }

    // Expected variance of momentum at equilibrium: m/β (slope 1/β through origin).
    let momentum_scan_expected: Vec<(f64, f64)> = (1..=params.num_masses)
        .map(|i| {
            let mass = i as f64 * params.mass_upper_bound / params.num_masses as f64;
            (mass, mass / params.beta)
        })
        .collect();

    on_stage(Stage::Histogram);

    let (plot_lo, plot_hi) = params.potential.default_range();
    let xs: Vec<f64> = (0..200)
        .map(|ix| plot_lo + (plot_hi - plot_lo) * ix as f64 / 199.0)
        .collect();
    let ys: Vec<f64> = xs.iter().map(|x| potential.eval(&[*x])).collect();
    let max_potential = ys.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let potential_curve: Vec<(f64, f64)> =
        xs.iter().cloned().zip(ys.iter().cloned()).collect();

    let histogram_bins = position_histogram_bins(&samples_1d, 50, max_potential);
    let mh_samples_1d: Vec<f64> = mh_samples.iter().map(|x| x[0]).collect();
    let mh_histogram_bins = position_histogram_bins(&mh_samples_1d, 50, max_potential);

    on_stage(Stage::Done);

    ExperimentResult {
        params: params.clone(),
        hmc_samples,
        mh_samples,
        samples_1d,
        momentum_scan_hmc,
        momentum_scan_mh,
        momentum_scan_expected,
        hmc_momenta_per_mass,
        mh_momenta_per_mass,
        potential_curve,
        histogram_bins,
        mh_histogram_bins,
        max_potential,
    }
}

/// Slice `chain[lo%..hi%]`, with bounds clamped and ordered so misuse can't
/// panic. Returns an empty slice if the resulting window is degenerate.
fn window_slice(chain: &[Vec<f64>], lo_pct: f64, hi_pct: f64) -> &[Vec<f64>] {
    let n = chain.len();
    if n == 0 {
        return &[];
    }
    let (lo_pct, hi_pct) = if lo_pct <= hi_pct {
        (lo_pct, hi_pct)
    } else {
        (hi_pct, lo_pct)
    };
    let to_idx =
        |pct: f64| (((pct.clamp(0.0, 100.0) / 100.0) * n as f64).round() as usize).min(n);
    let lo = to_idx(lo_pct);
    let hi = to_idx(hi_pct).max(lo);
    &chain[lo..hi]
}

/// 50-bin (by default) stair-step histogram, with bar heights scaled so the
/// tallest bar matches `target_peak` — keeps the bars visually comparable
/// to the potential overlay.
fn position_histogram_bins(samples: &[f64], num_bins: usize, target_peak: f64) -> Vec<(f64, f64)> {
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
    let max_count = *counts.iter().max().unwrap_or(&1) as f64;
    let scale = if max_count > 0.0 {
        target_peak / max_count
    } else {
        1.0
    };
    let mut out = Vec::with_capacity(num_bins * 4);
    for (i, &c) in counts.iter().enumerate() {
        let xs = lo + i as f64 * bin_w;
        let xe = lo + (i as f64 + 1.0) * bin_w;
        let y = c as f64 * scale;
        out.push((xs, 0.0));
        out.push((xs, y));
        out.push((xe, y));
        out.push((xe, 0.0));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn run_experiment_default_params_produces_finite_samples() {
        let mut params = SimParams::default();
        params.iterations = 50;
        params.check_iterations = 20;
        params.num_masses = 5;
        let result = run_experiment(&params, |_| {});
        assert!(!result.hmc_samples.is_empty());
        assert!(!result.mh_samples.is_empty());
        assert_eq!(result.momentum_scan_hmc.len(), params.num_masses);
        for (m, s) in &result.momentum_scan_hmc {
            assert!(m.is_finite());
            assert!(s.is_finite());
            assert!(*s >= 0.0);
        }
    }
}
