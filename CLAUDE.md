# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

Numerical investigation of a momentum-equipartition convergence diagnostic for
MCMC samplers, described in `theory.typ` / `theory.pdf`. Given a Markov chain
targeting `pi(x) ∝ exp(-beta V(x))`, the diagnostic draws positions from the
chain, attaches fresh canonical momenta, evolves under Hamiltonian dynamics,
and checks whether the resulting momentum variance matches the equipartition
prediction `m / beta`. Persistent deviation is evidence the chain has not
converged.

The reference experiment compares a Poisson-step HMC sampler (`PoissonHMC`)
against a Metropolis–Hastings baseline on a 1-D double-well potential, sweeps
particle mass, and renders three PDFs.

## Common commands

- `cargo run --release -- [flags]` — run the experiment. Release matters;
  debug builds are dramatically slower for the chain loops.
- `cargo test` — runs the `finite_difference` tests in `src/lib.rs`.
- `cargo build --release` — build only.
- Full flag list lives on the `Cli` struct in `src/main.rs` (mass, beta,
  dimensions, step_size, avg_sim_time, iterations, num_chains,
  acceptance_temp, mass_upper_bound, num_masses, check_iterations,
  proposal_std, burn_in_percent). Outputs land in the repo root as
  `hmc_chain.pdf`, `hmc_histogram.pdf`, `momentum_std_vs_mass.pdf`.
- `typst compile theory.typ` — rebuild the theory report (needs the `typst`
  CLI; not part of the cargo build).

## Architecture

- `src/lib.rs` — shared numerical primitives: `kahan_summation`,
  `avg_and_std`, `linear_regression_basic`, `finite_difference` (central
  differences with hardcoded `STEP_SIZE = 1e-7`), and the `leapfrog`
  integrator (half-kick / drift / half-kick, with gradient recomputed at
  every drift). All samplers depend on these.
- `src/potential_sampler.rs` — `PoissonHMC`. The number of leapfrog steps
  per proposal is `Poisson(avg_sim_time / epsilon)`. Acceptance uses a
  Fermi-Dirac form `p_accept = 1 / (1 + exp(temp · ΔH))`, not the usual
  `min(1, exp(-ΔH))`. **Note the inversion in `apply_transition`**: when
  the uniform draw is `<= p_accept`, the code restores the *initial*
  `(x, p)` — that branch is the reject branch. Be careful when editing.
- `src/metropolis_hastings.rs` — Gaussian random-walk MH with the same
  Fermi-Dirac acceptance form. Baseline for the diagnostic comparison.
- `src/momentum_checker.rs` — the diagnostic itself. `run` evaluates
  momentum variance at the current mass; `scan_masses` sweeps masses and
  **scales `epsilon` linearly with mass** (`epsilon * mass / self.mass`)
  so the leapfrog drift rate `epsilon / mass` stays constant across the
  scan. The Poisson step sampler is re-derived per mass for the same
  reason.
- `src/sampler.rs` — generic Vose/Walker `AliasSampler` (`O(1)` draws
  from a discrete distribution). Currently unused by `main.rs`; keep it.
- `src/main.rs` — wiring: parses CLI, defines the double-well potential
  `(x - x_left)^2 (x - x_right)^2` inline, runs HMC + MH chains, plots
  momentum-std-vs-mass (after dropping `burn_in_percent` of samples),
  and writes the three PDFs via `kuva`.

## Conventions worth knowing

- The potential is always `Fn(&[f64]) -> f64`, so multi-dimensional
  extensions come for free. The current hardcoded potential only reads
  `x[0]` despite supporting `--dimensions`.
- `beta` is the physical inverse temperature. `acceptance_temp` is a
  *separate* knob (a tempering parameter on the Fermi-Dirac acceptance).
  Don't conflate them.
- `STEP_SIZE = 1e-7` in `src/lib.rs` is the finite-difference step for
  the gradient, **not** the leapfrog step `epsilon`. Easy to misread.
- `kuva` is the plotting crate; PDFs are produced via
  `render_to_pdf(plots, layout)` and written with `std::fs::write`.
- Rust edition is 2024.

## Reference material

- `theory.typ` (+ rendered `theory.pdf`) — derivation and motivation for
  the momentum-checker diagnostic. Read this before changing
  `MomentumChecker` semantics or before writing about *why* the code
  does what it does.
