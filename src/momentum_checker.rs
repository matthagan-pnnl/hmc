use rand::RngExt;
use rand_distr::{Normal, Poisson, Uniform};

use crate::leapfrog;

pub struct MomentumChecker<F>
where
    F: Fn(&[f64]) -> f64,
{
    potential_energy: F,
    beta: f64,
    mass: f64,
    epsilon: f64,
    average_simulation_time: f64,
    step_size_sampler: Poisson<f64>,
    dimensions: usize,
}

impl<F> MomentumChecker<F>
where
    F: Fn(&[f64]) -> f64,
{
    pub fn new(
        potential_energy: F,
        beta: f64,
        mass: f64,
        epsilon: f64,
        average_simulation_time: f64,
        dimensions: usize,
    ) -> Self {
        let poisson_average = average_simulation_time / epsilon;
        Self {
            potential_energy,
            beta,
            mass,
            epsilon,
            average_simulation_time,
            step_size_sampler: Poisson::new(poisson_average).unwrap(),
            dimensions,
        }
    }

    /// Runs the momentum variance check against `position_samples`.
    /// Returns per-dimension sample variances of the evolved momenta.
    /// Compare each element to `self.expected_variance()`.
    pub fn run(&self, position_samples: &[Vec<f64>], n_iterations: usize) -> Vec<f64> {
        let mut rng = rand::rng();
        let momentum_dist = Normal::new(0.0, (self.mass / self.beta).sqrt()).unwrap();
        let index_dist = Uniform::new(0, position_samples.len()).unwrap();

        let mut p: Vec<f64> = (0..self.dimensions)
            .map(|_| rng.sample(momentum_dist))
            .collect();
        let mut recorded_momenta: Vec<Vec<f64>> = Vec::with_capacity(n_iterations);

        for _ in 0..n_iterations {
            let n_steps = (rng.sample(self.step_size_sampler).round() as usize).max(1);
            let mut x = position_samples[rng.sample(index_dist)].clone();
            leapfrog(
                &mut x,
                &mut p,
                self.mass,
                n_steps,
                self.epsilon,
                &self.potential_energy,
            );
            recorded_momenta.push(p.clone());
        }

        (0..self.dimensions)
            .map(|dim| {
                let values: Vec<f64> = recorded_momenta.iter().map(|p| p[dim]).collect();
                let (_mean, std) = crate::avg_and_std(&values);
                std * std
            })
            .collect()
    }

    /// Sweeps over `num_masses` evenly-spaced mass values up to `mass_upper_bound`,
    /// running the momentum chain for each, and returns `(mass, std)` pairs.
    pub fn scan_masses(
        &self,
        mass_upper_bound: f64,
        num_masses: usize,
        position_samples: &[Vec<f64>],
        n_iterations: usize,
    ) -> Vec<(f64, f64)> {
        let index_dist = Uniform::new(0, position_samples.len()).unwrap();
        (1..=num_masses)
            .map(|i| {
                let mass = i as f64 * mass_upper_bound / num_masses as f64;
                // Scale epsilon linearly with mass so that epsilon/mass stays constant,
                // keeping the leapfrog position-update rate stable across the mass scan.
                let epsilon = self.epsilon * mass / self.mass;
                let step_sampler =
                    Poisson::new(self.average_simulation_time / epsilon).unwrap();
                let momentum_dist = Normal::new(0.0, (mass / self.beta).sqrt()).unwrap();
                let mut rng = rand::rng();
                let mut p: Vec<f64> = (0..self.dimensions)
                    .map(|_| rng.sample(momentum_dist))
                    .collect();
                let mut recorded_momenta: Vec<Vec<f64>> = Vec::with_capacity(n_iterations);
                for _ in 0..n_iterations {
                    let n_steps = (rng.sample(step_sampler).round() as usize).max(1);
                    let mut x = position_samples[rng.sample(index_dist)].clone();
                    leapfrog(&mut x, &mut p, mass, n_steps, epsilon, &self.potential_energy);
                    recorded_momenta.push(p.clone());
                }
                let mean_var = (0..self.dimensions)
                    .map(|dim| {
                        let values: Vec<f64> = recorded_momenta.iter().map(|p| p[dim]).collect();
                        let (_mean, std) = crate::avg_and_std(&values);
                        std * std
                    })
                    .sum::<f64>()
                    / self.dimensions as f64;
                (mass, mean_var.sqrt())
            })
            .collect()
    }

    /// Expected per-dimension momentum variance under equilibrium: mass / beta.
    pub fn expected_variance(&self) -> f64 {
        self.mass / self.beta
    }
}
