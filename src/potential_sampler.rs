use rand::{Rng, RngExt};
use rand_distr::{Normal, Poisson, Uniform};

use crate::leapfrog;

pub struct PoissonHMC<F>
where
    F: Fn(&[f64]) -> f64,
{
    potential_energy: F,
    mass: f64,
    beta: f64,
    epsilon: f64,
    step_size_sampler: Poisson<f64>,
    acceptance_temperature: f64,
    position_bounds: Vec<(f64, f64)>,
    dimensions: usize,
}

impl<F> PoissonHMC<F>
where
    F: Fn(&[f64]) -> f64,
{
    /// sets all mass dimensions equal for now
    /// Resorts to an acceptance temperature of 1.0.
    pub fn new(
        potential_energy: F,
        mass: f64,
        beta: f64,
        dimensions: usize,
        average_simulation_time: f64,
        step_size: f64,
        acceptance_temperature: Option<f64>,
        position_bounds: Option<Vec<(f64, f64)>>,
    ) -> Self {
        let poisson_average = average_simulation_time / step_size;
        let acceptance_temperature = acceptance_temperature.unwrap_or(1.0);
        let position_bounds =
            position_bounds.unwrap_or((0..dimensions).map(|_| (-1.0, 1.0)).collect());
        Self {
            potential_energy,
            mass,
            beta,
            dimensions,
            epsilon: step_size,
            step_size_sampler: Poisson::new(poisson_average).unwrap(),
            acceptance_temperature,
            position_bounds,
        }
    }

    pub fn apply_transition<R: Rng>(&self, x: &mut [f64], p: &mut [f64], rng: &mut R) {
        let num_steps = rng.sample(self.step_size_sampler).round() as usize;
        let num_steps = num_steps.max(1);
        let energy_before = (self.potential_energy)(x)
            + p.iter()
                .map(|p_i| p_i * p_i / (2.0 * self.mass))
                .sum::<f64>();
        let initial_x = x.to_vec();
        let initial_p = p.to_vec();
        leapfrog(
            x,
            p,
            self.mass,
            num_steps,
            self.epsilon,
            &self.potential_energy,
        );
        let energy_after = (self.potential_energy)(x)
            + p.iter()
                .map(|p_i| p_i * p_i / (2.0 * self.mass))
                .sum::<f64>();
        let delta = self.acceptance_temperature * (energy_after - energy_before);
        let acceptance_probability = 1.0 / (1.0 + delta.exp());
        let sample = rng.sample(Uniform::new(0.0, 1.0).unwrap());
        if sample <= acceptance_probability {
            for ix in 0..x.len() {
                x[ix] = initial_x[ix];
                p[ix] = initial_p[ix];
            }
        }
    }

    pub fn run_chain(&self, num_iterations: usize) -> Vec<Vec<f64>> {
        let mut x = Vec::new();
        let mut p = Vec::new();
        let mut rng = rand::rng();
        for (lower, upper) in self.position_bounds.iter() {
            let sample = rng.sample(Uniform::new_inclusive(lower, upper).unwrap());
            x.push(sample);
        }

        let momentum_distribution = Normal::new(0.0, (self.mass / self.beta).sqrt()).unwrap();
        for _ in 0..self.dimensions {
            p.push(rng.sample(momentum_distribution));
        }
        let mut chain_history = Vec::new();
        chain_history.push(x.clone());
        let progress_interval = (num_iterations / 20).max(1);
        for iteration in 0..num_iterations {
            self.apply_transition(&mut x, &mut p, &mut rng);
            chain_history.push(x.to_vec());
            p.clear();
            for _ in 0..self.dimensions {
                p.push(rng.sample(momentum_distribution));
            }
            if (iteration + 1) % progress_interval == 0 {
                let percentage = ((iteration + 1) as f64 / num_iterations as f64) * 100.0;
                println!("Progress: {:.0}%", percentage);
            }
        }
        chain_history
    }
}
