use rand::{Rng, RngExt};
use rand_distr::{Normal, Uniform};

pub struct MetropolisHastings<F>
where
    F: Fn(&[f64]) -> f64,
{
    potential_energy: F,
    acceptance_temperature: f64,
    proposal_std: f64,
    dimensions: usize,
    position_bounds: Vec<(f64, f64)>,
}

impl<F> MetropolisHastings<F>
where
    F: Fn(&[f64]) -> f64,
{
    pub fn new(
        potential_energy: F,
        acceptance_temperature: f64,
        dimensions: usize,
        proposal_std: f64,
        position_bounds: Option<Vec<(f64, f64)>>,
    ) -> Self {
        let position_bounds =
            position_bounds.unwrap_or((0..dimensions).map(|_| (-1.0, 1.0)).collect());
        Self {
            potential_energy,
            acceptance_temperature,
            proposal_std,
            dimensions,
            position_bounds,
        }
    }

    pub fn apply_transition<R: Rng>(&self, x: &mut [f64], rng: &mut R) {
        let proposal_dist = Normal::new(0.0, self.proposal_std).unwrap();
        let mut x_proposed = x.to_vec();
        for ix in 0..self.dimensions {
            x_proposed[ix] += rng.sample(proposal_dist);
        }
        let delta = self.acceptance_temperature * ((self.potential_energy)(&x_proposed) - (self.potential_energy)(x));
        let acceptance_probability = 1.0 / (1.0 + delta.exp());
        if rng.sample(Uniform::new(0.0, 1.0).unwrap()) < acceptance_probability {
            x.copy_from_slice(&x_proposed);
        }
    }

    pub fn run_chain(&self, num_iterations: usize) -> Vec<Vec<f64>> {
        self.run_chain_with_progress(num_iterations, |_, _| {})
    }

    /// Like `run_chain`, but invokes `on_progress(completed, total)` periodically.
    pub fn run_chain_with_progress<P>(
        &self,
        num_iterations: usize,
        mut on_progress: P,
    ) -> Vec<Vec<f64>>
    where
        P: FnMut(usize, usize),
    {
        let mut rng = rand::rng();
        let mut x: Vec<f64> = self
            .position_bounds
            .iter()
            .map(|(lo, hi)| rng.sample(Uniform::new_inclusive(lo, hi).unwrap()))
            .collect();
        let mut chain_history = vec![x.clone()];
        let progress_interval = (num_iterations / 20).max(1);
        for iteration in 0..num_iterations {
            self.apply_transition(&mut x, &mut rng);
            chain_history.push(x.clone());
            if (iteration + 1) % progress_interval == 0 {
                on_progress(iteration + 1, num_iterations);
            }
        }
        chain_history
    }
}
