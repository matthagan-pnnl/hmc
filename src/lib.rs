use rand::{Rng, RngExt};
use rand_distr::{Normal, Poisson, Uniform};

const STEP_SIZE: f64 = 1e-7;

pub mod sampler;

pub fn kahan_summation(data: &[f64]) -> f64 {
    let mut tot = 0.0;
    let mut carry = 0.0;
    for x in data.iter() {
        let y = x - carry;
        let t = tot + y;
        carry = (t - tot) - y;
        tot = t;
    }
    tot
}

pub fn leapfrog<F>(
    x: &mut [f64],
    p: &mut [f64],
    mass: f64,
    n_steps: usize,
    delta_t: f64,
    potential_energy: &F,
) where
    F: Fn(&[f64]) -> f64,
{
    assert!(n_steps > 0);
    let grad = finite_difference(&x, potential_energy);

    for ix in 0..p.len() {
        p[ix] -= 0.5 * delta_t * grad[ix];
    }

    for n_step in 0..n_steps {
        for ix in 0..x.len() {
            x[ix] += (delta_t / mass) * p[ix];
        }

        let grad = finite_difference(&x, potential_energy);
        if n_step == n_steps - 1 {
            for ix in 0..p.len() {
                p[ix] -= 0.5 * delta_t * grad[ix];
            }
        } else {
            for ix in 0..p.len() {
                p[ix] -= delta_t * grad[ix];
            }
        }
    }
}

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

        let momentum_distribution = Normal::new(0.0, self.mass / self.beta.sqrt()).unwrap();
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

pub fn finite_difference<F>(x: &[f64], f: F) -> Vec<f64>
where
    F: Fn(&[f64]) -> f64,
{
    let mut gradient = vec![0.0; x.len()];
    let mut x_prime = x.to_vec();
    for dim_ix in 0..x.len() {
        let original = x_prime[dim_ix];
        x_prime[dim_ix] = original + 0.5 * STEP_SIZE;
        let f_plus = f(&x_prime);
        x_prime[dim_ix] = original - 0.5 * STEP_SIZE;
        let f_minus = f(&x_prime);
        x_prime[dim_ix] = original;
        let numerator = f_plus - f_minus;
        gradient[dim_ix] = numerator / STEP_SIZE;
    }
    gradient
}

#[cfg(test)]
mod tests {
    use crate::finite_difference;

    #[test]
    fn test_finite_difference_power() {
        let f = |x: &[f64]| x.iter().map(|x_i| x_i.powi(7)).sum::<f64>();
        let x = [1.0];
        let grad = finite_difference(&x, f);
        let expected = 7.0;
        assert!(
            (grad[0] - expected).abs() < 1e-6,
            "Expected {}, got {}",
            expected,
            grad[0]
        );
    }

    #[test]
    fn test_finite_difference_multivariate() {
        // f(x, y) = x^2 + y^3, grad f = [2x, 3y^2]
        let f = |x: &[f64]| x[0].powi(2) + x[1].powi(3);
        let x = [2.0, 3.0];
        let grad = finite_difference(&x, f);
        let expected = [4.0, 27.0];

        for i in 0..x.len() {
            assert!(
                (grad[i] - expected[i]).abs() < 1e-6,
                "At index {}, expected {}, got {}",
                i,
                expected[i],
                grad[i]
            );
        }
    }
}
