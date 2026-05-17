use rand::{Rng, RngExt};
use rand_distr::{Normal, Poisson, Uniform};

const STEP_SIZE: f64 = 1e-7;

pub mod metropolis_hastings;
pub mod momentum_checker;
pub mod potential_sampler;
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

/// Returns `(avg, std)`, computed using kahan summation and the standard deviation
/// as estimated from the population variance.
pub fn avg_and_std(data: &[f64]) -> (f64, f64) {
    let avg = kahan_summation(data) / data.len() as f64;
    let diffs: Vec<f64> = data.iter().map(|x| (*x - avg).powi(2)).collect();
    let var = kahan_summation(&diffs) / (data.len() - 1) as f64;
    (avg, var.sqrt())
}
/// Determines (slope, intercept, RMS error) of data `(x, y)`
pub fn linear_regression_basic(data: &[(f64, f64)]) -> (f64, f64, f64) {
    if data.len() < 2 {
        panic!("Cannot least squares with less than 2 points. You should feel bad.");
    }
    let x_data: Vec<f64> = data.iter().map(|(x, _y)| *x).collect();
    let y_data: Vec<f64> = data.iter().map(|(_x, y)| *y).collect();
    let x_mean = kahan_summation(&x_data) / data.len() as f64;
    let y_mean = kahan_summation(&y_data) / data.len() as f64;
    let numerator_points: Vec<f64> = x_data
        .iter()
        .zip(y_data.iter())
        .map(|(x, y)| (*x - x_mean) * (*y - y_mean))
        .collect();
    let denominator_points: Vec<f64> = x_data
        .iter()
        .map(|x| (*x - x_mean) * (*x - x_mean))
        .collect();
    let slope = kahan_summation(&numerator_points) / kahan_summation(&denominator_points);
    let intercept = y_mean - slope * x_mean;
    let sum_of_errors = x_data
        .iter()
        .zip(y_data.iter())
        .map(|(x, y)| (y - (*x * slope + intercept)).powi(2))
        .sum::<f64>();
    let rms_error = (sum_of_errors / data.len() as f64).sqrt();
    (slope, intercept, rms_error)
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
