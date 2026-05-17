/// A small set of potential-energy shapes the experiments can target.
///
/// `eval` is dimension-agnostic. The shipped variants only read `x[0]`,
/// matching the historical CLI behavior — extending to higher dims is a
/// matter of adding new variants.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Potential {
    /// Symmetric quartic double well: `(x - x_left)^2 * (x - x_right)^2`.
    DoubleWell { x_left: f64, x_right: f64 },
    /// Harmonic well: `0.5 * k * x^2`.
    Harmonic { k: f64 },
}

impl Potential {
    pub fn eval(&self, x: &[f64]) -> f64 {
        match *self {
            Potential::DoubleWell { x_left, x_right } => {
                let a = x[0] - x_left;
                let b = x[0] - x_right;
                a * a * b * b
            }
            Potential::Harmonic { k } => 0.5 * k * x[0] * x[0],
        }
    }

    /// A reasonable x-range for plotting / sampling initialization.
    pub fn default_range(&self) -> (f64, f64) {
        match *self {
            Potential::DoubleWell { x_left, x_right } => {
                let pad = 0.5 * (x_right - x_left).abs().max(1.0);
                (x_left - pad, x_right + pad)
            }
            Potential::Harmonic { .. } => (-3.0, 3.0),
        }
    }
}

impl Default for Potential {
    fn default() -> Self {
        Potential::DoubleWell {
            x_left: -1.0,
            x_right: 1.0,
        }
    }
}
