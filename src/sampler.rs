use rand::RngExt;

#[derive(Debug, Clone)]
pub struct AliasSampler {
    og_table: Vec<f64>,
    alias_table: Vec<usize>,
}

impl AliasSampler {
    pub fn new(weights: &[f64]) -> Self {
        let n = weights.len();
        let norm = crate::kahan_summation(weights);
        let mut og_table = vec![0.0f64; n];
        let mut alias_table = vec![0usize; n];
        let mut small: Vec<usize> = Vec::new();
        let mut large: Vec<usize> = Vec::new();

        for (ix, weight) in weights.iter().enumerate() {
            let q = *weight / norm * n as f64;
            og_table[ix] = q;
            if q < 1.0 {
                small.push(ix);
            } else {
                large.push(ix);
            }
        }

        while !small.is_empty() && !large.is_empty() {
            let l = small.pop().unwrap();
            let g = *large.last().unwrap();
            alias_table[l] = g;
            og_table[g] -= 1.0 - og_table[l];
            if og_table[g] < 1.0 {
                large.pop();
                small.push(g);
            }
        }
        for &idx in small.iter().chain(large.iter()) {
            og_table[idx] = 1.0;
        }

        Self {
            og_table,
            alias_table,
        }
    }

    /// Draw one sample. Returns the index of the sampled weight.
    pub fn sample<R: rand::Rng>(&self, rng: &mut R) -> usize {
        let i = rng.random_range(0..self.og_table.len());
        if rng.random::<f64>() < self.og_table[i] {
            i
        } else {
            self.alias_table[i]
        }
    }
}
