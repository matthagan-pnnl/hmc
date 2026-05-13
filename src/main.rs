use clap::Parser;
use kuva::prelude::*;
use hmc::PoissonHMC;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Cli {
    /// Mass parameter for HMC
    #[arg(short, long, default_value = "1.0")]
    mass: f64,

    /// Inverse temperature (beta)
    #[arg(short, long, default_value = "1.0")]
    beta: f64,

    /// Number of dimensions
    #[arg(short, long, default_value = "1")]
    dimensions: usize,

    /// Average simulation time for Poisson step sampler
    #[arg(long, default_value = "1.0")]
    avg_sim_time: f64,

    /// Step size for leapfrog integration
    #[arg(short, long, default_value = "0.1")]
    step_size: f64,

    /// Number of chain iterations
    #[arg(short, long, default_value = "1000")]
    iterations: usize,

    /// Number of independent chains to run
    #[arg(long, default_value = "1")]
    num_chains: usize,

    /// Acceptance temperature (optional)
    #[arg(long)]
    acceptance_temp: Option<f64>,
}

fn main() {
    let cli = Cli::parse();

    let x_left = -1.0;
    let x_right = 1.0;
    let potential = |x: &[f64]| (x[0] - x_left) * (x[0] - x_left) * (x[0] - x_right) * (x[0] - x_right);

    let sampler = PoissonHMC::new(
        potential,
        cli.mass,
        cli.beta,
        cli.dimensions,
        cli.avg_sim_time,
        cli.step_size,
        cli.acceptance_temp,
        None,
    );

    println!("Running {} HMC chain(s) with {} iterations each...", cli.num_chains, cli.iterations);
    let mut all_samples: Vec<Vec<f64>> = Vec::new();
    let mut all_samples_1d: Vec<f64> = Vec::new();

    for chain_id in 0..cli.num_chains {
        println!("Starting chain {}/{}", chain_id + 1, cli.num_chains);
        let chain = sampler.run_chain(cli.iterations);
        let samples_1d: Vec<f64> = chain.iter().map(|x| x[0]).collect();
        all_samples_1d.extend(samples_1d);
        all_samples.extend(chain);
    }

    let samples_1d = all_samples_1d;

    // Create scatter plot of chain samples with potential energy
    let sample_data: Vec<(f64, f64)> = all_samples
        .iter()
        .map(|x| (x[0], potential(x)))
        .collect();

    let samples_plot = ScatterPlot::new()
        .with_data(sample_data)
        .with_color("steelblue")
        .with_size(3.0);

    // Create plot of potential function
    let xs = (0..200)
        .map(|ix| -2.0 + 0.02 * ix as f64)
        .collect::<Vec<f64>>();
    let ys = xs.iter().map(|x| potential(&[*x])).collect::<Vec<f64>>();
    let max_potential = ys.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let potential_data: Vec<(f64, f64)> = xs.into_iter().zip(ys.into_iter()).collect();

    let potential_plot = LinePlot::new()
        .with_data(potential_data.clone())
        .with_color("coral");

    // Create normalized histogram plot
    let num_bins = 50;
    let min_sample = samples_1d.iter().cloned().fold(f64::INFINITY, f64::min);
    let max_sample = samples_1d.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let bin_width = (max_sample - min_sample) / num_bins as f64;

    let mut bin_counts = vec![0usize; num_bins];
    for &sample in &samples_1d {
        let bin_idx = ((sample - min_sample) / bin_width).floor() as usize;
        if bin_idx < num_bins {
            bin_counts[bin_idx] += 1;
        }
    }

    let max_count = *bin_counts.iter().max().unwrap_or(&1) as f64;
    let scale_factor = max_potential / max_count;

    let mut histogram_data: Vec<(f64, f64)> = Vec::new();
    for (i, &count) in bin_counts.iter().enumerate() {
        let x_start = min_sample + i as f64 * bin_width;
        let x_end = min_sample + (i as f64 + 1.0) * bin_width;
        let y = count as f64 * scale_factor;
        histogram_data.push((x_start, 0.0));
        histogram_data.push((x_start, y));
        histogram_data.push((x_end, y));
        histogram_data.push((x_end, 0.0));
    }

    let histogram_line = LinePlot::new()
        .with_data(histogram_data)
        .with_color("steelblue");

    let potential_line = LinePlot::new()
        .with_data(potential_data)
        .with_color("coral");

    let plots: Vec<Plot> = vec![
        samples_plot.into(),
        potential_plot.into(),
    ];

    let histogram_plots: Vec<Plot> = vec![
        histogram_line.into(),
        potential_line.into(),
    ];

    let layout = Layout::auto_from_plots(&plots)
        .with_title("HMC Chain Samples and Potential Function")
        .with_x_label("Position")
        .with_y_label("Potential Energy");

    let histogram_layout = Layout::auto_from_plots(&histogram_plots)
        .with_title("Histogram of Chain Samples with Potential Function")
        .with_x_label("Position")
        .with_y_label("Count / Potential Energy");

    println!("Rendering plots...");
    let pdf = render_to_pdf(plots, layout).unwrap();
    std::fs::write("hmc_chain.pdf", pdf).unwrap();
    println!("Plot saved to hmc_chain.pdf");

    let histogram_pdf = render_to_pdf(histogram_plots, histogram_layout).unwrap();
    std::fs::write("hmc_histogram.pdf", histogram_pdf).unwrap();
    println!("Histogram saved to hmc_histogram.pdf");
}
