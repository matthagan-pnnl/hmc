use clap::Parser;
use hmc::export::{write_pdfs, ExportPaths};
use hmc::potentials::Potential;
use hmc::run::{run_experiment, SimParams, Stage};

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

    /// Upper bound for mass scan in momentum checker
    #[arg(long, default_value = "2.0")]
    mass_upper_bound: f64,

    /// Number of masses to scan in momentum checker
    #[arg(long, default_value = "20")]
    num_masses: usize,

    /// Iterations per mass for momentum checker
    #[arg(long, default_value = "500")]
    check_iterations: usize,

    /// Proposal std for Metropolis-Hastings
    #[arg(long, default_value = "0.5")]
    proposal_std: f64,

    /// Lower edge (percent) of the sample window fed to the momentum checker
    #[arg(long, default_value = "10")]
    burn_in_lo_percent: f64,

    /// Upper edge (percent) of the sample window fed to the momentum checker
    #[arg(long, default_value = "100")]
    burn_in_hi_percent: f64,

    /// Left well center for the default double-well potential
    #[arg(long, default_value = "-1.0")]
    x_left: f64,

    /// Right well center for the default double-well potential
    #[arg(long, default_value = "1.0")]
    x_right: f64,

    /// Lower edge of the initial-position bounding box (broadcast over dims)
    #[arg(long, default_value = "-1.0")]
    init_box_lo: f64,

    /// Upper edge of the initial-position bounding box (broadcast over dims)
    #[arg(long, default_value = "1.0")]
    init_box_hi: f64,
}

impl From<&Cli> for SimParams {
    fn from(cli: &Cli) -> Self {
        SimParams {
            mass: cli.mass,
            beta: cli.beta,
            dimensions: cli.dimensions,
            avg_sim_time: cli.avg_sim_time,
            step_size: cli.step_size,
            iterations: cli.iterations,
            num_chains: cli.num_chains,
            acceptance_temp: cli.acceptance_temp,
            mass_upper_bound: cli.mass_upper_bound,
            num_masses: cli.num_masses,
            check_iterations: cli.check_iterations,
            proposal_std: cli.proposal_std,
            burn_in_lo_percent: cli.burn_in_lo_percent,
            burn_in_hi_percent: cli.burn_in_hi_percent,
            potential: Potential::DoubleWell {
                x_left: cli.x_left,
                x_right: cli.x_right,
            },
            init_box: (cli.init_box_lo, cli.init_box_hi),
        }
    }
}

fn main() {
    let cli = Cli::parse();
    let params: SimParams = (&cli).into();

    println!(
        "Running {} HMC + MH chain pair(s) with {} iterations each...",
        params.num_chains, params.iterations
    );

    let result = run_experiment(&params, |stage| match stage {
        Stage::Hmc { chain, total_chains, done, total } => {
            let pct = (done as f64 / total as f64) * 100.0;
            println!("HMC chain {}/{}: {:.0}%", chain + 1, total_chains, pct);
        }
        Stage::Mh { chain, total_chains, done, total } => {
            let pct = (done as f64 / total as f64) * 100.0;
            println!("MH chain {}/{}: {:.0}%", chain + 1, total_chains, pct);
        }
        Stage::MomentumScanHmc => println!("Scanning masses (HMC)..."),
        Stage::MomentumScanMh => println!("Scanning masses (MH)..."),
        Stage::Histogram => println!("Building histograms..."),
        Stage::Done => println!("Simulation complete; rendering PDFs..."),
    });

    let paths = ExportPaths::default();
    write_pdfs(&result, &paths).expect("write PDFs");
    println!("Wrote {}", paths.position_histogram);
    println!("Wrote {}", paths.mh_position_histogram);
    println!("Wrote {}", paths.momentum_histogram);
    println!("Wrote {}", paths.momentum_scan);
}
