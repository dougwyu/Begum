mod dna;
mod filter;
mod sort;

use anyhow::Result;
use clap::{Parser, Subcommand};

#[derive(Parser)]
#[command(
    name = "begum",
    version = "0.1",
    about = "Metabarcoding and eDNA sequence preprocessing tool"
)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    Sort(SortArgs),
    Filter(FilterArgs),
}

#[derive(Parser, Debug)]
pub struct SortArgs {
    #[arg(short = 'p', long, value_name = "PrimerFile")]
    pub primers: Option<String>,

    #[arg(long = "p1", value_name = "FwdPrimer")]
    pub fwd_primer: Option<String>,

    #[arg(long = "p2", value_name = "RevPrimer")]
    pub rev_primer: Option<String>,

    #[arg(short = 't', long, value_name = "TagFile", required = true)]
    pub tags: String,

    #[arg(short = 's', long = "sampleInfo", value_name = "SampleFile", required = true)]
    pub sample_info: String,

    #[arg(short = 'l', long = "pool", value_name = "PoolFile", required = true)]
    pub pool: String,

    #[arg(short = 'm', long, default_value_t = false)]
    pub allow_multiple_primers: bool,

    #[arg(long = "pm", value_name = "PrimerMismatches", default_value_t = 0)]
    pub primer_mismatches: usize,

    #[arg(long = "tm", value_name = "TagMismatches", default_value_t = 0)]
    pub tag_mismatches: usize,

    #[arg(short = 'd', long, value_name = "OutDirectory", default_value = ".")]
    pub output_directory: String,

    #[arg(short = 'o', long, value_name = "OutPrefix", default_value = "")]
    pub output_prefix: String,
}

#[derive(Parser, Debug)]
pub struct FilterArgs {
    #[arg(short = 'i', long = "inputPrefix", value_name = "InputPrefix", required = true)]
    pub input_prefix: String,

    #[arg(short = 's', long = "sampleInfo", value_name = "SampleFile", required = true)]
    pub sample_info: String,

    #[arg(short = 'p', long, value_name = "propPCRs", default_value_t = 1.0)]
    pub prop_pcrs: f64,

    #[arg(short = 'm', long, value_name = "minTimes", default_value_t = 1)]
    pub min_occurrence: usize,

    #[arg(short = 'l', long, value_name = "minLength", default_value_t = 1)]
    pub min_length: usize,

    #[arg(short = 'd', long, value_name = "OutDirectory", default_value = ".")]
    pub output_directory: String,

    #[arg(short = 'o', long, value_name = "OutPrefix", default_value = "Filtered")]
    pub output_prefix: String,
}

fn main() -> Result<()> {
    env_logger::init();
    let cli = Cli::parse();
    match cli.command {
        Commands::Sort(args) => sort::run(args),
        Commands::Filter(args) => filter::run(args),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn test_cli_config_valid() {
        Cli::command().debug_assert();
    }

    #[test]
    fn test_sort_subcommand_parses() {
        let args = Cli::try_parse_from([
            "begum", "sort",
            "-t", "tags.txt",
            "-s", "samples.txt",
            "-l", "pools.txt",
            "-p", "primers.txt",
        ]).unwrap();
        match args.command {
            Commands::Sort(s) => {
                assert_eq!(s.tags, "tags.txt");
                assert_eq!(s.primers, Some("primers.txt".to_string()));
                assert_eq!(s.primer_mismatches, 0);
                assert_eq!(s.tag_mismatches, 0);
                assert!(!s.allow_multiple_primers);
                assert_eq!(s.output_directory, ".");
                assert_eq!(s.output_prefix, "");
            }
            _ => panic!("Expected Sort"),
        }
    }

    #[test]
    fn test_filter_subcommand_parses() {
        let args = Cli::try_parse_from([
            "begum", "filter",
            "-i", "out/test",
            "-s", "samples.txt",
            "-p", "0.5",
            "-m", "2",
            "-l", "10",
        ]).unwrap();
        match args.command {
            Commands::Filter(f) => {
                assert_eq!(f.input_prefix, "out/test");
                assert!((f.prop_pcrs - 0.5).abs() < 1e-9);
                assert_eq!(f.min_occurrence, 2);
                assert_eq!(f.min_length, 10);
                assert_eq!(f.output_prefix, "Filtered");
            }
            _ => panic!("Expected Filter"),
        }
    }
}
