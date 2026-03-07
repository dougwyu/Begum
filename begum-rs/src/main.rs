mod args;
mod dna;
mod filter;
mod sort;

use anyhow::Result;
use args::{FilterArgs, SortArgs};
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
