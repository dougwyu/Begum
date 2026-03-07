use clap::Parser;

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
