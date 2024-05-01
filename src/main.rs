mod catalog;
mod find;
mod reference;
mod region;
mod types;
mod workflow;

use anyhow::{Ok, Result};
use catalog::write_locus_catalog;
use clap::{command, Parser};
use std::{collections::HashMap, path::PathBuf};
use workflow::process_regions;

use crate::{catalog::LocusNode, workflow::load_region_catalog_and_extend};

/// Find exact tandem repeat expansion regions from ExpansionHunter denovo output
#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// Path to reference genome
    reference: PathBuf,

    /// Path to input catalog file (tab-delimited)
    input_catalog: PathBuf,

    /// Path to outout catalog file (json)
    output_catalog: PathBuf,

    /// Number of threads
    #[clap(short, long, default_value_t = 1)]
    threads: usize,

    /// Extend the search region by this many bases
    #[clap(short, long, default_value_t = 1000)]
    extend: u64,

    /// Minimum number of consecutive pure repeat units
    /// required to be in the reference genome
    #[clap(short = 'u', long, default_value_t = 2)]
    min_units: usize,

    /// Maximum distance between two pure repeat tracts
    /// to be merged into one prior to running DP algorithm
    #[clap(short = 'd', long, default_value_t = 20)]
    max_distance: usize,

    /// Minimum repeat purity score for a repeat variant
    /// to be outputted in the catalog
    #[clap(short = 'p', long, default_value = "0.90")]
    min_purity: f64,

    /// Maximum number of variants
    #[clap(short = 'n', long, default_value_t = usize::MAX)]
    max_variants: usize,

    /// Maximum length of locus structure field in the
    /// output catalog
    #[clap(short = 'l', long, default_value_t = usize::MAX)]
    max_locus_length: usize,

    /// Output both the original and reverse complement repeat regions.
    /// Otherwise, only output the direction with the least variants. If
    /// there is a tie, output the original direction.
    #[clap(long)]
    output_rev_comp: bool,
}

fn main() -> Result<()> {
    let args = Args::parse();

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .unwrap();

    let regions = load_region_catalog_and_extend(&args)?;
    // .into_iter()
    // .skip(10)
    // .take(1)
    // .collect();

    let max_nodes = args.max_variants * 2 - 1;
    let loci = process_regions(&regions, &args)?
        .into_iter()
        .filter(|locus| {
            let locus_structure_len = locus.structure().len();
            locus.nodes.len() > 0
                && locus.nodes.len() <= max_nodes
                && locus_structure_len <= args.max_locus_length
        })
        .collect();

    write_locus_catalog(&loci, &args.output_catalog)?;

    eprintln!("Done!");

    // print summary
    let n_regions = regions.len();
    let n_loci = loci.len();
    let pct_loci = (n_loci as f64 / n_regions as f64) * 100.0;

    let mut n_nodes: HashMap<usize, usize> = HashMap::new();
    for locus in loci.iter() {
        let n_variant_nodes = locus
            .nodes
            .iter()
            .filter(|node| {
                if let LocusNode::Variant(_) = node {
                    true
                } else {
                    false
                }
            })
            .count();
        let entry = n_nodes.entry(n_variant_nodes).or_insert(0);
        *entry += 1;
    }

    eprintln!(
        "Mapped {} / {} ({:.0}%) EHdn regions to repeat loci in reference genome",
        n_loci, n_regions, pct_loci
    );

    let mut sorted_entries: Vec<(&usize, &usize)> = n_nodes.iter().collect();
    sorted_entries.sort_by_key(|&(k, _)| *k);

    let mut sorted_entries: Vec<(&usize, &usize)> = n_nodes.iter().collect();
    sorted_entries.sort_by_key(|&(k, _)| *k);

    // calculate average length of nodes
    let mut total_num_variant_nodes = 0;
    let avg_length: f64 = loci
        .iter()
        .map(|locus| {
            let mut num_variant_nodes = 0;
            let total_length: usize = locus
                .nodes
                .iter()
                .map(|node| match node {
                    LocusNode::NonVariant(_) => 0,
                    LocusNode::Variant(var) => {
                        num_variant_nodes += 1;
                        (var.region.end - var.region.start) as usize
                    }
                })
                .sum();
            total_num_variant_nodes += num_variant_nodes;
            total_length as f64 / num_variant_nodes as f64
        })
        .sum::<f64>()
        / total_num_variant_nodes as f64;

    eprintln!(
        "Avg. length of mapped repeat tracts found in reference genome: {:.2}",
        avg_length
    );

    println!("+---------------------+-------------+---------+");
    println!("| No. of repeat nodes | No. of loci |    Pct. |");
    println!("+---------------------+-------------+---------+");

    for (key, value) in sorted_entries {
        let percentage = (*value as f64 / n_regions as f64) * 100.0;
        println!("| {:19} | {:11} | {:6.2}% |", key, value, percentage);
    }

    println!("+---------------------+-------------+---------+");

    Ok(())
}
