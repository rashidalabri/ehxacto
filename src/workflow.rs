use crate::catalog::load_region_catalog;
use crate::find::find_repeat_variants;
use crate::reference::load_reference;
use crate::Args;
use crate::{
    catalog::{LocusNode, RepeatLocus, RepeatRegion, Variant, VariantNodeDefinition, VariantType},
    reference::Reference,
    region::ReferenceRegion,
    types::Sequence,
};
use anyhow::{Ok, Result};
use indicatif::{ParallelProgressIterator, ProgressStyle};
use rayon::prelude::*;

const EH_EXTENSION_LENGTH: i64 = 1000;
const MAX_NUM_FLANK_NS: usize = 10;

pub fn load_region_catalog_and_extend(args: &Args) -> Result<Vec<RepeatRegion>> {
    let regions = load_region_catalog(&args.input_catalog)?;
    let regions_ext = regions
        .into_iter()
        .map(|mut region| {
            region.start -= args.extend as i64;
            region.end += args.extend as i64;
            region
        })
        .collect();
    Ok(regions_ext)
}

fn passes_n_around_variant_check(reference: &mut Reference, variant: &Variant) -> Result<bool> {
    let contig = &variant.region.contig;
    let (left_start, left_end) = (
        variant.region.start - EH_EXTENSION_LENGTH,
        variant.region.start,
    );
    let (right_start, right_end) = (variant.region.end, variant.region.end + EH_EXTENSION_LENGTH);

    let mut left_flank = Sequence::from("");
    reference.fetch(contig, left_start as u64, left_end as u64)?;
    reference.read(&mut left_flank)?;

    if left_flank.iter().filter(|&&c| c == b'N').count() >= MAX_NUM_FLANK_NS {
        return Ok(false);
    }

    let mut right_flank = Sequence::from("");
    reference.fetch(contig, right_start as u64, right_end as u64)?;
    reference.read(&mut right_flank)?;

    if right_flank.iter().filter(|&&c| c == b'N').count() >= MAX_NUM_FLANK_NS {
        return Ok(false);
    }

    Ok(true)
}

fn fetch_reference_sequence(reference: &mut Reference, region: &RepeatRegion) -> Result<Sequence> {
    let seq_start = region.start as u64;
    let seq_end = region.end as u64;

    let mut ref_seq = Sequence::from("");
    reference.fetch(&region.contig, seq_start, seq_end)?;
    reference.read(&mut ref_seq)?;

    Ok(ref_seq)
}

fn create_variant_node(
    region: &RepeatRegion,
    motif: &Sequence,
    tract: &(usize, usize),
    variant_id: &str,
) -> LocusNode {
    let id = format!("{}/{}", region.id, variant_id);
    let region = ReferenceRegion {
        contig: region.contig.clone(),
        start: region.start as i64 + tract.0 as i64,
        end: region.start as i64 + tract.1 as i64,
    };
    let definition = VariantNodeDefinition::ShortTandemRepeatZeroOrMore(motif.clone());
    let variant_type = VariantType::Repeat;

    let variant = Variant {
        id,
        region,
        definition,
        variant_type,
    };

    LocusNode::Variant(variant)
}

fn create_nonvariant_node(ref_seq: &Sequence, start: usize, end: usize) -> LocusNode {
    let nonvariant_seq: Sequence = ref_seq[start..end].into();
    LocusNode::NonVariant(nonvariant_seq)
}

fn create_repeat_locus(
    region: &RepeatRegion,
    nodes: Vec<LocusNode>,
    add_rev_comp: bool,
) -> RepeatLocus {
    let id = if add_rev_comp {
        format!("{}/rc", region.id)
    } else {
        region.id.clone()
    };
    RepeatLocus { id, nodes }
}

fn process_region_helper(
    reference: &mut Reference,
    region: &RepeatRegion,
    reverse_complement: bool,
    args: &Args,
) -> Result<RepeatLocus> {
    let motif = if reverse_complement {
        region.motif.reverse_complement()
    } else {
        region.motif.clone()
    };

    let ref_seq = fetch_reference_sequence(reference, region)?;

    let tracts = find_repeat_variants(
        &motif,
        &ref_seq,
        args.min_units,
        args.min_purity,
        args.max_distance,
    );

    let mut nodes: Vec<LocusNode> = Vec::new();

    for (i, tract) in tracts.iter().enumerate() {
        let variant_id = i.to_string();
        let variant_node = create_variant_node(region, &motif, tract, &variant_id);
        nodes.push(variant_node);

        if i != tracts.len() - 1 {
            let start = tract.1;
            let end = tracts[i + 1].0;
            let nonvariant_node = create_nonvariant_node(&ref_seq, start, end);
            nodes.push(nonvariant_node);
        }
    }

    for node in nodes.iter() {
        if let LocusNode::Variant(variant) = node {
            if !passes_n_around_variant_check(reference, variant)? {
                return Ok(create_repeat_locus(
                    region,
                    vec![],
                    args.output_rev_comp && reverse_complement,
                ));
            }
        }
    }

    let locus = create_repeat_locus(region, nodes, args.output_rev_comp && reverse_complement);

    Ok(locus)
}

fn process_region(
    reference: &mut Reference,
    region: &RepeatRegion,
    args: &Args,
) -> Result<(Option<RepeatLocus>, Option<RepeatLocus>)> {
    // {
    //     let mut seq = Sequence::from("");
    //     reference.fetch(&"chr2".to_string(), 6253300, 6253395)?;
    //     reference.read(&mut seq)?;

    //     let score = purity_score_hamming(&seq, b"ATATATATACGTGT");
    //     println!("SCORE={}", score);
    //     panic!();
    // }

    let original = process_region_helper(reference, region, false, &args)?;
    let reverse_complement = process_region_helper(reference, region, true, &args)?;

    // Return both directions
    if args.output_rev_comp {
        return Ok((Some(original), Some(reverse_complement)));
    }

    if original.nodes.is_empty() && reverse_complement.nodes.is_empty() {
        return Ok((None, None));
    } else if original.nodes.is_empty() {
        return Ok((None, Some(reverse_complement)));
    } else if reverse_complement.nodes.is_empty() {
        return Ok((Some(original), None));
    } else if original.nodes.len() <= reverse_complement.nodes.len() {
        return Ok((Some(original), None));
    } else {
        return Ok((None, Some(reverse_complement)));
    }
}

pub fn process_regions(regions: &Vec<RepeatRegion>, args: &Args) -> Result<Vec<RepeatLocus>> {
    let style = ProgressStyle::default_bar().template(
        "{spinner:.green} [{elapsed_precise}] [{bar:60.cyan/blue}] {pos}/{len} ({eta})",
    )?;

    let loci: Vec<(Option<RepeatLocus>, Option<RepeatLocus>)> = regions
        .par_iter()
        .progress_with_style(style)
        .map_init(
            || load_reference(&args.reference),
            |reference, region| process_region(reference, region, args),
        )
        .collect::<Result<Vec<(Option<RepeatLocus>, Option<RepeatLocus>)>>>()?;

    let loci: Vec<RepeatLocus> = loci
        .into_iter()
        .flat_map(|(a, b)| vec![a, b])
        .filter_map(|locus| locus)
        .collect();

    Ok(loci)
}
