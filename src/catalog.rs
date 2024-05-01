use crate::region::ReferenceRegion;
use crate::types::Sequence;
use anyhow::Result;
use serde::ser;
use serde::ser::SerializeStruct;
use serde::{Deserialize, Serialize, Serializer};
use serde_json::to_string_pretty;
use std::fmt;
use std::fs::File;
use std::io::Write;
use std::path::Path;

pub type VariantId = String;

#[derive(Debug, Clone)]
pub enum VariantNodeDefinition {
    ShortTandemRepeatZeroOrMore(Sequence),
}

impl fmt::Display for VariantNodeDefinition {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            VariantNodeDefinition::ShortTandemRepeatZeroOrMore(seq) => write!(f, "({})*", seq),
        }
    }
}

impl serde::Serialize for VariantNodeDefinition {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let variant_str = format!("{}", self);
        serializer.serialize_str(&variant_str)
    }
}

#[derive(Debug, Clone, Serialize)]
pub enum VariantType {
    Repeat,
}

#[derive(Debug, Clone, Serialize)]
pub struct Variant {
    pub id: VariantId,
    pub region: ReferenceRegion,
    pub definition: VariantNodeDefinition,
    pub variant_type: VariantType,
}

pub type LocusId = String;

#[derive(Debug, Clone)]
pub enum LocusNode {
    NonVariant(Sequence),
    Variant(Variant),
}

#[derive(Debug, Clone)]
pub struct RepeatLocus {
    pub id: LocusId,
    pub nodes: Vec<LocusNode>,
}

impl RepeatLocus {
    pub fn structure(&self) -> String {
        self.nodes
            .iter()
            .map(|node| match node {
                LocusNode::NonVariant(seq) => format!("{}", seq),
                LocusNode::Variant(var) => format!("{}", var.definition),
            })
            .collect::<Vec<_>>()
            .join("")
    }
}

impl ser::Serialize for RepeatLocus {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("Locus", 4)?;

        // Serialize the Locus ID
        state.serialize_field("LocusId", &self.id)?;

        // Concatenate all node representations for LocusStructure
        let locus_structure = self
            .nodes
            .iter()
            .map(|node| match node {
                LocusNode::NonVariant(seq) => format!("{}", seq),
                LocusNode::Variant(var) => format!("{}", var.definition),
            })
            .collect::<Vec<_>>()
            .join("");
        state.serialize_field("LocusStructure", &locus_structure)?;

        // Collect and serialize ReferenceRegions and VariantTypes from Variants only
        let (ids, (regions, types)): (Vec<_>, (Vec<_>, Vec<_>)) = self
            .nodes
            .iter()
            .filter_map(|node| {
                if let LocusNode::Variant(var) = node {
                    Some((
                        var.id.to_string(),
                        (format!("{}", var.region), format!("{:?}", var.variant_type)),
                    ))
                } else {
                    None
                }
            })
            .unzip();

        if self.nodes.len() == 1 {
            state.serialize_field("ReferenceRegion", &regions[0])?;
            state.serialize_field("VariantId", &ids[0])?;
            state.serialize_field("VariantType", &types[0])?;
        } else {
            state.serialize_field("ReferenceRegion", &regions)?;
            state.serialize_field("VariantId", &ids)?;
            state.serialize_field("VariantType", &types)?;
        }

        state.end()
    }
}

#[derive(Debug, Deserialize, Clone)]
pub struct RepeatRegion {
    pub id: String,
    pub contig: String,
    pub start: i64,
    pub end: i64,
    pub motif: Sequence,
}

pub fn load_region_catalog(path: &Path) -> Result<Vec<RepeatRegion>> {
    let mut rdr = csv::ReaderBuilder::new().delimiter(b'\t').from_path(path)?;
    rdr.deserialize()
        .map(|result| result.map_err(Into::into))
        .collect()
}

pub fn write_locus_catalog(loci: &Vec<RepeatLocus>, file_path: &Path) -> Result<()> {
    let json = to_string_pretty(&loci)?;
    let mut file = File::create(file_path)?;
    file.write_all(json.as_bytes())?;
    Ok(())
}
