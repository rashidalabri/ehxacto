use std::fmt;

use serde::{Serialize, Serializer};

#[derive(Debug, Clone)]
pub struct ReferenceRegion {
    pub contig: String,
    pub start: i64,
    pub end: i64,
}

impl fmt::Display for ReferenceRegion {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}:{}-{}", self.contig, self.start, self.end)
    }
}

impl Serialize for ReferenceRegion {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let formatted_string = format!("{}:{}-{}", self.contig, self.start, self.end);
        serializer.serialize_str(&formatted_string)
    }
}
