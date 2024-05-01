use serde::de::{self, Visitor};
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::convert::From;
use std::fmt;
use std::ops::Deref;
use std::ops::DerefMut;

#[derive(PartialEq, Eq, Debug, Clone)]
pub struct Sequence(Vec<u8>);

impl Sequence {
    // Helper function to convert to String
    pub fn to_string(&self) -> String {
        String::from_utf8_lossy(&self.0).into_owned()
    }

    pub fn reverse_complement(&self) -> Self {
        let revcomp = self.0.iter().rev().map(|&c| match c {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => c,
        });

        Sequence(revcomp.collect())
    }
}

impl From<&[u8]> for Sequence {
    fn from(item: &[u8]) -> Self {
        Sequence(item.to_vec())
    }
}

impl From<Vec<u8>> for Sequence {
    fn from(item: Vec<u8>) -> Self {
        Sequence(item)
    }
}

impl From<&Vec<u8>> for Sequence {
    fn from(item: &Vec<u8>) -> Self {
        Sequence(item.clone())
    }
}

impl From<&str> for Sequence {
    fn from(item: &str) -> Self {
        Sequence(item.as_bytes().to_vec())
    }
}

impl From<String> for Sequence {
    fn from(item: String) -> Self {
        Sequence(item.into_bytes())
    }
}

impl AsRef<str> for Sequence {
    fn as_ref(&self) -> &str {
        std::str::from_utf8(&self.0).unwrap_or_default()
    }
}

impl AsRef<[u8]> for Sequence {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}

impl Deref for Sequence {
    type Target = Vec<u8>;

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl DerefMut for Sequence {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.0
    }
}

impl fmt::Display for Sequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_string())
    }
}

impl Serialize for Sequence {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        serializer.serialize_str(self.as_ref())
    }
}

struct SequenceVisitor;

impl<'de> Visitor<'de> for SequenceVisitor {
    type Value = Sequence;

    fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
        formatter.write_str("a DNA sequence as a string")
    }

    fn visit_str<E>(self, value: &str) -> Result<Self::Value, E>
    where
        E: de::Error,
    {
        Ok(Sequence::from(value))
    }
}

impl<'de> Deserialize<'de> for Sequence {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        deserializer.deserialize_str(SequenceVisitor)
    }
}

impl IntoIterator for Sequence {
    type Item = u8;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.into_iter()
    }
}

impl<'a> IntoIterator for &'a Sequence {
    type Item = &'a u8;
    type IntoIter = std::slice::Iter<'a, u8>;

    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

// Implement from iterator <u8> for Sequence
impl From<std::vec::IntoIter<u8>> for Sequence {
    fn from(item: std::vec::IntoIter<u8>) -> Self {
        Sequence(item.collect())
    }
}

impl<const N: usize> From<&[u8; N]> for Sequence {
    fn from(item: &[u8; N]) -> Self {
        Sequence(item.to_vec())
    }
}
