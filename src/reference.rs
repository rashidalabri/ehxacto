use std::{fs::File, path::Path};

use bio::io::fasta;

pub type Reference = fasta::IndexedReader<File>;

pub fn load_reference(path: &Path) -> fasta::IndexedReader<File> {
    fasta::IndexedReader::from_file(&path).unwrap()
}
