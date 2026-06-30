//! FASTA parser — streaming parser for large genome sequences.

use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// A single FASTA record.
#[derive(Clone, Debug, PartialEq)]
pub struct FastaRecord {
    /// Sequence identifier (without the leading '>').
    pub id: String,
    /// DNA sequence (uppercase bytes).
    pub sequence: Vec<u8>,
    /// Optional description line content.
    pub description: Option<String>,
}

/// Streaming FASTA reader for large genomes.
///
/// Reads one record at a time without loading the entire file into memory.
pub struct FastaReader {
    reader: BufReader<File>,
    peeked_line: Option<String>,
}

impl FastaReader {
    /// Create a new FASTA reader for the given file path.
    pub fn new(path: &Path) -> Result<Self> {
        let file = File::open(path)
            .with_context(|| format!("Failed to open FASTA file: {}", path.display()))?;
        Ok(Self {
            reader: BufReader::new(file),
            peeked_line: None,
        })
    }

    /// Read the next FASTA record.
    ///
    /// Returns `Ok(None)` when EOF is reached.
    pub fn next_record(&mut self) -> Result<Option<FastaRecord>> {
        // Read header line (starts with '>')
        let header = match self.read_next_non_empty_line()? {
            Some(line) => line,
            None => return Ok(None),
        };

        if !header.starts_with('>') {
            anyhow::bail!("Expected FASTA header line starting with '>', got: {}", header);
        }

        // Parse ID and description from header
        let header_content = &header[1..];
        let (id, description) = match header_content.find(' ') {
            Some(idx) => {
                let id = header_content[..idx].trim().to_string();
                let desc = header_content[idx + 1..].trim().to_string();
                (id, Some(desc))
            }
            None => (header_content.trim().to_string(), None),
        };

        // Read sequence lines until next header or EOF
        let mut sequence = Vec::new();
        loop {
            match self.peek_line()? {
                Some(line) => {
                    if line.starts_with('>') {
                        break;
                    }
                    let line = self.read_next_line()?.unwrap();
                    // Convert to uppercase and filter to ACGT only
                    sequence.extend(
                        line.bytes()
                            .filter(|b| matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't' | b'N' | b'n'))
                            .map(|b| b.to_ascii_uppercase()),
                    );
                }
                None => break,
            }
        }

        Ok(Some(FastaRecord {
            id,
            sequence,
            description,
        }))
    }

    /// Read the next non-empty line.
    fn read_next_non_empty_line(&mut self) -> Result<Option<String>> {
        loop {
            match self.read_next_line()? {
                Some(line) => {
                    if !line.trim().is_empty() {
                        return Ok(Some(line));
                    }
                }
                None => return Ok(None),
            }
        }
    }

    /// Read the next line (or return peeked line).
    fn read_next_line(&mut self) -> Result<Option<String>> {
        if let Some(line) = self.peeked_line.take() {
            return Ok(Some(line));
        }
        let mut buf = String::new();
        let n = self
            .reader
            .read_line(&mut buf)
            .context("Failed to read line from FASTA file")?;
        if n == 0 {
            Ok(None)
        } else {
            if buf.ends_with('\n') {
                buf.pop();
                if buf.ends_with('\r') {
                    buf.pop();
                }
            }
            Ok(Some(buf))
        }
    }

    /// Peek at the next line without consuming it.
    fn peek_line(&mut self) -> Result<Option<&str>> {
        if self.peeked_line.is_none() {
            self.peeked_line = self.read_next_line()?;
        }
        Ok(self.peeked_line.as_deref())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    #[test]
    fn test_fasta_reader_single_record() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.fasta");
        {
            let mut f = File::create(&path).unwrap();
            writeln!(f, ">seq1 description here").unwrap();
            writeln!(f, "ATCGATCG").unwrap();
            writeln!(f, "ATCGATCG").unwrap();
        }

        let mut reader = FastaReader::new(&path).unwrap();
        let record = reader.next_record().unwrap().unwrap();
        assert_eq!(record.id, "seq1");
        assert_eq!(record.description, Some("description here".to_string()));
        assert_eq!(record.sequence, b"ATCGATCGATCGATCG");
    }

    #[test]
    fn test_fasta_reader_multiple_records() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.fasta");
        {
            let mut f = File::create(&path).unwrap();
            writeln!(f, ">seq1").unwrap();
            writeln!(f, "AAA").unwrap();
            writeln!(f, ">seq2").unwrap();
            writeln!(f, "TTT").unwrap();
        }

        let mut reader = FastaReader::new(&path).unwrap();
        let r1 = reader.next_record().unwrap().unwrap();
        assert_eq!(r1.id, "seq1");
        assert_eq!(r1.sequence, b"AAA");

        let r2 = reader.next_record().unwrap().unwrap();
        assert_eq!(r2.id, "seq2");
        assert_eq!(r2.sequence, b"TTT");

        assert!(reader.next_record().unwrap().is_none());
    }

    #[test]
    fn test_fasta_reader_empty() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("empty.fasta");
        File::create(&path).unwrap();

        let mut reader = FastaReader::new(&path).unwrap();
        assert!(reader.next_record().unwrap().is_none());
    }
}
