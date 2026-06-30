//! TGT reader — text and binary input for TGT records.
//!
//! Provides `TgtReader` for reading `TgtRecord` instances from either
//! human-readable text format or compact binary format.
//!
//! The text format parser can parse the output of `TgtWriter::write_record`,
//! enabling round-trip conversion between text and binary.

use crate::enzyme::enzyme::EnzymeType;
use crate::tgt::gap::Gap;
use crate::tgt::record::TgtRecord;
use crate::tgt::tag::{Strand, Tag};
use anyhow::{bail, Context, Result};
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

/// Reader for TGT records (text or binary format).
pub struct TgtReader {
    reader: BufReader<File>,
    /// Buffer holding the next header line, if peeked.
    peeked_line: Option<String>,
}

impl TgtReader {
    /// Create a new TGT reader that reads from the given file path.
    pub fn new(path: &Path) -> Result<Self> {
        let file = File::open(path)
            .with_context(|| format!("Failed to open TGT input file: {}", path.display()))?;
        Ok(Self {
            reader: BufReader::new(file),
            peeked_line: None,
        })
    }

    /// Read the next TGT record from the file in text format.
    ///
    /// Returns `Ok(None)` when EOF is reached.
    ///
    /// # Format
    /// ```text
    /// >genome_id|length=4641652
    /// BcgI:ATCG... -1313- GCTA... -1298- ...
    /// AlfI:CGAT... -892- ...
    /// ```
    pub fn read_record(&mut self) -> Result<Option<TgtRecord>> {
        // Read the header line
        let header_line = match self.read_next_non_empty_line()? {
            Some(line) => line,
            None => return Ok(None),
        };

        // Parse header: >genome_id|length=NNN
        let (genome_id, total_length) = parse_header(&header_line)?;
        let mut record = TgtRecord::new(&genome_id, total_length);

        // Read tag lines until we hit the next header or EOF
        loop {
            match self.peek_line()? {
                Some(line) => {
                    if line.starts_with('>') {
                        // Next record starts
                        break;
                    }
                    // Consume the line
                    let line = self.read_next_line()?.unwrap();
                    let trimmed = line.trim();
                    if trimmed.is_empty() {
                        continue;
                    }
                    parse_tag_line(trimmed, &mut record)?;
                }
                None => break, // EOF
            }
        }

        Ok(Some(record))
    }

    /// Read a single TGT record in binary format.
    ///
    /// Binary format layout (per the spec):
    /// - Header (32 bytes): magic, version, genome length, tag count, enzyme count
    /// - Tag table (N x 48 bytes each)
    /// - Gap table ((N-1) x 4 bytes each)
    ///
    /// Returns `Ok(None)` when EOF is reached (no more records to read).
    pub fn read_binary(&mut self) -> Result<Option<TgtRecord>> {
        let mut header_buf = [0u8; 32];
        match self.reader.read_exact(&mut header_buf) {
            Ok(()) => {}
            Err(e) => {
                // Check if we just hit EOF (no more records)
                if e.kind() == std::io::ErrorKind::UnexpectedEof {
                    return Ok(None);
                }
                return Err(e).context("Failed to read binary TGT header");
            }
        }

        // Verify magic bytes
        if &header_buf[0..4] != b"TGT\x01" {
            bail!(
                "Invalid binary TGT magic bytes: expected b'TGT\\x01', got {:?}",
                &header_buf[0..4]
            );
        }

        // Parse version
        let version = u32::from_le_bytes([
            header_buf[4], header_buf[5], header_buf[6], header_buf[7],
        ]);
        if version != 1 {
            bail!("Unsupported binary TGT version: {} (expected 1)", version);
        }

        // Parse genome length
        let total_length = u64::from_le_bytes([
            header_buf[8], header_buf[9], header_buf[10], header_buf[11],
            header_buf[12], header_buf[13], header_buf[14], header_buf[15],
        ]);

        // Parse tag count
        let tag_count = u32::from_le_bytes([
            header_buf[16], header_buf[17], header_buf[18], header_buf[19],
        ]);

        // Enzyme count is present but not needed for reading (bytes 20..22)
        // Reserved bytes 22..32 are ignored

        // Create record — we don't know the genome_id from binary format,
        // so we use a default empty string. The caller can set it if needed.
        let mut record = TgtRecord::new("", total_length);

        // Read tag table
        for _ in 0..tag_count {
            let mut tag_buf = [0u8; 48];
            self.reader
                .read_exact(&mut tag_buf)
                .context("Failed to read binary tag entry")?;

            // Parse sequence (bytes 0..32)
            let mut sequence = [0u8; 32];
            sequence.copy_from_slice(&tag_buf[0..32]);

            // Parse position (bytes 32..40)
            let position = u64::from_le_bytes([
                tag_buf[32], tag_buf[33], tag_buf[34], tag_buf[35],
                tag_buf[36], tag_buf[37], tag_buf[38], tag_buf[39],
            ]);

            // Parse enzyme index (byte 40)
            let enzyme = EnzymeType::from_index(tag_buf[40])
                .with_context(|| format!("Invalid enzyme index: {}", tag_buf[40]))?;

            // Parse strand (byte 41)
            let strand = Strand::from_u8(tag_buf[41])
                .with_context(|| format!("Invalid strand value: {}", tag_buf[41]))?;

            // Reserved bytes 42..48 are ignored

            let tag = Tag::new(sequence, position, enzyme, strand);
            record.add_tag(tag);
        }

        // Read gap table (tag_count - 1 gaps)
        let gap_count = if tag_count > 0 {
            tag_count - 1
        } else {
            0
        };
        // Collect stored gaps for verification against auto-computed gaps
        let mut stored_gaps = Vec::with_capacity(gap_count as usize);
        for _ in 0..gap_count {
            let mut gap_buf = [0u8; 4];
            self.reader
                .read_exact(&mut gap_buf)
                .context("Failed to read binary gap entry")?;
            let size = u32::from_le_bytes([gap_buf[0], gap_buf[1], gap_buf[2], gap_buf[3]]);
            stored_gaps.push(size);
        }

        // Verify stored gaps match auto-computed gaps from add_tag()
        if stored_gaps.len() != record.gaps.len() {
            bail!(
                "Gap count mismatch: stored={}, computed={}",
                stored_gaps.len(),
                record.gaps.len()
            );
        }
        for (i, (stored, computed)) in stored_gaps.iter().zip(record.gaps.iter()).enumerate() {
            if *stored != computed.size {
                bail!(
                    "Gap mismatch at index {}: stored={}, computed={}",
                    i, stored, computed.size
                );
            }
        }

        Ok(Some(record))
    }

    // --- Helper methods ---

    /// Read the next non-empty line from the file.
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

    /// Read the next line from the file (or peeked buffer).
    fn read_next_line(&mut self) -> Result<Option<String>> {
        if let Some(line) = self.peeked_line.take() {
            return Ok(Some(line));
        }
        let mut buf = String::new();
        let n = self
            .reader
            .read_line(&mut buf)
            .context("Failed to read line from TGT file")?;
        if n == 0 {
            Ok(None)
        } else {
            // Trim trailing newline but preserve other content
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

/// Parse a header line of the form `>genome_id|length=NNN`
fn parse_header(line: &str) -> Result<(String, u64)> {
    if !line.starts_with('>') {
        bail!("Expected header line starting with '>', got: {}", line);
    }
    let content = &line[1..]; // strip leading '>'

    // Split on '|'
    let parts: Vec<&str> = content.splitn(2, '|').collect();
    let genome_id = parts[0].trim().to_string();

    let mut total_length = 0u64;
    if parts.len() > 1 {
        let attr_part = parts[1];
        // Parse length=NNN
        for attr in attr_part.split(',') {
            let attr = attr.trim();
            if let Some(val) = attr.strip_prefix("length=") {
                total_length = val
                    .parse::<u64>()
                    .with_context(|| format!("Invalid length value: {}", val))?;
            }
        }
    }

    Ok((genome_id, total_length))
}

/// Parse a tag line of the form `Enzyme:SEQ [-gap- Enzyme:SEQ]*`
///
/// Each tag is prefixed with its enzyme type. Gaps between consecutive tags
/// are verified against the parsed gap tokens.
fn parse_tag_line(line: &str, record: &mut TgtRecord) -> Result<()> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    // Expected: [Enzyme:SEQ, -NNN-, Enzyme:SEQ, -NNN-, ...]

    let mut i = 0;
    while i < tokens.len() {
        // Each tag token must be of the form "Enzyme:SEQ"
        let tag_token = tokens[i];
        let colon_idx = tag_token
            .find(':')
            .with_context(|| format!("Expected enzyme prefix (Enzyme:SEQ), got: {}", tag_token))?;
        let enzyme_name = &tag_token[..colon_idx];
        let enzyme = parse_enzyme_type(enzyme_name)?;
        let seq_str = &tag_token[colon_idx + 1..];
        let sequence = seq_to_bytes(seq_str)?;

        // Compute position: for parsed records, we reconstruct positions from gaps.
        let position = if record.tags.is_empty() {
            0
        } else {
            let last_pos = record.tags.last().unwrap().position;
            let last_gap = record.gaps.last().map(|g| g.size as u64).unwrap_or(0);
            last_pos + last_gap
        };

        let tag = Tag::new(sequence, position, enzyme, Strand::Forward);
        record.add_tag(tag);

        i += 1;

        // Check for gap token after this tag
        if i < tokens.len() {
            let gap_token = tokens[i];
            if gap_token.starts_with('-') && gap_token.ends_with('-') {
                // Parse and validate the gap value
                let gap_inner = &gap_token[1..gap_token.len() - 1];
                let parsed_gap: u32 = gap_inner
                    .parse()
                    .with_context(|| format!("Invalid gap value: {}", gap_token))?;

                // The gap was already computed by add_tag(). Verify it matches.
                let gap_idx = record.gaps.len().saturating_sub(1);
                if gap_idx < record.gaps.len() {
                    let computed_gap = record.gaps[gap_idx].size;
                    if computed_gap != parsed_gap {
                        bail!(
                            "Gap mismatch at index {}: computed={}, parsed={}",
                            gap_idx, computed_gap, parsed_gap
                        );
                    }
                }
                i += 1;
            }
        }
    }

    Ok(())
}

/// Convert an enzyme name string to EnzymeType.
fn parse_enzyme_type(name: &str) -> Result<EnzymeType> {
    match name.trim() {
        "BcgI" => Ok(EnzymeType::BcgI),
        "AlfI" => Ok(EnzymeType::AlfI),
        "AloI" => Ok(EnzymeType::AloI),
        "BaeI" => Ok(EnzymeType::BaeI),
        "BplI" => Ok(EnzymeType::BplI),
        "BsaXI" => Ok(EnzymeType::BsaXI),
        "BslFI" => Ok(EnzymeType::BslFI),
        "Bsp24I" => Ok(EnzymeType::Bsp24I),
        "CjeI" => Ok(EnzymeType::CjeI),
        "CjePI" => Ok(EnzymeType::CjePI),
        "CspCI" => Ok(EnzymeType::CspCI),
        "FalI" => Ok(EnzymeType::FalI),
        "HaeIV" => Ok(EnzymeType::HaeIV),
        "Hin4I" => Ok(EnzymeType::Hin4I),
        "PpiI" => Ok(EnzymeType::PpiI),
        "PsrI" => Ok(EnzymeType::PsrI),
        other => bail!("Unknown enzyme type: {}", other),
    }
}

/// Convert a DNA sequence string to a 32-byte padded array.
fn seq_to_bytes(seq: &str) -> Result<[u8; 32]> {
    if seq.len() > 32 {
        bail!("Sequence too long: {} (max 32)", seq.len());
    }
    let mut arr = [0u8; 32];
    let bytes = seq.as_bytes();
    // Validate: only A, C, G, T allowed
    for &b in bytes {
        if !matches!(b, b'A' | b'C' | b'G' | b'T' | b'a' | b'c' | b'g' | b't') {
            bail!("Invalid DNA base in sequence: {}", b as char);
        }
    }
    arr[..bytes.len()].copy_from_slice(bytes);
    Ok(arr)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    fn make_seq(s: &str) -> [u8; 32] {
        let mut arr = [0u8; 32];
        let bytes = s.as_bytes();
        arr[..bytes.len()].copy_from_slice(bytes);
        arr
    }

    fn make_tag(seq: &str, pos: u64, enzyme: EnzymeType) -> Tag {
        Tag::new(make_seq(seq), pos, enzyme, Strand::Forward)
    }

    #[test]
    fn test_parse_header() {
        let (id, len) = parse_header(">NC_000913|length=4641652").unwrap();
        assert_eq!(id, "NC_000913");
        assert_eq!(len, 4_641_652);
    }

    #[test]
    fn test_parse_header_no_length() {
        let (id, len) = parse_header(">genome001").unwrap();
        assert_eq!(id, "genome001");
        assert_eq!(len, 0);
    }

    #[test]
    fn test_parse_header_invalid() {
        assert!(parse_header("NC_000913|length=4641652").is_err());
    }

    #[test]
    fn test_parse_enzyme_type() {
        assert_eq!(parse_enzyme_type("BcgI").unwrap(), EnzymeType::BcgI);
        assert_eq!(parse_enzyme_type("AlfI").unwrap(), EnzymeType::AlfI);
        assert_eq!(parse_enzyme_type("PsrI").unwrap(), EnzymeType::PsrI);
        assert!(parse_enzyme_type("Unknown").is_err());
    }

    #[test]
    fn test_seq_to_bytes() {
        let arr = seq_to_bytes("ATCG").unwrap();
        assert_eq!(&arr[..4], b"ATCG");
        assert_eq!(arr[4..].iter().all(|&b| b == 0), true);
    }

    #[test]
    fn test_seq_to_bytes_too_long() {
        assert!(seq_to_bytes("A".repeat(33).as_str()).is_err());
    }

    #[test]
    fn test_read_record_single() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.tgt");
        {
            let mut f = File::create(&path).unwrap();
            writeln!(f, ">G001|length=10000").unwrap();
            writeln!(f, "BcgI:ATCG -500- BcgI:GCTA").unwrap();
        }

        let mut reader = TgtReader::new(&path).unwrap();
        let record = reader.read_record().unwrap().unwrap();
        assert_eq!(record.genome_id, "G001");
        assert_eq!(record.total_length, 10000);
        assert_eq!(record.tag_count(), 2);
        assert_eq!(record.tags[0].sequence_str(), "ATCG");
        assert_eq!(record.tags[0].enzyme, EnzymeType::BcgI);
        assert_eq!(record.tags[1].sequence_str(), "GCTA");
        assert_eq!(record.tags[1].enzyme, EnzymeType::BcgI);
        assert_eq!(record.gaps.len(), 1);
        assert_eq!(record.gaps[0].size, 500);
    }

    #[test]
    fn test_read_record_multiple() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.tgt");
        {
            let mut f = File::create(&path).unwrap();
            writeln!(f, ">G001|length=5000").unwrap();
            writeln!(f, "BcgI:ATCG -500- BcgI:GCTA").unwrap();
            writeln!(f, ">G002|length=6000").unwrap();
            writeln!(f, "AlfI:TTAA -300- AlfI:CCGG").unwrap();
        }

        let mut reader = TgtReader::new(&path).unwrap();

        let record1 = reader.read_record().unwrap().unwrap();
        assert_eq!(record1.genome_id, "G001");
        assert_eq!(record1.tag_count(), 2);

        let record2 = reader.read_record().unwrap().unwrap();
        assert_eq!(record2.genome_id, "G002");
        assert_eq!(record2.tag_count(), 2);

        assert!(reader.read_record().unwrap().is_none());
    }

    #[test]
    fn test_read_record_eof() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("empty.tgt");
        File::create(&path).unwrap();

        let mut reader = TgtReader::new(&path).unwrap();
        assert!(reader.read_record().unwrap().is_none());
    }
}
