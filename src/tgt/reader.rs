//! TGT reader — text and binary input for TGT records.
//!
//! Provides `TgtReader` for reading `TgtRecord` instances from either
//! human-readable text format or compact binary format.
//!
//! The text format parser can parse the output of `TgtWriter::write_record`,
//! enabling round-trip conversion between text and binary.

use crate::enzyme::enzyme::EnzymeType;
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

        // Check for optional contig metadata comment line
        if let Some(line) = self.peek_line()? {
            if line.starts_with("#contigs=") {
                let line = self.read_next_line()?.unwrap();
                parse_contig_comment(&line, &mut record)?;
            }
        }

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

    /// Read a single TGT record in binary format (v2).
    ///
    /// Binary format v2 layout:
    /// - Header (48 bytes): magic "TGT\x02", version, genome length, tag count, enzyme count, contig count
    /// - Genome ID (variable): u16 length + bytes
    /// - Tag table (N x 48 bytes each)
    /// - Gap table ((N-1) x 4 bytes each)
    /// - Contig name table (variable): for each contig, u16 name_len + name bytes
    ///
    /// Returns `Ok(None)` when EOF is reached (no more records to read).
    pub fn read_binary(&mut self) -> Result<Option<TgtRecord>> {
        let mut header_buf = [0u8; 48];
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
        if &header_buf[0..4] == b"TGT\x01" {
            bail!(
                "Detected obsolete binary TGT v1 format (magic b'TGT\\x01'). \
                 Please convert to v2 or regenerate with the current toolchain."
            );
        }
        if &header_buf[0..4] != b"TGT\x02" {
            bail!(
                "Invalid binary TGT magic bytes: expected b'TGT\\x02', got {:?}",
                &header_buf[0..4]
            );
        }

        // Parse version
        let version = u32::from_le_bytes([
            header_buf[4], header_buf[5], header_buf[6], header_buf[7],
        ]);
        if version != 2 {
            bail!("Unsupported binary TGT version: {} (expected 2)", version);
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
        let contig_count = u16::from_le_bytes([header_buf[22], header_buf[23]]);

        // Reserved bytes 24..48 are ignored

        // Read genome_id
        let mut id_len_buf = [0u8; 2];
        self.reader.read_exact(&mut id_len_buf)
            .context("Failed to read genome_id length")?;
        let id_len = u16::from_le_bytes(id_len_buf) as usize;
        
        let mut genome_id_bytes = vec![0u8; id_len];
        self.reader.read_exact(&mut genome_id_bytes)
            .context("Failed to read genome_id bytes")?;
        let genome_id = String::from_utf8(genome_id_bytes)
            .unwrap_or_default();
        
        let mut record = TgtRecord::new(&genome_id, total_length);

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
            let strand = Strand::from_u8(tag_buf[41]);

            // Parse contig_id (bytes 42..44)
            let contig_id = u16::from_le_bytes([tag_buf[42], tag_buf[43]]);

            // Reserved bytes 44..48 are ignored

            let tag = Tag::new(sequence, position, enzyme, strand, contig_id);
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

        // Read contig name table
        for _ in 0..contig_count {
            let mut name_len_buf = [0u8; 2];
            self.reader.read_exact(&mut name_len_buf)
                .context("Failed to read contig name length")?;
            let name_len = u16::from_le_bytes(name_len_buf) as usize;
            let mut name_bytes = vec![0u8; name_len];
            self.reader.read_exact(&mut name_bytes)
                .context("Failed to read contig name bytes")?;
            let name = String::from_utf8(name_bytes)
                .unwrap_or_default();
            record.contig_names.push(name);
        }
        // contig_offsets are not stored in the binary format; leave empty

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

/// Parse a contig metadata comment line of the form `#contigs=name1:len1;name2:len2;...`
fn parse_contig_comment(line: &str, record: &mut TgtRecord) -> Result<()> {
    let content = line.trim();
    if let Some(val) = content.strip_prefix("#contigs=") {
        let mut offset = 0u64;
        for part in val.split(';') {
            let part = part.trim();
            if part.is_empty() { continue; }
            let mut split = part.splitn(2, ':');
            let name = split.next().unwrap_or("").trim().to_string();
            let len = split.next().unwrap_or("0").trim().parse::<u64>().unwrap_or(0);
            record.contig_names.push(name);
            record.contig_offsets.push(offset);
            offset += len;
        }
    }
    Ok(())
}

/// Parse a tag line of the form `Enzyme:SEQ@POS[:contig][/strand] [-gap- ...]*`
///
/// Each tag is prefixed with its enzyme type and position. Gaps between
/// consecutive tags are verified against the parsed gap tokens.
fn parse_tag_line(line: &str, record: &mut TgtRecord) -> Result<()> {
    let tokens: Vec<&str> = line.split_whitespace().collect();
    // Expected: [Enzyme:SEQ@POS, -NNN-, Enzyme:SEQ@POS, -NNN-, ...]

    let mut i = 0;
    let mut pending_gap: Option<u32> = None;

    while i < tokens.len() {
        // If current token is a gap, validate the previous gap (if any)
        if tokens[i].starts_with('-') && tokens[i].ends_with('-') {
            let gap_token = tokens[i];
            let gap_inner = &gap_token[1..gap_token.len() - 1];
            let parsed_gap: u32 = gap_inner
                .parse()
                .with_context(|| format!("Invalid gap value: {}", gap_token))?;

            if let Some(expected) = pending_gap.take() {
                if let Some(last_gap) = record.gaps.last() {
                    if last_gap.size != expected {
                        bail!(
                            "Gap mismatch at index {}: computed={}, expected={}",
                            record.gaps.len() - 1, last_gap.size, expected
                        );
                    }
                }
            }
            pending_gap = Some(parsed_gap);
            i += 1;
            continue;
        }

        // Strip the optional trailing strand suffix ("/+" or "/-"). Files written
        // before the field existed simply lack it and default to Forward, so old
        // .tgt files still parse. Only a one-character "+"/"-" suffix is accepted,
        // which keeps contig names containing '/' unambiguous.
        let (tag_token, strand) = match tokens[i].rsplit_once('/') {
            Some((head, "+")) => (head, Strand::Forward),
            Some((head, "-")) => (head, Strand::Reverse),
            _ => (tokens[i], Strand::Forward),
        };

        // Parse enzyme: "Enzyme:SEQ@POS"
        let colon_idx = tag_token
            .find(':')
            .with_context(|| format!("Expected enzyme prefix (Enzyme:SEQ@POS), got: {}", tag_token))?;
        let enzyme_name = &tag_token[..colon_idx];
        let enzyme = parse_enzyme_type(enzyme_name)?;

        // Parse sequence and position: "SEQ@POS"
        let seq_pos_str = &tag_token[colon_idx + 1..];
        let at_idx = seq_pos_str
            .find('@')
            .with_context(|| format!("Expected position (@POS) in tag token, got: {}", tag_token))?;
        let seq_str = &seq_pos_str[..at_idx];
        // Parse sequence and position: "SEQ@POS[:contig_name]"
        let pos_contig_str = &seq_pos_str[at_idx + 1..];
        let (pos_str, contig_name) = if let Some(cidx) = pos_contig_str.find(':') {
            (&pos_contig_str[..cidx], Some(&pos_contig_str[cidx + 1..]))
        } else {
            (pos_contig_str, None)
        };
        let position: u64 = pos_str
            .parse()
            .with_context(|| format!("Invalid position value: {}", tag_token))?;
        let sequence = seq_to_bytes(seq_str)?;

        // Assign contig_id
        let contig_id = if let Some(name) = contig_name {
            if let Some(idx) = record.contig_names.iter().position(|n| n == name) {
                (idx + 1) as u16
            } else {
                record.contig_names.push(name.to_string());
                record.contig_names.len() as u16
            }
        } else {
            0
        };

        let tag = Tag::new(sequence, position, enzyme, strand, contig_id);
        record.add_tag(tag);

        // Validate the pending gap (from previous tag) after add_tag computes the gap
        if let Some(expected) = pending_gap.take() {
            if let Some(last_gap) = record.gaps.last() {
                if last_gap.size != expected {
                    bail!(
                        "Gap mismatch at index {}: computed={}, expected={}",
                        record.gaps.len() - 1, last_gap.size, expected
                    );
                }
            }
        }

        i += 1;
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
        Tag::new(make_seq(seq), pos, enzyme, Strand::Forward, 0)
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
            writeln!(f, "BcgI:ATCG@100 -500- BcgI:GCTA@600").unwrap();
        }

        let mut reader = TgtReader::new(&path).unwrap();
        let record = reader.read_record().unwrap().unwrap();
        assert_eq!(record.genome_id, "G001");
        assert_eq!(record.total_length, 10000);
        assert_eq!(record.tag_count(), 2);
        assert_eq!(record.tags[0].sequence_str(), "ATCG");
        assert_eq!(record.tags[0].enzyme, EnzymeType::BcgI);
        assert_eq!(record.tags[0].position, 100);
        assert_eq!(record.tags[1].sequence_str(), "GCTA");
        assert_eq!(record.tags[1].enzyme, EnzymeType::BcgI);
        assert_eq!(record.tags[1].position, 600);
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
            writeln!(f, "BcgI:ATCG@100 -500- BcgI:GCTA@600").unwrap();
            writeln!(f, ">G002|length=6000").unwrap();
            writeln!(f, "AlfI:TTAA@200 -300- AlfI:CCGG@500").unwrap();
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

    #[test]
    fn test_read_binary_v2_with_contigs() {
        use crate::tgt::writer::TgtWriter;

        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.btgt");

        let mut record = TgtRecord::new("test_genome", 10000);
        record.add_tag(make_tag("AAAA", 100, EnzymeType::BcgI));
        record.add_tag(make_tag("TTTT", 600, EnzymeType::AlfI));
        record.contig_names = vec!["chr1".to_string(), "chr2".to_string()];

        {
            let mut writer = TgtWriter::new(&path).unwrap();
            writer.write_binary(&record).unwrap();
        }

        let mut reader = TgtReader::new(&path).unwrap();
        let read = reader.read_binary().unwrap().unwrap();
        assert_eq!(read.genome_id, "test_genome");
        assert_eq!(read.total_length, 10000);
        assert_eq!(read.tag_count(), 2);
        assert_eq!(read.contig_names, vec!["chr1", "chr2"]);
    }

    #[test]
    fn test_read_binary_rejects_v1() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("v1.btgt");
        {
            let mut f = File::create(&path).unwrap();
            // Write a 48-byte v1-style header so read_exact succeeds and magic check triggers
            f.write_all(b"TGT\x01").unwrap();               // magic (4 bytes)
            f.write_all(&1u32.to_le_bytes()).unwrap();        // version 1 (4 bytes)
            f.write_all(&10000u64.to_le_bytes()).unwrap();    // genome length (8 bytes)
            f.write_all(&0u32.to_le_bytes()).unwrap();        // tag count (4 bytes)
            f.write_all(&0u16.to_le_bytes()).unwrap();        // enzyme count (2 bytes)
            f.write_all(&0u16.to_le_bytes()).unwrap();        // contig count (2 bytes)
            f.write_all(&[0u8; 24]).unwrap();                 // reserved (24 bytes)
            f.write_all(b"TGT\x01").unwrap();
            f.write_all(&1u32.to_le_bytes()).unwrap(); // version 1
            f.write_all(&10000u64.to_le_bytes()).unwrap(); // genome length
            f.write_all(&0u32.to_le_bytes()).unwrap(); // tag count
            f.write_all(&0u16.to_le_bytes()).unwrap(); // enzyme count
            f.write_all(&[0u8; 10]).unwrap(); // reserved
        }

        let mut reader = TgtReader::new(&path).unwrap();
        let result = reader.read_binary();
        assert!(result.is_err());
        let err_msg = format!("{}", result.unwrap_err());
        assert!(
            err_msg.contains("obsolete") || err_msg.contains("v1"),
            "Error should mention obsolete v1 format, got: {}", err_msg
        );
    }
}
