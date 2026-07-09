//! TGT writer — text and binary output for TGT records.
//!
//! Provides `TgtWriter` for writing `TgtRecord` instances in either
//! human-readable text format or compact binary format.

use crate::tgt::record::TgtRecord;
use anyhow::{Context, Result};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

/// Writer for TGT records (text or binary format).
pub struct TgtWriter {
    writer: BufWriter<File>,
}

impl TgtWriter {
    /// Create a new TGT writer that writes to the given file path.
    ///
    /// The file is created if it doesn't exist, or truncated if it does.
    pub fn new(path: &Path) -> Result<Self> {
        let file = File::create(path)
            .with_context(|| format!("Failed to create TGT output file: {}", path.display()))?;
        Ok(Self {
            writer: BufWriter::new(file),
        })
    }

    /// Write a single TGT record in text format.
    ///
    /// Uses the `Display` implementation of `TgtRecord`, which produces:
    /// ```text
    /// >genome_id|length=4641652
    /// BcgI:ATCG... -1313- GCTA... -1298- TTAA...
    /// ```
    pub fn write_record(&mut self, record: &TgtRecord) -> Result<()> {
        writeln!(self.writer, "{}", record)
            .context("Failed to write TGT record in text format")?;
        Ok(())
    }

    /// Write a single TGT record in binary format.
    ///
    /// Binary format layout:
    /// - Header (32 bytes): magic "TGT\x01", version, genome length, tag count, enzyme count
    /// - Tag table (N x 48 bytes each): sequence, position, enzyme index, strand
    /// - Gap table ((N-1) x 4 bytes each): gap sizes
    pub fn write_binary(&mut self, record: &TgtRecord) -> Result<()> {
        // --- Header (32 bytes) ---
        let magic = b"TGT\x01";
        self.writer.write_all(magic)?; // bytes 0..4

        let version: u32 = 1;
        self.writer.write_all(&version.to_le_bytes())?; // bytes 4..8

        self.writer
            .write_all(&record.total_length.to_le_bytes())?; // bytes 8..16

        let tag_count = record.tags.len() as u32;
        self.writer.write_all(&tag_count.to_le_bytes())?; // bytes 16..20

        let enzyme_count = record.enzyme_count() as u16;
        self.writer
            .write_all(&enzyme_count.to_le_bytes())?; // bytes 20..22

        // Reserved bytes 22..32 (10 bytes)
        self.writer.write_all(&[0u8; 10])?;

        // --- Genome ID (variable length) ---
        let genome_id_bytes = record.genome_id.as_bytes();
        let id_len = genome_id_bytes.len() as u16;
        self.writer.write_all(&id_len.to_le_bytes())?;
        self.writer.write_all(genome_id_bytes)?;

        // --- Tag Table (N x 48 bytes) ---
        for tag in &record.tags {
            // Sequence: 32 bytes
            self.writer.write_all(&tag.sequence)?; // bytes 0..32

            // Position: 8 bytes (u64)
            self.writer.write_all(&tag.position.to_le_bytes())?; // bytes 32..40

            // Enzyme index: 1 byte (u8)
            self.writer.write_all(&[tag.enzyme.index()])?; // byte 40

            // Strand: 1 byte (0=FWD, 1=REV)
            self.writer.write_all(&[tag.strand.to_u8()])?; // byte 41

            // Contig ID: 2 bytes (u16)
            self.writer.write_all(&tag.contig_id.to_le_bytes())?; // bytes 42..44

            // Reserved: 4 bytes
            self.writer.write_all(&[0u8; 4])?; // bytes 44..48
        }

        // --- Gap Table ((N-1) x 4 bytes) ---
        for gap in &record.gaps {
            self.writer.write_all(&gap.size.to_le_bytes())?;
        }

        self.writer.flush().context("Failed to flush binary TGT data")?;
        Ok(())
    }

    /// Flush any buffered data to the underlying file.
    pub fn flush(&mut self) -> Result<()> {
        self.writer.flush().context("Failed to flush TGT writer")?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::enzyme::enzyme::EnzymeType;
    use crate::tgt::tag::{Strand, Tag};
    use crate::tgt::record::TgtRecord;
    use std::io::Read;

    fn make_seq(s: &str) -> [u8; 32] {
        let mut arr = [0u8; 32];
        let bytes = s.as_bytes();
        arr[..bytes.len()].copy_from_slice(bytes);
        arr
    }

    fn make_tag(seq: &str, pos: u64, enzyme: EnzymeType) -> Tag {
        Tag::new(make_seq(seq), pos, enzyme, Strand::Forward, 0)
    }

    fn make_test_record() -> TgtRecord {
        let mut record = TgtRecord::new("NC_000913", 4_641_652);
        record.add_tag(make_tag("ATCGATCGATCGATCG", 100, EnzymeType::BcgI));
        record.add_tag(make_tag("GCTAGCTAGCTAGCTA", 1413, EnzymeType::BcgI));
        record.add_tag(make_tag("CGATCGATCGATCGAT", 2711, EnzymeType::BcgI));
        record
    }

    #[test]
    fn test_write_text_record() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.tgt");
        let record = make_test_record();

        {
            let mut writer = TgtWriter::new(&path).unwrap();
            writer.write_record(&record).unwrap();
        }

        let mut contents = String::new();
        File::open(&path).unwrap().read_to_string(&mut contents).unwrap();

        assert!(contents.starts_with(">NC_000913|length=4641652"));
        assert!(contents.contains("BcgI:ATCGATCGATCGATCG"));
        assert!(contents.contains("-1313-"));
        assert!(contents.contains("-1298-"));
        assert!(contents.contains("BcgI:GCTAGCTAGCTAGCTA"));
        assert!(contents.contains("BcgI:CGATCGATCGATCGAT"));
    }

    #[test]
    fn test_write_binary_header() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.btgt");
        let record = make_test_record();

        {
            let mut writer = TgtWriter::new(&path).unwrap();
            writer.write_binary(&record).unwrap();
        }

        let mut buf = Vec::new();
        File::open(&path).unwrap().read_to_end(&mut buf).unwrap();

        // Check magic
        assert_eq!(&buf[0..4], b"TGT\x01");

        // Check version
        assert_eq!(u32::from_le_bytes([buf[4], buf[5], buf[6], buf[7]]), 1);

        // Check genome length
        assert_eq!(
            u64::from_le_bytes([
                buf[8], buf[9], buf[10], buf[11], buf[12], buf[13], buf[14], buf[15]
            ]),
            4_641_652
        );

        // Check tag count (3)
        assert_eq!(u32::from_le_bytes([buf[16], buf[17], buf[18], buf[19]]), 3);

        // Check enzyme count (1)
        assert_eq!(u16::from_le_bytes([buf[20], buf[21]]), 1);
    }

    #[test]
    fn test_write_binary_tag_table() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.btgt");
        let record = make_test_record();

        {
            let mut writer = TgtWriter::new(&path).unwrap();
            writer.write_binary(&record).unwrap();
        }

        let mut buf = Vec::new();
        File::open(&path).unwrap().read_to_end(&mut buf).unwrap();

        // Tag table starts after header (32) + genome_id length (2) + genome_id (9)
        let tag0_offset = 32 + 2 + 9;

        // Check first tag sequence (first 4 bytes should be "ATCG")
        assert_eq!(&buf[tag0_offset..tag0_offset + 4], b"ATCG");

        // Check first tag position (byte offset 32+32=64)
        let pos_offset = tag0_offset + 32;
        let pos = u64::from_le_bytes([
            buf[pos_offset],
            buf[pos_offset + 1],
            buf[pos_offset + 2],
            buf[pos_offset + 3],
            buf[pos_offset + 4],
            buf[pos_offset + 5],
            buf[pos_offset + 6],
            buf[pos_offset + 7],
        ]);
        assert_eq!(pos, 100);

        // Check enzyme index (byte offset 32+40=72)
        assert_eq!(buf[tag0_offset + 40], 0); // BcgI = index 0

        // Check strand (byte offset 32+41=73)
        assert_eq!(buf[tag0_offset + 41], 0); // Forward
    }

    #[test]
    fn test_write_binary_gap_table() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("test.btgt");
        let record = make_test_record();

        {
            let mut writer = TgtWriter::new(&path).unwrap();
            writer.write_binary(&record).unwrap();
        }

        let mut buf = Vec::new();
        File::open(&path).unwrap().read_to_end(&mut buf).unwrap();

        // Gap table starts after header (32) + genome_id (2+9) + tag table (3*48 = 144) = 187
        let gap_offset = 32 + 2 + 9 + 3 * 48;

        // Check first gap (1313)
        let gap0 = u32::from_le_bytes([buf[gap_offset], buf[gap_offset + 1], buf[gap_offset + 2], buf[gap_offset + 3]]);
        assert_eq!(gap0, 1313);

        // Check second gap (1298)
        let gap1 = u32::from_le_bytes([buf[gap_offset + 4], buf[gap_offset + 5], buf[gap_offset + 6], buf[gap_offset + 7]]);
        assert_eq!(gap1, 1298);
    }

    #[test]
    fn test_write_multiple_records() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("multi.tgt");

        let mut record1 = TgtRecord::new("G001", 10000);
        record1.add_tag(make_tag("AAAA", 100, EnzymeType::BcgI));
        record1.add_tag(make_tag("TTTT", 600, EnzymeType::BcgI));

        let mut record2 = TgtRecord::new("G002", 12000);
        record2.add_tag(make_tag("CCCC", 200, EnzymeType::AlfI));
        record2.add_tag(make_tag("GGGG", 900, EnzymeType::AlfI));

        {
            let mut writer = TgtWriter::new(&path).unwrap();
            writer.write_record(&record1).unwrap();
            writer.write_record(&record2).unwrap();
        }

        let mut contents = String::new();
        File::open(&path).unwrap().read_to_string(&mut contents).unwrap();

        assert!(contents.contains(">G001|length=10000"));
        assert!(contents.contains(">G002|length=12000"));
        assert!(contents.contains("BcgI:"));
        assert!(contents.contains("AlfI:"));
    }
}
