//! Enzyme type definitions for Type IIB restriction enzymes.
//!
//! Defines the 16 Type IIB restriction enzymes used in Syn2b for
//! in silico genome digestion. Based on the 2bRAD-M paper and
//! Fast2bRAD-M implementation.

use serde::{Deserialize, Serialize};
use std::fmt;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum EnzymeType {
    BcgI, AlfI, AloI, BaeI, BplI, BsaXI, BslFI, Bsp24I,
    CjeI, CjePI, CspCI, FalI, HaeIV, Hin4I, PpiI, PsrI,
    /// Not a restriction enzyme — the marker for a landmark selected by
    /// FracMinHash rather than by digestion. It lives in this enum because the
    /// per-tag `enzyme` field is what already travels through both TGT formats,
    /// so recording the landmark source here needs no format change and no
    /// record-level flag that could disagree with the tags it describes.
    ///
    /// It is deliberately excluded from `all()` (which means "all enzymes") and
    /// rejected by the enzyme-name parsers, so it can never be requested as one.
    /// `properties()` returns a zero-length, pattern-free descriptor: there is no
    /// recognition site to match, and `digest_genome_contig` must never be called
    /// with it.
    FracMinHash,
}

impl EnzymeType {
    pub fn all() -> &'static [EnzymeType] {
        &[
            EnzymeType::BcgI, EnzymeType::AlfI, EnzymeType::AloI, EnzymeType::BaeI,
            EnzymeType::BplI, EnzymeType::BsaXI, EnzymeType::BslFI, EnzymeType::Bsp24I,
            EnzymeType::CjeI, EnzymeType::CjePI, EnzymeType::CspCI, EnzymeType::FalI,
            EnzymeType::HaeIV, EnzymeType::Hin4I, EnzymeType::PpiI, EnzymeType::PsrI,
        ]
    }
    pub fn properties(&self) -> Enzyme { Enzyme::properties(*self) }
    pub fn index(&self) -> u8 {
        match self {
            EnzymeType::BcgI=>0, EnzymeType::AlfI=>1, EnzymeType::AloI=>2, EnzymeType::BaeI=>3,
            EnzymeType::BplI=>4, EnzymeType::BsaXI=>5, EnzymeType::BslFI=>6, EnzymeType::Bsp24I=>7,
            EnzymeType::CjeI=>8, EnzymeType::CjePI=>9, EnzymeType::CspCI=>10, EnzymeType::FalI=>11,
            EnzymeType::HaeIV=>12, EnzymeType::Hin4I=>13, EnzymeType::PpiI=>14, EnzymeType::PsrI=>15,
            EnzymeType::FracMinHash=>16,
        }
    }
    pub fn from_index(i: u8) -> Option<Self> {
        match i {
            0=>Some(EnzymeType::BcgI), 1=>Some(EnzymeType::AlfI), 2=>Some(EnzymeType::AloI), 3=>Some(EnzymeType::BaeI),
            4=>Some(EnzymeType::BplI), 5=>Some(EnzymeType::BsaXI), 6=>Some(EnzymeType::BslFI), 7=>Some(EnzymeType::Bsp24I),
            8=>Some(EnzymeType::CjeI), 9=>Some(EnzymeType::CjePI), 10=>Some(EnzymeType::CspCI), 11=>Some(EnzymeType::FalI),
            12=>Some(EnzymeType::HaeIV), 13=>Some(EnzymeType::Hin4I), 14=>Some(EnzymeType::PpiI), 15=>Some(EnzymeType::PsrI),
            16=>Some(EnzymeType::FracMinHash),
            _=>None,
        }
    }
}

impl fmt::Display for EnzymeType { fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result { write!(f, "{:?}", self) } }

// ── Anchor / IUPAC / Pattern / Enzyme ──────────────────────────────────────

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Anchor { pub offset: usize, pub motif: &'static [u8] }

/// IUPAC degenerate base constraint.
/// `allowed` is a bitmask: bit0=A, bit1=T, bit2=C, bit3=G.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct IupacConstraint { pub offset: usize, pub allowed: u8 }

const BASE_MASK: [u8; 256] = {
    let mut t = [0u8; 256];
    t[b'A' as usize]=1; t[b'a' as usize]=1; t[b'T' as usize]=2; t[b't' as usize]=2;
    t[b'C' as usize]=4; t[b'c' as usize]=4; t[b'G' as usize]=8; t[b'g' as usize]=8;
    t
};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Pattern { pub anchors: &'static [Anchor], pub iupac: &'static [IupacConstraint] }

/// For `EnzymeType::FracMinHash`, which recognises nothing.
static NO_PATTERNS: [Pattern; 0] = [];

impl Pattern {
    pub fn matches(&self, window: &[u8]) -> bool {
        let anchors_ok = self.anchors.iter().all(|a| {
            let e = a.offset + a.motif.len();
            e <= window.len() && &window[a.offset..e] == a.motif
        });
        if !anchors_ok { return false; }
        self.iupac.iter().all(|c| c.offset < window.len() && (BASE_MASK[window[c.offset] as usize] & c.allowed) != 0)
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Enzyme { pub enzyme_type: EnzymeType, pub tag_length: u8, pub patterns: &'static [Pattern] }

impl Enzyme {
    pub fn properties(t: EnzymeType) -> Self {
        match t {
            EnzymeType::BcgI => Enzyme { enzyme_type:t, tag_length:32, patterns:&BCGI_PATTERNS },
            EnzymeType::AlfI => Enzyme { enzyme_type:t, tag_length:32, patterns:&ALFI_PATTERNS },
            EnzymeType::AloI => Enzyme { enzyme_type:t, tag_length:27, patterns:&ALOI_PATTERNS },
            EnzymeType::BaeI => Enzyme { enzyme_type:t, tag_length:28, patterns:&BAEI_PATTERNS },
            EnzymeType::BplI => Enzyme { enzyme_type:t, tag_length:27, patterns:&BPLI_PATTERNS },
            EnzymeType::BsaXI => Enzyme { enzyme_type:t, tag_length:27, patterns:&BSAXI_PATTERNS },
            EnzymeType::BslFI => Enzyme { enzyme_type:t, tag_length:25, patterns:&BSLFI_PATTERNS },
            EnzymeType::Bsp24I => Enzyme { enzyme_type:t, tag_length:27, patterns:&BSP24I_PATTERNS },
            EnzymeType::CjeI => Enzyme { enzyme_type:t, tag_length:28, patterns:&CJEI_PATTERNS },
            EnzymeType::CjePI => Enzyme { enzyme_type:t, tag_length:27, patterns:&CJEPI_PATTERNS },
            EnzymeType::CspCI => Enzyme { enzyme_type:t, tag_length:33, patterns:&CSPCI_PATTERNS },
            EnzymeType::FalI => Enzyme { enzyme_type:t, tag_length:27, patterns:&FALI_PATTERNS },
            EnzymeType::HaeIV => Enzyme { enzyme_type:t, tag_length:27, patterns:&HAEIV_PATTERNS },
            EnzymeType::Hin4I => Enzyme { enzyme_type:t, tag_length:27, patterns:&HIN4I_PATTERNS },
            EnzymeType::PpiI => Enzyme { enzyme_type:t, tag_length:27, patterns:&PPII_PATTERNS },
            EnzymeType::PsrI => Enzyme { enzyme_type:t, tag_length:27, patterns:&PSRI_PATTERNS },
            // No recognition site: FracMinHash selects on a hash of the k-mer, not
            // on sequence content. `tag_length` is set by the sketch's k, not here.
            EnzymeType::FracMinHash => Enzyme { enzyme_type:t, tag_length:0, patterns:&NO_PATTERNS },
        }
    }
}

// ── 16 Enzyme definitions ─────────────────────────────────────────────────

// 1. BcgI (32)
const BCGI_F1: Anchor = Anchor{offset:10, motif:b"CGA"};
const BCGI_F2: Anchor = Anchor{offset:19, motif:b"TGC"};
const BCGI_R1: Anchor = Anchor{offset:10, motif:b"GCA"};
const BCGI_R2: Anchor = Anchor{offset:19, motif:b"TCG"};
const BCGI_PATTERNS: [Pattern;2] = [
    Pattern{anchors:&[BCGI_F1,BCGI_F2], iupac:&[]},
    Pattern{anchors:&[BCGI_R1,BCGI_R2], iupac:&[]},
];

// 2. AlfI (32, palindrome)
const ALFI_A1: Anchor = Anchor{offset:10, motif:b"GCA"};
const ALFI_A2: Anchor = Anchor{offset:19, motif:b"TGC"};
const ALFI_PATTERNS: [Pattern;1] = [Pattern{anchors:&[ALFI_A1,ALFI_A2], iupac:&[]}];

// 3. AloI (27)
const ALOI_F1: Anchor = Anchor{offset:7, motif:b"GAAC"};
const ALOI_F2: Anchor = Anchor{offset:17, motif:b"TCC"};
const ALOI_R1: Anchor = Anchor{offset:7, motif:b"GGA"};
const ALOI_R2: Anchor = Anchor{offset:16, motif:b"GTTC"};
const ALOI_PATTERNS: [Pattern;2] = [
    Pattern{anchors:&[ALOI_F1,ALOI_F2], iupac:&[]},
    Pattern{anchors:&[ALOI_R1,ALOI_R2], iupac:&[]},
];

// 4. BaeI (28, degenerate)
// fwd [CT]@19 (Y=6)  rev [AG]@8 (R=9)
const BAEI_F1: Anchor = Anchor{offset:10, motif:b"AC"};
const BAEI_F2: Anchor = Anchor{offset:16, motif:b"GTA"};
const BAEI_R1: Anchor = Anchor{offset:7, motif:b"G"};
const BAEI_R2: Anchor = Anchor{offset:9, motif:b"TAC"};
const BAEI_FWD_IUPAC: [IupacConstraint;1] = [IupacConstraint{offset:19, allowed:6}]; // Y=[CT]
const BAEI_REV_IUPAC: [IupacConstraint;1] = [IupacConstraint{offset:8, allowed:9}];  // R=[AG]
const BAEI_PATTERNS: [Pattern;2] = [
    Pattern{anchors:&[BAEI_F1,BAEI_F2], iupac:&BAEI_FWD_IUPAC},
    Pattern{anchors:&[BAEI_R1,BAEI_R2], iupac:&BAEI_REV_IUPAC},
];

// 5. BplI (27, palindrome)
const BPLI_A1: Anchor = Anchor{offset:8, motif:b"GAG"};
const BPLI_A2: Anchor = Anchor{offset:16, motif:b"CTC"};
const BPLI_PATTERNS: [Pattern;1] = [Pattern{anchors:&[BPLI_A1,BPLI_A2], iupac:&[]}];

// 6. BsaXI (27)
const BSAXI_F1: Anchor = Anchor{offset:9, motif:b"AC"};
const BSAXI_F2: Anchor = Anchor{offset:16, motif:b"CTCC"};
const BSAXI_R1: Anchor = Anchor{offset:7, motif:b"GGAG"};
const BSAXI_R2: Anchor = Anchor{offset:16, motif:b"GT"};
const BSAXI_PATTERNS: [Pattern;2] = [
    Pattern{anchors:&[BSAXI_F1,BSAXI_F2], iupac:&[]},
    Pattern{anchors:&[BSAXI_R1,BSAXI_R2], iupac:&[]},
];

// 7. BslFI (25)
const BSLFI_F1: Anchor = Anchor{offset:6, motif:b"GGGAC"};
const BSLFI_R1: Anchor = Anchor{offset:14, motif:b"GTCCC"};
const BSLFI_PATTERNS: [Pattern;2] = [
    Pattern{anchors:&[BSLFI_F1], iupac:&[]},
    Pattern{anchors:&[BSLFI_R1], iupac:&[]},
];

// 8. Bsp24I (27)
const BSP24I_F1: Anchor = Anchor{offset:8, motif:b"GAC"};
const BSP24I_F2: Anchor = Anchor{offset:17, motif:b"TGG"};
const BSP24I_R1: Anchor = Anchor{offset:7, motif:b"CCA"};
const BSP24I_R2: Anchor = Anchor{offset:16, motif:b"GTC"};
const BSP24I_PATTERNS: [Pattern;2] = [
    Pattern{anchors:&[BSP24I_F1,BSP24I_F2], iupac:&[]},
    Pattern{anchors:&[BSP24I_R1,BSP24I_R2], iupac:&[]},
];

// 9. CjeI (28)
const CJEI_F1: Anchor = Anchor{offset:8, motif:b"CCA"};
const CJEI_F2: Anchor = Anchor{offset:17, motif:b"GT"};
const CJEI_R1: Anchor = Anchor{offset:9, motif:b"AC"};
const CJEI_R2: Anchor = Anchor{offset:17, motif:b"TGG"};
const CJEI_PATTERNS: [Pattern;2] = [
    Pattern{anchors:&[CJEI_F1,CJEI_F2], iupac:&[]},
    Pattern{anchors:&[CJEI_R1,CJEI_R2], iupac:&[]},
];

// 10. CjePI (27)
const CJEPI_F1: Anchor = Anchor{offset:7, motif:b"CCA"};
const CJEPI_F2: Anchor = Anchor{offset:17, motif:b"TC"};
const CJEPI_R1: Anchor = Anchor{offset:8, motif:b"GA"};
const CJEPI_R2: Anchor = Anchor{offset:17, motif:b"TGG"};
const CJEPI_PATTERNS: [Pattern;2] = [
    Pattern{anchors:&[CJEPI_F1,CJEPI_F2], iupac:&[]},
    Pattern{anchors:&[CJEPI_R1,CJEPI_R2], iupac:&[]},
];

// 11. CspCI (33)
const CSPCI_F1: Anchor = Anchor{offset:11, motif:b"CAA"};
const CSPCI_F2: Anchor = Anchor{offset:19, motif:b"GTGG"};
const CSPCI_R1: Anchor = Anchor{offset:10, motif:b"CCAC"};
const CSPCI_R2: Anchor = Anchor{offset:19, motif:b"TTG"};
const CSPCI_PATTERNS: [Pattern;2] = [
    Pattern{anchors:&[CSPCI_F1,CSPCI_F2], iupac:&[]},
    Pattern{anchors:&[CSPCI_R1,CSPCI_R2], iupac:&[]},
];

// 12. FalI (27, palindrome)
const FALI_A1: Anchor = Anchor{offset:8, motif:b"AAG"};
const FALI_A2: Anchor = Anchor{offset:16, motif:b"CTT"};
const FALI_PATTERNS: [Pattern;1] = [Pattern{anchors:&[FALI_A1,FALI_A2], iupac:&[]}];

// 13. HaeIV (27, degenerate)
// fwd Y@9 R@15  rev Y@11 R@17
const HAEIV_F1: Anchor = Anchor{offset:7, motif:b"GA"};
const HAEIV_R1: Anchor = Anchor{offset:9, motif:b"GA"};
const HAEIV_FWD_IUPAC: [IupacConstraint;2] = [
    IupacConstraint{offset:9, allowed:6},   // Y=[CT]
    IupacConstraint{offset:15, allowed:9},  // R=[AG]
];
const HAEIV_REV_IUPAC: [IupacConstraint;2] = [
    IupacConstraint{offset:11, allowed:6},  // Y=[CT]
    IupacConstraint{offset:17, allowed:9},  // R=[AG]
];
const HAEIV_PATTERNS: [Pattern;2] = [
    Pattern{anchors:&[HAEIV_F1], iupac:&HAEIV_FWD_IUPAC},
    Pattern{anchors:&[HAEIV_R1], iupac:&HAEIV_REV_IUPAC},
];

// 14. Hin4I (27, degenerate)
// fwd Y@10 [GAC]@16  rev [CTG]@10 R@16
const HIN4I_F1: Anchor = Anchor{offset:8, motif:b"GA"};
const HIN4I_R1: Anchor = Anchor{offset:8, motif:b"GA"};
const HIN4I_FWD_IUPAC: [IupacConstraint;2] = [
    IupacConstraint{offset:10, allowed:6},   // Y=[CT]
    IupacConstraint{offset:16, allowed:13},  // [GAC]
];
const HIN4I_REV_IUPAC: [IupacConstraint;2] = [
    IupacConstraint{offset:10, allowed:14},  // [CTG]
    IupacConstraint{offset:16, allowed:9},   // R=[AG]
];
const HIN4I_PATTERNS: [Pattern;2] = [
    Pattern{anchors:&[HIN4I_F1], iupac:&HIN4I_FWD_IUPAC},
    Pattern{anchors:&[HIN4I_R1], iupac:&HIN4I_REV_IUPAC},
];

// 15. PpiI (27)
const PPII_F1: Anchor = Anchor{offset:7, motif:b"GAAC"};
const PPII_F2: Anchor = Anchor{offset:16, motif:b"CTC"};
const PPII_R1: Anchor = Anchor{offset:8, motif:b"GAG"};
const PPII_R2: Anchor = Anchor{offset:16, motif:b"GTTC"};
const PPII_PATTERNS: [Pattern;2] = [
    Pattern{anchors:&[PPII_F1,PPII_F2], iupac:&[]},
    Pattern{anchors:&[PPII_R1,PPII_R2], iupac:&[]},
];

// 16. PsrI (27)
const PSRI_F1: Anchor = Anchor{offset:7, motif:b"GAAC"};
const PSRI_F2: Anchor = Anchor{offset:17, motif:b"TAC"};
const PSRI_R1: Anchor = Anchor{offset:7, motif:b"GTA"};
const PSRI_R2: Anchor = Anchor{offset:16, motif:b"GTTC"};
const PSRI_PATTERNS: [Pattern;2] = [
    Pattern{anchors:&[PSRI_F1,PSRI_F2], iupac:&[]},
    Pattern{anchors:&[PSRI_R1,PSRI_R2], iupac:&[]},
];

// ── Tests ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test] fn test_all_count() { assert_eq!(EnzymeType::all().len(), 16); }

    #[test] fn test_index_roundtrip() {
        for e in EnzymeType::all() { assert_eq!(EnzymeType::from_index(e.index()), Some(*e)); }
    }

    #[test] fn test_tag_lengths() {
        assert_eq!(Enzyme::properties(EnzymeType::BcgI).tag_length, 32);
        assert_eq!(Enzyme::properties(EnzymeType::AlfI).tag_length, 32);
        assert_eq!(Enzyme::properties(EnzymeType::AloI).tag_length, 27);
        assert_eq!(Enzyme::properties(EnzymeType::BaeI).tag_length, 28);
        assert_eq!(Enzyme::properties(EnzymeType::BplI).tag_length, 27);
        assert_eq!(Enzyme::properties(EnzymeType::BsaXI).tag_length, 27);
        assert_eq!(Enzyme::properties(EnzymeType::BslFI).tag_length, 25);
        assert_eq!(Enzyme::properties(EnzymeType::Bsp24I).tag_length, 27);
        assert_eq!(Enzyme::properties(EnzymeType::CjeI).tag_length, 28);
        assert_eq!(Enzyme::properties(EnzymeType::CjePI).tag_length, 27);
        assert_eq!(Enzyme::properties(EnzymeType::CspCI).tag_length, 33);
        assert_eq!(Enzyme::properties(EnzymeType::FalI).tag_length, 27);
        assert_eq!(Enzyme::properties(EnzymeType::HaeIV).tag_length, 27);
        assert_eq!(Enzyme::properties(EnzymeType::Hin4I).tag_length, 27);
        assert_eq!(Enzyme::properties(EnzymeType::PpiI).tag_length, 27);
        assert_eq!(Enzyme::properties(EnzymeType::PsrI).tag_length, 27);
    }

    #[test] fn test_baeI_iupac() {
        // BaeI fwd: offset 19 must be C or T (Y=bitmask 6)
        let pat = &BAEI_PATTERNS[0];
        assert!(pat.matches(b"AAAAAAAAAAACAAAAGTACCAAAAAAA"));  // C@19
        assert!(pat.matches(b"AAAAAAAAAAACAAAAGTACTAAAAAAA"));  // T@19
        assert!(!pat.matches(b"AAAAAAAAAAACAAAAGTAACAAAAAAA")); // A@19
        assert!(!pat.matches(b"AAAAAAAAAAACAAAAGTAGCAAAAAAA")); // G@19
    }

    #[test] fn test_baeI_rev_iupac() {
        // BaeI rev: anchor1 @7=G, anchor2 @9=TAC, iupac @8=R=[AG]=9
        let pat = &BAEI_PATTERNS[1];
        assert!(pat.matches(b"AAAAAAAGATACAAAAAAAAAAAAAAAA"));  // A@8
        assert!(pat.matches(b"AAAAAAAGGTACAAAAAAAAAAAAAAAA"));  // G@8
        assert!(!pat.matches(b"AAAAAAAGCTACAAAAAAAAAAAAAAAA")); // C@8
        assert!(!pat.matches(b"AAAAAAAGTTACAAAAAAAAAAAAAAAA")); // T@8
    }

    #[test] fn test_haeIV_iupac() {
        // HaeIV fwd: anchor @7=GA, iupac @9=Y=[CT], @15=R=[AG]
        let pat = &HAEIV_PATTERNS[0];
        assert!(pat.matches(b"AAAAAAAGACTAAAAGATCAAAAAAAAAA")); // C@9, A@15
        assert!(pat.matches(b"AAAAAAAGATTAAAAGCTCAAAAAAAAAA")); // T@9, G@15
        assert!(!pat.matches(b"AAAAAAAGAATAAAAGATCAAAAAAAAAA")); // A@9 (bad)
        assert!(!pat.matches(b"AAAAAAAGACTAAAATTCAAAAAAAAAA"));  // T@15 (bad)
    }

    #[test] fn test_haeIV_rev_iupac() {
        // HaeIV rev: anchor @9=GA, iupac @11=Y=[CT], @17=R=[AG]
        let pat = &HAEIV_PATTERNS[1];
        assert!(pat.matches(b"AAAAAAAAAGACTAAAAAGAAAAAAAAA"));  // C@11, G@17
        assert!(pat.matches(b"AAAAAAAAAGATTAAAAAGAAAAAAAAA"));  // T@11, G@17
        assert!(!pat.matches(b"AAAAAAAAAGAATAAAAAGAAAAAAAAA")); // A@11 (bad)
        assert!(!pat.matches(b"AAAAAAAAAGACTAAAATAAAAAAAAAA")); // T@17 (bad)
    }

    #[test] fn test_hin4I_iupac() {
        // Hin4I fwd: anchor @8=GA, iupac @10=Y=[CT], @16=[GAC]=13
        let pat = &HIN4I_PATTERNS[0];
        assert!(pat.matches(b"AAAAAAAAGACAAAAAGAAAAAAAAAAA"));  // C@10, G@16
        assert!(pat.matches(b"AAAAAAAAGATAAAAAAAAAAAAAAAAA"));  // T@10, A@16
        assert!(!pat.matches(b"AAAAAAAAGAAAAAAAAGAAAAAAAAAA")); // A@10 (bad)
        assert!(!pat.matches(b"AAAAAAAAGACAAAAATAAAAAAAAAAA")); // T@16 (bad)
    }

    #[test] fn test_hin4I_rev_iupac() {
        // Hin4I rev: anchor @8=GA, iupac @10=[CTG]=14, @16=R=[AG]=9
        let pat = &HIN4I_PATTERNS[1];
        assert!(pat.matches(b"AAAAAAAAGACAAAAAGAAAAAAAAAAA"));  // C@10, G@16
        assert!(pat.matches(b"AAAAAAAAGATAAAAAGAAAAAAAAAAA"));  // T@10, G@16
        assert!(!pat.matches(b"AAAAAAAAGAAAAAAAAGAAAAAAAAAA")); // A@10 (bad)
        assert!(!pat.matches(b"AAAAAAAAGACAAAAATAAAAAAAAAAA")); // T@16 (bad)
    }

    #[test] fn test_ppiI() {
        let ppi = Enzyme::properties(EnzymeType::PpiI);
        assert_eq!(ppi.tag_length, 27);
        let pat = &PPII_PATTERNS[0];
        let seq = b"AAAAAAAGAACAAAAACTCAAAAAAAA"; // 27 bp
        assert!(pat.matches(seq));
    }

    #[test] fn test_bcgI_pattern() {
        let seq = b"AAAAAAAAAACGAAAAAAATGCAAAAAAAA"; // 32 bp
        assert!(BCGI_PATTERNS[0].matches(seq));
        assert!(!BCGI_PATTERNS[1].matches(seq));
    }
}
