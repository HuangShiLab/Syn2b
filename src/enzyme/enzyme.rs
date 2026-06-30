//! Enzyme type definitions for Type IIB restriction enzymes.
//!
//! Defines the 16 Type IIB restriction enzymes used in 2bSyn for
//! in silico genome digestion, along with their recognition sites
//! and cutting properties.

use serde::{Deserialize, Serialize};
use std::fmt;

/// A Type IIB restriction enzyme type.
///
/// Type IIB enzymes cut on both sides of their recognition site, producing
/// small, defined tag sequences useful for RAD-seq applications.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum EnzymeType {
    BcgI,
    AlfI,
    AloI,
    BaeI,
    BplI,
    BsaXI,
    BslFI,
    Bsp24I,
    CjeI,
    CjePI,
    CspCI,
    FalI,
    HaeIV,
    Hin4I,
    PpiI,
    PsrI,
}

impl EnzymeType {
    /// Return all 16 enzyme types in a static slice.
    pub fn all() -> &'static [EnzymeType] {
        &[
            EnzymeType::BcgI,
            EnzymeType::AlfI,
            EnzymeType::AloI,
            EnzymeType::BaeI,
            EnzymeType::BplI,
            EnzymeType::BsaXI,
            EnzymeType::BslFI,
            EnzymeType::Bsp24I,
            EnzymeType::CjeI,
            EnzymeType::CjePI,
            EnzymeType::CspCI,
            EnzymeType::FalI,
            EnzymeType::HaeIV,
            EnzymeType::Hin4I,
            EnzymeType::PpiI,
            EnzymeType::PsrI,
        ]
    }

    /// Return the 0-based index of this enzyme type (0..16).
    pub fn index(&self) -> u8 {
        match self {
            EnzymeType::BcgI => 0,
            EnzymeType::AlfI => 1,
            EnzymeType::AloI => 2,
            EnzymeType::BaeI => 3,
            EnzymeType::BplI => 4,
            EnzymeType::BsaXI => 5,
            EnzymeType::BslFI => 6,
            EnzymeType::Bsp24I => 7,
            EnzymeType::CjeI => 8,
            EnzymeType::CjePI => 9,
            EnzymeType::CspCI => 10,
            EnzymeType::FalI => 11,
            EnzymeType::HaeIV => 12,
            EnzymeType::Hin4I => 13,
            EnzymeType::PpiI => 14,
            EnzymeType::PsrI => 15,
        }
    }

    /// Convert a 0-based index back to an EnzymeType.
    ///
    /// Returns None if the index is out of range.
    pub fn from_index(index: u8) -> Option<Self> {
        match index {
            0 => Some(EnzymeType::BcgI),
            1 => Some(EnzymeType::AlfI),
            2 => Some(EnzymeType::AloI),
            3 => Some(EnzymeType::BaeI),
            4 => Some(EnzymeType::BplI),
            5 => Some(EnzymeType::BsaXI),
            6 => Some(EnzymeType::BslFI),
            7 => Some(EnzymeType::Bsp24I),
            8 => Some(EnzymeType::CjeI),
            9 => Some(EnzymeType::CjePI),
            10 => Some(EnzymeType::CspCI),
            11 => Some(EnzymeType::FalI),
            12 => Some(EnzymeType::HaeIV),
            13 => Some(EnzymeType::Hin4I),
            14 => Some(EnzymeType::PpiI),
            15 => Some(EnzymeType::PsrI),
            _ => None,
        }
    }
}

impl fmt::Display for EnzymeType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

/// Properties of a restriction enzyme.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Enzyme {
    pub enzyme_type: EnzymeType,
    pub recognition_site: &'static [u8],
    pub tag_length: u8,
    pub cut_offset_5: i8,
    pub cut_offset_3: i8,
}

impl Enzyme {
    /// Return the enzyme properties for a given enzyme type.
    pub fn properties(enzyme_type: EnzymeType) -> Self {
        match enzyme_type {
            EnzymeType::BcgI => Enzyme {
                enzyme_type,
                recognition_site: b"CGANNNNNNTGC",
                tag_length: 32,
                cut_offset_5: -10,
                cut_offset_3: 10,
            },
            EnzymeType::AlfI => Enzyme {
                enzyme_type,
                recognition_site: b"GCANNNNNNTGC",
                tag_length: 32,
                cut_offset_5: -9,
                cut_offset_3: 9,
            },
            EnzymeType::AloI => Enzyme {
                enzyme_type,
                recognition_site: b"GAACNNNNNNTCC",
                tag_length: 32,
                cut_offset_5: -11,
                cut_offset_3: 11,
            },
            EnzymeType::BaeI => Enzyme {
                enzyme_type,
                recognition_site: b"ACNNNNGTAYC",
                tag_length: 32,
                cut_offset_5: -12,
                cut_offset_3: 12,
            },
            EnzymeType::BplI => Enzyme {
                enzyme_type,
                recognition_site: b"GAGNNNNNCTC",
                tag_length: 32,
                cut_offset_5: -8,
                cut_offset_3: 8,
            },
            EnzymeType::BsaXI => Enzyme {
                enzyme_type,
                recognition_site: b"ACNNNNNCTCC",
                tag_length: 32,
                cut_offset_5: -10,
                cut_offset_3: 10,
            },
            EnzymeType::BslFI => Enzyme {
                enzyme_type,
                recognition_site: b"GGGAC",
                tag_length: 28,
                cut_offset_5: -7,
                cut_offset_3: 7,
            },
            EnzymeType::Bsp24I => Enzyme {
                enzyme_type,
                recognition_site: b"GACNNNNNNTGG",
                tag_length: 32,
                cut_offset_5: -11,
                cut_offset_3: 11,
            },
            EnzymeType::CjeI => Enzyme {
                enzyme_type,
                recognition_site: b"RYANNNNNNCTC",
                tag_length: 32,
                cut_offset_5: -10,
                cut_offset_3: 10,
            },
            EnzymeType::CjePI => Enzyme {
                enzyme_type,
                recognition_site: b"GCANNNNNNGTG",
                tag_length: 32,
                cut_offset_5: -10,
                cut_offset_3: 10,
            },
            EnzymeType::CspCI => Enzyme {
                enzyme_type,
                recognition_site: b"CAANNNNNGTGG",
                tag_length: 32,
                cut_offset_5: -11,
                cut_offset_3: 11,
            },
            EnzymeType::FalI => Enzyme {
                enzyme_type,
                recognition_site: b"AAGNNNNNCTT",
                tag_length: 32,
                cut_offset_5: -9,
                cut_offset_3: 9,
            },
            EnzymeType::HaeIV => Enzyme {
                enzyme_type,
                recognition_site: b"GAYNNNNNRTC",
                tag_length: 32,
                cut_offset_5: -9,
                cut_offset_3: 9,
            },
            EnzymeType::Hin4I => Enzyme {
                enzyme_type,
                recognition_site: b"GAYNNNNNVTC",
                tag_length: 32,
                cut_offset_5: -10,
                cut_offset_3: 10,
            },
            EnzymeType::PpiI => Enzyme {
                enzyme_type,
                recognition_site: b"GAACNNNNNCTC",
                tag_length: 32,
                cut_offset_5: -10,
                cut_offset_3: 10,
            },
            EnzymeType::PsrI => Enzyme {
                enzyme_type,
                recognition_site: b"GAACNNNNNNTAC",
                tag_length: 32,
                cut_offset_5: -11,
                cut_offset_3: 11,
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_all_enzymes_count() {
        assert_eq!(EnzymeType::all().len(), 16);
    }

    #[test]
    fn test_enzyme_index_roundtrip() {
        for enzyme in EnzymeType::all() {
            let idx = enzyme.index();
            assert_eq!(EnzymeType::from_index(idx), Some(*enzyme));
        }
    }

    #[test]
    fn test_enzyme_index_out_of_range() {
        assert_eq!(EnzymeType::from_index(16), None);
        assert_eq!(EnzymeType::from_index(255), None);
    }

    #[test]
    fn test_enzyme_display() {
        assert_eq!(format!("{}", EnzymeType::BcgI), "BcgI");
        assert_eq!(format!("{}", EnzymeType::PsrI), "PsrI");
    }

    #[test]
    fn test_enzyme_properties() {
        let bcgi = Enzyme::properties(EnzymeType::BcgI);
        assert_eq!(bcgi.enzyme_type, EnzymeType::BcgI);
        assert_eq!(bcgi.tag_length, 32);
    }
}
