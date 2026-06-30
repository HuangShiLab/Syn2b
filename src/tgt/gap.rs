//! Gap data structure — inter-tag distance

use serde::{Deserialize, Serialize};
use std::fmt;

/// Distance in base pairs between adjacent tags
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Gap {
    pub size: u32,
}

impl Gap {
    /// Create a new Gap
    pub fn new(size: u32) -> Self {
        Self { size }
    }
}

impl fmt::Display for Gap {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "-{}-", self.size)
    }
}
