//! TGT (Tag-Gap-Tag) module — core data structures for 2bRAD genome representation

pub mod gap;
pub mod record;
pub mod tag;
pub mod reader;
pub mod writer;

pub use gap::Gap;
pub use record::TgtRecord;
pub use tag::{Strand, Tag};
pub use reader::TgtReader;
pub use writer::TgtWriter;
