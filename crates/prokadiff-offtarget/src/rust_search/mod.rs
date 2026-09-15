pub mod bulge;
pub mod exact;

pub use bulge::{scan_contig_bulge, scan_genome_bulge, BulgeSearchOptions};
pub use exact::{scan_contig, scan_genome};
