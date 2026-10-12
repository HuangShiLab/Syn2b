//! Syn2b-Cluster: functional gene-cluster ingestion, indexing, and search.

pub mod features;
pub mod index;
pub mod model;
pub mod scoring;
pub mod search;

pub use features::{parse_features, parse_gff, write_features, Feature, FEATURE_COLUMNS};
pub use index::{
    build_cluster_index, collect_cluster_files, ClusterIndex, CLUSTER_INDEX_MAGIC,
    CLUSTER_INDEX_VERSION,
};
pub use model::{
    build_clusters, read_cluster_file, write_cluster_file, GeneCluster, CLUSTER_FILE_MAGIC,
};
pub use scoring::{score_pair, ClusterScores};
pub use search::{
    search_clusters, search_feature_table, write_cluster_hits, ClusterHit, ClusterSearchOptions,
};
