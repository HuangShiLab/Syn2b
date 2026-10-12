//! Deterministic greedy/DP chaining for ordered landmark matches.

use serde::{Deserialize, Serialize};

/// One shared canonical landmark between a query contig and a reference contig.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct TagMatch {
    pub query_position: u64,
    pub reference_position: u64,
    pub tag_id: u32,
    /// Retained per match so a genome-level chain can still report the best
    /// reference contig and count contig-boundary transitions.
    pub reference_contig_id: u32,
}

/// A maximal ordered run of landmarks in one orientation.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct ChainedBlock {
    pub query_start: u64,
    pub query_end: u64,
    pub reference_start: u64,
    pub reference_end: u64,
    pub orientation: bool,
    pub reference_contig_id: u32,
    pub n_tags: u64,
    pub n_adjacencies: u64,
}

/// Weighted longest ordered chain, evaluated independently in each orientation.
///
/// Input need not be sorted. Repeated reference landmarks are allowed; the DP
/// will not join a repeated position into a false increasing run. The tie break
/// is explicit and deterministic, preferring the lexicographically smaller
/// predecessor path.
pub fn chain_matches(matches: &[TagMatch], orientation: bool) -> Vec<TagMatch> {
    let mut matches = matches.to_vec();
    matches.sort_by_key(|m| (m.query_position, m.reference_position, m.tag_id));
    if matches.is_empty() {
        return matches;
    }

    let mut score = vec![1u32; matches.len()];
    let mut previous = vec![usize::MAX; matches.len()];
    for i in 0..matches.len() {
        for j in 0..i {
            // Query coordinates are monotonically increasing after sorting.
            let ordered = if orientation {
                matches[j].reference_position < matches[i].reference_position
            } else {
                matches[j].reference_position > matches[i].reference_position
            };
            if !ordered {
                continue;
            }
            let candidate = score[j].saturating_add(1);
            if candidate > score[i]
                || (candidate == score[i] && previous[i] != usize::MAX && j < previous[i])
            {
                score[i] = candidate;
                previous[i] = j;
            }
        }
    }

    let mut best = 0;
    for i in 1..score.len() {
        if score[i] > score[best] {
            best = i;
        }
    }
    let mut chain = Vec::new();
    let mut cursor = best;
    while cursor != usize::MAX {
        chain.push(matches[cursor]);
        cursor = previous[cursor];
    }
    chain.reverse();
    chain
}

#[cfg(test)]
mod tests {
    use super::*;

    fn m(q: u64, r: u64, id: u32) -> TagMatch {
        TagMatch {
            query_position: q,
            reference_position: r,
            tag_id: id,
            reference_contig_id: 1,
        }
    }

    #[test]
    fn chains_forward_and_reverse_matches() {
        let forward = [m(10, 100, 1), m(20, 200, 2), m(30, 150, 3), m(40, 300, 4)];
        assert_eq!(chain_matches(&forward, true).len(), 3);

        let reverse = [m(10, 300, 1), m(20, 200, 2), m(30, 150, 3), m(40, 100, 4)];
        let chain = chain_matches(&reverse, false);
        assert_eq!(chain.len(), 4);
        assert_eq!(chain[0], m(10, 300, 1));
    }

    #[test]
    fn empty_input_has_empty_chain() {
        assert!(chain_matches(&[], true).is_empty());
    }
}
