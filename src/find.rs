use std::collections::HashMap;

use crate::types::Sequence;
use bio::alignment::distance::simd::{bounded_levenshtein, hamming};

const MAX_LDIST_LEN: usize = 5000;

pub fn find_motif_matches_no_rot(motif: &Sequence, seq: &Sequence) -> Vec<(usize, usize)> {
    let mut results = Vec::new();
    let seq_len = seq.len();
    let motif_len = motif.len();

    let mut idx = 0;
    while idx <= seq_len - motif_len {
        if seq[idx..idx + motif_len] == **motif {
            results.push((idx, idx + motif_len));
            idx += motif_len; // Move past the current match
        } else {
            idx += 1;
        }
    }

    results
}

pub fn find_motif_matches(motif: &Sequence, seq: &Sequence) -> Vec<(usize, usize)> {
    let mut motif = motif.clone();
    let mut results = Vec::new();

    // Do for all rotations
    for _ in 0..motif.len() {
        let rot_results = find_motif_matches_no_rot(&motif, seq);
        let rot_results = merge_overlapping_intervals(&rot_results);
        results.extend(rot_results);
        motif.rotate_left(1);
    }

    // sort by start position
    results.sort_by_key(|t| t.0);

    results
}

pub fn merge_nearby_intervals(
    intervals: &[(usize, usize)],
    distance: usize,
) -> Vec<(usize, usize)> {
    if intervals.is_empty() {
        return vec![];
    }

    let mut merged = Vec::new();
    let mut start = intervals[0].0;
    let mut end = intervals[0].1;

    for &(next_start, next_end) in intervals.iter().skip(1) {
        if next_start <= end + distance {
            end = end.max(next_end);
        } else {
            merged.push((start, end));
            start = next_start;
            end = next_end;
        }
    }

    merged.push((start, end));
    merged
}

pub fn merge_overlapping_intervals(intervals: &[(usize, usize)]) -> Vec<(usize, usize)> {
    merge_nearby_intervals(intervals, 0)
}

pub fn purity_score(seq: &[u8], motif: &[u8], threshold: f64) -> f64 {
    if seq.len() < MAX_LDIST_LEN {
        purity_score_ldist(seq, motif, threshold)
    } else {
        purity_score_hamming(seq, motif)
    }
}

pub fn purity_score_hamming(seq: &[u8], motif: &[u8]) -> f64 {
    let pure_repeat = motif
        .iter()
        .cycle()
        .take(seq.len() + motif.len() - 1)
        .cloned()
        .collect::<Vec<u8>>();
    let mut best_score: u64 = u64::MAX;
    for i in 0..motif.len() {
        let slice = &pure_repeat[i..i + seq.len()];
        let score = hamming(seq, slice);
        best_score = best_score.min(score);
    }
    1.0 - (best_score as f64 / seq.len() as f64)
}

pub fn purity_score_ldist(seq: &[u8], motif: &[u8], threshold: f64) -> f64 {
    let pure_repeat: Sequence = motif
        .iter()
        .cycle()
        .take(seq.len() + motif.len() - 1)
        .cloned()
        .collect::<Vec<u8>>()
        .into();

    let motif_len = motif.len() as u32;

    let edits_bound = (seq.len() as f64 * (1.0 - threshold)) as u32 + motif_len + 1;
    let ldist = bounded_levenshtein(seq, &pure_repeat, edits_bound).unwrap_or(edits_bound)
        - (motif_len - 1);
    // let ldist = levenshtein(&seq, &pure_repeat) - (motif_len - 1);
    let norm_ldist = ldist as f64 / seq.len() as f64;
    1.0 - norm_ldist
}

pub fn filter_intervals_by_length(
    intervals: &[(usize, usize)],
    min_length: usize,
) -> Vec<(usize, usize)> {
    intervals
        .iter()
        .filter(|t| t.1 - t.0 >= min_length)
        .cloned()
        .collect()
}

pub fn filter_intervals_by_purity(
    intervals: &[(usize, usize)],
    threshold: f64,
    motif: &Sequence,
    seq: &Sequence,
) -> Vec<(usize, usize)> {
    intervals
        .iter()
        .filter(|t| {
            let seq = &seq[t.0..t.1];
            purity_score(seq, motif, threshold) >= threshold
        })
        .cloned()
        .collect()
}

pub fn find_repeat_variants(
    motif: &Sequence,
    seq: &Sequence,
    min_units: usize,
    threshold: f64,
    distance: usize,
) -> Vec<(usize, usize)> {
    // Find all intervals where the motif matches
    let intervals = find_motif_matches(motif, seq);

    // Merge intervals that are nearby to each other
    let intervals: Vec<(usize, usize)> = merge_nearby_intervals(&intervals, distance);

    // Filter out intervals with low purity
    let intervals = filter_intervals_by_purity(&intervals, threshold, motif, seq);

    // Use DP algorithm to find the longest set of tracts while maintaining purity
    let intervals = merge_intervals_by_purity(intervals, threshold, motif, seq);

    // Ensure purity is maintained after merging
    assert_eq!(
        intervals,
        filter_intervals_by_purity(&intervals, threshold, motif, seq)
    );

    // Merge intervals that are close to each other, again
    let intervals: Vec<(usize, usize)> = merge_nearby_intervals(&intervals, distance);

    // Remove intervals that are too short
    let min_length = min_units * motif.len();
    let intervals = filter_intervals_by_length(&intervals, min_length);

    // Filter out intervals that have low purity
    filter_intervals_by_purity(&intervals, threshold, motif, seq)
}

fn merge_intervals_by_purity(
    intervals: Vec<(usize, usize)>,
    threshold: f64,
    motif: &Sequence,
    seq: &Sequence,
) -> Vec<(usize, usize)> {
    let n = intervals.len();

    if n == 0 {
        return vec![];
    }

    let mut dp = vec![vec![0; n]; n];
    let mut merge = vec![vec![0; n]; n];

    // Initialize base cases
    for i in 0..n {
        dp[i][i] = intervals[i].1 - intervals[i].0;
    }

    // Fill the dp array
    for i in 0..n {
        for j in (i + 1)..n {
            // Consider merging interval i and j
            let merged_start = intervals[i].0.min(intervals[j].0);
            let merged_end = intervals[i].1.max(intervals[j].1);
            let merged_seq = &seq[merged_start..merged_end];
            let mut merged_length = 0;
            if purity_score(merged_seq, motif, threshold) >= threshold {
                merged_length = merged_end - merged_start;
            }

            // Consider other solutions between intervals i and j
            let mut best_mid = i;
            let mut best_mid_length = dp[i][i] + dp[i + 1][j];
            for mid in (i + 1)..j {
                let mid_length = dp[i][mid] + dp[mid + 1][j];
                if mid_length > best_mid_length {
                    best_mid = mid;
                    best_mid_length = mid_length;
                }
            }

            if best_mid_length > merged_length {
                dp[i][j] = best_mid_length;
                merge[i][j] = best_mid as i32;
            } else {
                dp[i][j] = merged_length;
                merge[i][j] = -1;
            }
        }
    }

    retrieve_intervals(&intervals, &merge, 0, n - 1)
}

fn retrieve_intervals(
    intervals: &Vec<(usize, usize)>,
    merge: &Vec<Vec<i32>>,
    i: usize,
    j: usize,
) -> Vec<(usize, usize)> {
    let mut cache = HashMap::new();
    retrieve_intervals_cached(intervals, merge, i, j, &mut cache)
}

fn retrieve_intervals_cached(
    intervals: &Vec<(usize, usize)>,
    merge: &Vec<Vec<i32>>,
    i: usize,
    j: usize,
    cache: &mut HashMap<(usize, usize), Vec<(usize, usize)>>,
) -> Vec<(usize, usize)> {
    if let Some(result) = cache.get(&(i, j)) {
        return result.clone();
    }

    let result = if i == j {
        vec![intervals[i]]
    } else if merge[i][j] == -1 {
        vec![(
            intervals[i].0.min(intervals[j].0),
            intervals[i].1.max(intervals[j].1),
        )]
    } else {
        let mid = merge[i][j] as usize;
        let mut left = retrieve_intervals_cached(intervals, merge, i, mid, cache);
        let mut right = retrieve_intervals_cached(intervals, merge, mid + 1, j, cache);
        left.append(&mut right);
        left
    };

    cache.insert((i, j), result.clone());
    result
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_purity_score_hamming() {
        // Test case where seq is a perfect repeat of motif
        let seq = b"ATATATAT";
        let motif = b"AT";
        let score = super::purity_score_hamming(seq, motif);
        assert_eq!(score, 1.0);

        // Test case where seq is a perfect repeat of motif but rotated
        let seq = b"TATATATA";
        let motif = b"AT";
        let score = super::purity_score_hamming(seq, motif);
        assert_eq!(score, 1.0);

        // Test case where seq is a one off
        let seq = b"ATAXATAT";
        let motif = b"AT";
        let score = super::purity_score_hamming(seq, motif);
        assert_eq!(score, 1.0 - 1.0 / 8.0);

        // Test case where seq is a completely off
        let seq = b"CGCGCGCG";
        let motif = b"AT";
        let score = super::purity_score_hamming(seq, motif);
        assert_eq!(score, 0.0);
    }

    #[test]
    fn test_purity_score_ldist() {
        // Test case where seq is a perfect repeat of motif
        let seq = b"ATATATAT";
        let motif = b"AT";
        let score = super::purity_score_ldist(seq, motif, 0.0);
        assert_eq!(score, 1.0);

        // Test case where seq is a perfect repeat of motif but rotated
        let seq = b"TATATATA";
        let motif = b"AT";
        let score = super::purity_score_ldist(seq, motif, 0.0);
        assert_eq!(score, 1.0);

        // Test case where seq is a one off
        let seq = b"ATAXATAT";
        let motif = b"AT";
        let score = super::purity_score_ldist(seq, motif, 0.0);
        assert_eq!(score, 1.0 - 1.0 / 8.0);

        // Test case where seq is a one off, but with insertion
        let seq = b"ATATXATA";
        let motif = b"AT";
        let score = super::purity_score_ldist(seq, motif, 0.0);
        assert_eq!(score, 1.0 - 1.0 / 8.0);

        // Test case where seq is a completely off
        let seq = b"CGCGCGCG";
        let motif = b"AT";
        let score = super::purity_score_ldist(seq, motif, 0.0);
        assert_eq!(score, 0.0);
    }

    #[test]
    fn test_find_motif_matches_no_rot() {
        // Test case where seq is a perfect repeat of motif
        let seq = b"ATATATAT";
        let motif = b"AT";
        let matches = super::find_motif_matches_no_rot(&motif.into(), &seq.into());
        assert_eq!(matches, vec![(0, 2), (2, 4), (4, 6), (6, 8)]);

        // Test case where seq is imperfect repeat of motif
        let seq = b"ATATAATAT";
        let motif = b"AT";
        let matches = super::find_motif_matches_no_rot(&motif.into(), &seq.into());
        assert_eq!(matches, vec![(0, 2), (2, 4), (5, 7), (7, 9)]);

        // Test case where seq has no matches
        let seq = b"CGCGCGCG";
        let motif = b"AT";
        let matches = super::find_motif_matches_no_rot(&motif.into(), &seq.into());
        assert_eq!(matches, vec![]);
    }

    #[test]
    fn test_find_motif_matches() {
        // Test case where seq is a perfect repeat of motif
        let seq = b"ATATATAT";
        let motif = b"AT";
        let matches = super::find_motif_matches(&motif.into(), &seq.into());
        assert_eq!(matches, vec![(0, 8), (1, 7)]);

        // Test case where seq is imperfect repeat of motif
        let seq = b"ATATAATAT";
        let motif = b"AT";
        let matches = super::find_motif_matches(&motif.into(), &seq.into());
        assert_eq!(matches, vec![(0, 4), (1, 5), (5, 9), (6, 8)]);

        // Test case where seq has no matches
        let seq = b"CGCGCGCG";
        let motif = b"AT";
        let matches = super::find_motif_matches(&motif.into(), &seq.into());
        assert_eq!(matches, vec![]);
    }

    #[test]
    fn test_merge_nearby_intervals() {
        // Test case where there are no intervals
        let intervals = vec![];
        let merged = super::merge_nearby_intervals(&intervals, 0);
        assert_eq!(merged, vec![]);

        // Test case where there is only one interval
        let intervals = vec![(0, 2)];
        let merged = super::merge_nearby_intervals(&intervals, 0);
        assert_eq!(merged, vec![(0, 2)]);

        // Test case where there are two intervals that are not touching
        let intervals = vec![(0, 2), (4, 6)];
        let merged = super::merge_nearby_intervals(&intervals, 0);
        assert_eq!(merged, vec![(0, 2), (4, 6)]);

        // Test case where there are two intervals that are touching
        let intervals = vec![(0, 2), (2, 4)];
        let merged = super::merge_nearby_intervals(&intervals, 0);
        assert_eq!(merged, vec![(0, 4)]);

        // Test case where there are three intervals that are touching
        let intervals = vec![(0, 2), (2, 4), (4, 6)];
        let merged = super::merge_nearby_intervals(&intervals, 0);
        assert_eq!(merged, vec![(0, 6)]);

        // Test case where there are three intervals that are close, but not touching
        let intervals = vec![(0, 2), (3, 5), (6, 8)];
        let merged = super::merge_nearby_intervals(&intervals, 0);
        assert_eq!(merged, vec![(0, 2), (3, 5), (6, 8)]);

        // Test case where there are two overlapping intervals
        let intervals = vec![(0, 4), (2, 6)];
        let merged = super::merge_nearby_intervals(&intervals, 0);
        assert_eq!(merged, vec![(0, 6)]);

        // Test case where there are three overlapping intervals
        let intervals = vec![(0, 4), (2, 6), (5, 8)];
        let merged = super::merge_nearby_intervals(&intervals, 0);
        assert_eq!(merged, vec![(0, 8)]);

        // Test case where there are two overlapping intervals and a non-overlapping interval
        let intervals = vec![(0, 4), (2, 6), (7, 9)];
        let merged = super::merge_nearby_intervals(&intervals, 0);
        assert_eq!(merged, vec![(0, 6), (7, 9)]);

        // Test case where there are two overlapping intervals and a non-overlapping touching interval
        let intervals = vec![(0, 4), (2, 6), (6, 9)];
        let merged = super::merge_nearby_intervals(&intervals, 0);
        assert_eq!(merged, vec![(0, 9)]);

        // Test case where there are two overlapping intervals but the second has a smaller end
        let intervals = vec![(0, 10), (2, 5)];
        let merged = super::merge_nearby_intervals(&intervals, 0);
        assert_eq!(merged, vec![(0, 10)]);

        // Test case where there are three overlapping intervals but the second two have a smaller end
        let intervals = vec![(0, 10), (2, 5), (6, 8)];
        let merged = super::merge_nearby_intervals(&intervals, 0);
        assert_eq!(merged, vec![(0, 10)]);
    }

    #[test]
    fn test_filter_intervals_by_length() {
        // Test case where there are no intervals
        let intervals = vec![];
        let filtered = super::filter_intervals_by_length(&intervals, 0);
        assert_eq!(filtered, vec![]);

        // Test case where there is only one interval
        let intervals = vec![(0, 2)];
        let filtered = super::filter_intervals_by_length(&intervals, 0);
        assert_eq!(filtered, vec![(0, 2)]);

        // Test case where there is only one interval that is too short
        let intervals = vec![(0, 2)];
        let filtered = super::filter_intervals_by_length(&intervals, 3);
        assert_eq!(filtered, vec![]);

        // Test case where there are two intervals, one of which is too short
        let intervals = vec![(0, 2), (3, 6)];
        let filtered = super::filter_intervals_by_length(&intervals, 3);
        assert_eq!(filtered, vec![(3, 6)]);

        // Test case where there are two intervals, both of which are too short
        let intervals = vec![(0, 2), (3, 5)];
        let filtered = super::filter_intervals_by_length(&intervals, 3);
        assert_eq!(filtered, vec![]);

        // Test case where there are two intervals, both of which are long enough
        let intervals = vec![(0, 5), (3, 10)];
        let filtered = super::filter_intervals_by_length(&intervals, 5);
        assert_eq!(filtered, vec![(0, 5), (3, 10)]);
    }

    #[test]
    fn test_filter_intervals_by_purity() {
        // Test case where there are no intervals
        let intervals = vec![];
        let seq = b"ATATATAT";
        let motif = b"AT";
        let filtered =
            super::filter_intervals_by_purity(&intervals, 0.9, &motif.into(), &seq.into());
        assert_eq!(filtered, vec![]);

        // Test case where there is only one interval
        let intervals = vec![(0, 2)];
        let seq = b"ATATATAT";
        let motif = b"AT";
        let filtered =
            super::filter_intervals_by_purity(&intervals, 0.9, &motif.into(), &seq.into());
        assert_eq!(filtered, vec![(0, 2)]);

        // Test case where there is only one interval that is not pure enough
        let intervals = vec![(0, 2)];
        let seq = b"ACATATAT";
        let motif = b"AT";
        let filtered =
            super::filter_intervals_by_purity(&intervals, 0.90, &motif.into(), &seq.into());
        assert_eq!(filtered, vec![]);

        // Test case where there are two intervals, one of which is not pure enough
        let intervals = vec![(0, 2), (3, 6)];
        let seq = b"ATACGTAT";
        let motif = b"AT";
        let filtered =
            super::filter_intervals_by_purity(&intervals, 0.95, &motif.into(), &seq.into());
        assert_eq!(filtered, vec![(0, 2)]);

        // Test case where there are two intervals, both of which are pure enough
        let intervals = vec![(0, 5), (1, 8)];
        let seq = b"ATATTTAT";
        let motif = b"AT";
        let filtered =
            super::filter_intervals_by_purity(&intervals, 0.85, &motif.into(), &seq.into());
        assert_eq!(filtered, vec![(0, 5), (1, 8)]);
    }

    #[test]
    fn test_merge_intervals_by_purity() {
        // Test case where there are no intervals
        let intervals = vec![];
        let seq = b"ATATATAT";
        let motif = b"AT";
        let merged = super::merge_intervals_by_purity(intervals, 0.9, &motif.into(), &seq.into());
        assert_eq!(merged, vec![]);

        // Test case where there is only one interval
        let intervals = vec![(0, 2)];
        let seq = b"ATATATAT";
        let motif = b"AT";
        let merged = super::merge_intervals_by_purity(intervals, 0.9, &motif.into(), &seq.into());
        assert_eq!(merged, vec![(0, 2)]);

        // Test case with multiple intervals
        let intervals = vec![(0, 3), (3, 6), (6, 9)];
        let seq: &[u8; 10] = b"ATATATATAC";
        let motif = b"AT";
        let merged = super::merge_intervals_by_purity(intervals, 0.90, &motif.into(), &seq.into());
        assert_eq!(merged, vec![(0, 9)]);
    }
}
