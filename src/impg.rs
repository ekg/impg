use std::collections::{HashMap, HashSet};
use coitrees::{BasicCOITree, Interval, IntervalTree};
use crate::paf::{PafRecord, ParseErr, Strand};
use crate::seqidx::SequenceIndex;
use xz2::write::XzEncoder;
use xz2::read::XzDecoder;
use serde::{Serialize, Deserialize};
use std::io::{Write, Read};
use rayon::prelude::*;

/// Parse a CIGAR string into a vector of CigarOp
// Note that the query_delta is negative for reverse strand alignments
#[derive(Clone, Debug)]
#[derive(PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct CigarOp {
    pub val: u32,
}

impl CigarOp {
    pub fn new(len: i32, op: char) -> Result<Self, String> {
        let val = match op {
            '=' => 0,
            'X' => 1,
            'I' => 2,
            'D' => 3,
            _ => return Err(format!("Invalid CIGAR operation: {}", op)),
        };
        Ok(Self { val: (val << 30) | (len as u32) })
    }

    pub fn op(&self) -> char {
        // two most significant bits in the val tell us the op
        match self.val >> 30 {
            0 => '=',
            1 => 'X',
            2 => 'I',
            3 => 'D',
            _ => panic!("Invalid CIGAR operation: {}", self.val >> 30),
        }
    }

    pub fn len(&self) -> i32 {
        (self.val & ((1 << 30) - 1)) as i32
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn target_delta(&self) -> i32 {
        match self.op() {
            '=' | 'X' | 'D' => self.len(),
            'I' => 0,
            _ => panic!("Invalid CIGAR operation: {}", self.op()),
        }
    }

    pub fn query_delta(&self, strand: Strand) -> i32 {
        match self.op() {
            '=' | 'X' | 'I' => if strand == Strand::Forward { self.len() } else { -self.len() },
            'D' => 0,
            _ => panic!("Invalid CIGAR operation: {}", self.op()),
        }
    }
}


#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct QueryMetadata {
    query_id: u32,
    compressed_cigar_ops: Vec<u8>,
    target_start: i32,
    target_end: i32,
    query_start: i32,
    query_end: i32,
    strand: Strand,
}

impl QueryMetadata {
    fn set_cigar_ops(&mut self, cigar_ops: &[CigarOp]) {
        let encoded_cigar_ops = bincode::serialize(cigar_ops).expect("Failed to serialize CIGAR ops");
        let mut encoder = XzEncoder::new(Vec::new(), 9);
        encoder.write_all(&encoded_cigar_ops).expect("Failed to compress CIGAR ops");
        self.compressed_cigar_ops = encoder.finish().expect("Failed to finish compression");
    }

    fn get_cigar_ops(&self) -> Vec<CigarOp> {
        let mut decoder = XzDecoder::new(&self.compressed_cigar_ops[..]);
        let mut decompressed_cigar_ops = Vec::new();
        decoder.read_to_end(&mut decompressed_cigar_ops).expect("Failed to decompress CIGAR ops");
        bincode::deserialize(&decompressed_cigar_ops).expect("Failed to deserialize CIGAR ops")
    }
}

pub type QueryInterval = Interval<u32>;
type TreeMap = HashMap<u32, BasicCOITree<QueryMetadata, u32>>;
pub type SerializableImpg = (HashMap<u32, Vec<SerializableInterval>>, SequenceIndex);

#[derive(Clone, Serialize, Deserialize)]
pub struct SerializableInterval {
    first: i32,
    last: i32,
    metadata: QueryMetadata,
}

#[derive(Clone)]
pub struct Impg {
    pub trees: TreeMap,
    pub seq_index: SequenceIndex,
}

impl Impg {
    pub fn from_paf_records(records: &[PafRecord]) -> Result<Self, ParseErr> {

        let mut seq_index = SequenceIndex::new();
        for record in records {
            seq_index.get_or_insert_id(&record.query_name);
            seq_index.get_or_insert_id(&record.target_name);
        }
        
        let intervals: HashMap<u32, Vec<Interval<QueryMetadata>>> = records.par_iter()
            .filter_map(|record| {
                let cigar_ops = record.cigar.as_ref().map(|x| parse_cigar_to_delta(x)).transpose().ok()?.unwrap_or_else(Vec::new);
                let query_id = seq_index.get_id(&record.query_name).expect("Query name not found in index");
                let target_id = seq_index.get_id(&record.target_name).expect("Target name not found in index");

                let mut query_metadata = QueryMetadata {
                    query_id,
                    compressed_cigar_ops: Vec::new(),
                    target_start: record.target_start as i32,
                    target_end: record.target_end as i32,
                    query_start: record.query_start as i32,
                    query_end: record.query_end as i32,
                    strand: record.strand,
                };
                query_metadata.set_cigar_ops(&cigar_ops);

                Some((target_id, Interval {
                    first: record.target_start as i32,
                    last: record.target_end as i32,
                    metadata: query_metadata,
                }))
            })  // Use fold and reduce to achieve grouping
            .fold(HashMap::new, |mut acc: HashMap<u32, Vec<Interval<QueryMetadata>>>, (target_id, interval)| {
                acc.entry(target_id).or_default().push(interval);
                acc
            })
            .reduce(HashMap::new, |mut acc, part| {
                for (key, value) in part {
                    acc.entry(key).or_default().extend(value);
                }
                acc
            });

        let trees: TreeMap = intervals.into_iter().map(|(target_id, interval_nodes)| {
            (target_id, BasicCOITree::new(interval_nodes.as_slice()))
        }).collect();

        Ok(Self { trees, seq_index })
    }

    pub fn to_serializable(&self) -> SerializableImpg {
        let serializable_trees = self.trees.iter().map(|(target_id, tree)| {
            let intervals = tree.iter().map(|interval| SerializableInterval {
                first: interval.first,
                last: interval.last,
                metadata: interval.metadata.clone(),
            }).collect();
            (*target_id, intervals)
        }).collect();
        (serializable_trees, self.seq_index.clone())
    }

    pub fn from_serializable(serializable: SerializableImpg) -> Self {
        let (serializable_trees, seq_index) = serializable;
        let trees = serializable_trees.into_iter().map(|(target_id, intervals)| {
            let tree = BasicCOITree::new(intervals.iter().map(|interval| Interval {
                first: interval.first,
                last: interval.last,
                metadata: interval.metadata.clone(),
            }).collect::<Vec<_>>().as_slice());
            (target_id, tree)
        }).collect();
        Self { trees, seq_index }
    }

    pub fn query(&self, target_id: u32, range_start: i32, range_end: i32) -> Vec<QueryInterval> {
        let mut results = Vec::new();
        // add the query interval to the results
        results.push(QueryInterval {
            first: range_start,
            last: range_end,
            metadata: target_id,
        });
        if let Some(tree) = self.trees.get(&target_id) {
            tree.query(range_start, range_end, |interval| {
                let metadata = &interval.metadata;
                let (adjusted_start, adjusted_end) = project_target_range_through_alignment(
                    (range_start, range_end),
                    (metadata.target_start, metadata.target_end, metadata.query_start, metadata.query_end, metadata.strand),
                    &metadata.get_cigar_ops()
                );

                let adjusted_interval = QueryInterval {
                    first: adjusted_start,
                    last: adjusted_end,
                    metadata: metadata.query_id,
                };
                results.push(adjusted_interval);
            });
        }
        results
    }

    pub fn query_transitive(&self, target_id: u32, range_start: i32, range_end: i32) -> Vec<QueryInterval> {
        let mut results = Vec::new();
        // add the query interval to the results
        results.push(QueryInterval {
            first: range_start,
            last: range_end,
            metadata: target_id,
        });
        let mut stack = vec![(target_id, range_start, range_end)];
        let mut visited = HashSet::new();

        while let Some((current_target, current_start, current_end)) = stack.pop() {
            if let Some(tree) = self.trees.get(&current_target) {
                tree.query(current_start, current_end, |interval| {
                    let metadata = &interval.metadata;
                    let (adjusted_start, adjusted_end) = project_target_range_through_alignment(
                        (current_start, current_end),
                        (metadata.target_start, metadata.target_end, metadata.query_start, metadata.query_end, metadata.strand),
                        &metadata.get_cigar_ops()
                    );

                    let adjusted_interval = QueryInterval {
                        first: adjusted_start,
                        last: adjusted_end,
                        metadata: metadata.query_id,
                    };
                    results.push(adjusted_interval);

                    if metadata.query_id != current_target {
                        let todo_range = (metadata.query_id, adjusted_start, adjusted_end);
                        if !visited.insert(todo_range) {
                            stack.push(todo_range);
                        }
                    }
                });
            }
        }

        results
    }
}

fn project_target_range_through_alignment(
    target_range: (i32, i32),
    record: (i32, i32, i32, i32, Strand),
    cigar_ops: &[CigarOp],
) -> (i32, i32) {
    let (target_start, _target_end, query_start, query_end, strand) = record;

    let mut target_pos = target_start;
    let mut query_pos = if strand == Strand::Forward { query_start } else { query_end };

    let mut projected_start: Option<i32> = None;
    let mut projected_end: Option<i32> = None;

    for cigar_op in cigar_ops {
        // If the target position is past the end of the range, we can stop
        if target_pos > target_range.1 {
            break;
        }
        match (cigar_op.target_delta(), cigar_op.query_delta(strand)) {
            (0, query_delta) => { // Insertion in query
                if target_pos >= target_range.0 && target_pos <= target_range.1 {
                    projected_start.get_or_insert(query_pos);
                    projected_end = Some(query_pos +
                                         if target_pos <= target_range.1 { 0 } else { query_delta });
                }
                query_pos += query_delta;
            },
            (target_delta, 0) => { // Deletion in target
                let overlap_start = target_pos.max(target_range.0);
                let overlap_end = (target_pos + target_delta).min(target_range.1);

                if overlap_start < overlap_end { // There's an overlap
                    projected_start.get_or_insert(query_pos);
                    projected_end = Some(query_pos); // Deletion does not advance query position
                }

                target_pos += target_delta;
            },
            (target_delta, query_delta) => { // Match or mismatch
                let overlap_start = target_pos.max(target_range.0);
                let overlap_end = (target_pos + target_delta).min(target_range.1);

                if overlap_start < overlap_end { // There's an overlap
                    let overlap_length = overlap_end - overlap_start;
                    let dir = if strand == Strand::Forward { 1 } else { -1 };
                    let query_overlap_start = query_pos + (overlap_start - target_pos) * dir;
                    let query_overlap_end = query_overlap_start + overlap_length * dir;

                    projected_start.get_or_insert(query_overlap_start);
                    projected_end = Some(query_overlap_end);
                }

                target_pos += target_delta;
                query_pos += query_delta;
            },
        }
    }

    if strand == Strand::Reverse {
        std::mem::swap(&mut projected_start, &mut projected_end);
    }
    (projected_start.unwrap_or(query_start), projected_end.unwrap_or(query_pos)) // Changed _query_end to query_pos
}

fn parse_cigar_to_delta(cigar: &str) -> Result<Vec<CigarOp>, ParseErr> {
    let mut ops = Vec::new();
    let mut num_buf = String::new();

    for c in cigar.chars() {
        if c.is_ascii_digit() {
            num_buf.push(c);
        } else {
            let len = num_buf.parse::<i32>().map_err(|_| ParseErr::InvalidCigarFormat)?;
            num_buf.clear(); // Reset the buffer for the next operation
            // raise any error from the cigar op parsing
            let op = CigarOp::new(len, c).map_err(|_| ParseErr::InvalidCigarFormat)?;
            ops.push(op);
        }
    }

    Ok(ops)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::BufReader;
    use crate::paf::parse_paf;

    #[test]
    fn test_project_target_range_through_alignment_forward() {
        let target_range = (100, 200);
        let record = (100, 200, 0, 100, Strand::Forward);
        let cigar_ops = vec![CigarOp::new(100, '=').unwrap()];
        let (start, end) = project_target_range_through_alignment(target_range, record, &cigar_ops);

        assert_eq!(start, 0);
        assert_eq!(end, 100);
    }

    #[test]
    fn test_project_target_range_through_alignment_reverse() {
        let target_range = (100, 200);
        let record = (100, 200, 0, 100, Strand::Reverse);
        let cigar_ops = vec![CigarOp::new(100, '=').unwrap()];
        let (start, end) = project_target_range_through_alignment(target_range, record, &cigar_ops);
        println!("Final result: ({}, {})", start, end);

        assert_eq!(start, 0);
        assert_eq!(end, 100);
    }

    #[test]
    fn test_project_target_range_through_alignment() {
        let cigar_ops = vec![
            CigarOp::new(10, '='), // 10, 60
            CigarOp::new(5, 'I'),  // 10, 65
            CigarOp::new(5, 'D'),  // 15, 65
            CigarOp::new(50, '='), // 65, 115
            CigarOp::new(50, 'I'), // 65, 165
            CigarOp::new(35, '='), // 100, 200
        ];
        // map unwrap over cigar_ops to get rid of the Result
        let cigar_ops: Vec<_> = cigar_ops.into_iter().map(|op| op.unwrap()).collect();
        let base = (0, 100, 50, 200, Strand::Forward);
        {
            let result = project_target_range_through_alignment((0, 100), base, &cigar_ops);
            assert_eq!(result, (50, 200));
        }
        {
            let result = project_target_range_through_alignment((50, 55), base, &cigar_ops);
            assert_eq!(result, (100, 105));
        }
        {
            let result = project_target_range_through_alignment((50, 64), base, &cigar_ops);
            assert_eq!(result, (100, 114));
        }
        {
            let result = project_target_range_through_alignment((65, 65), base, &cigar_ops);
            assert_eq!(result, (115, 115));
        }
        {
            let result = project_target_range_through_alignment((50, 65), base, &cigar_ops);
            assert_eq!(result, (100, 115));
        }
        {
            let result = project_target_range_through_alignment((50, 66), base, &cigar_ops);
            assert_eq!(result, (100, 166));
        }
        {
            let result = project_target_range_through_alignment((70, 95), base, &cigar_ops);
            assert_eq!(result, (170, 195));
        }
    }

    // 1. Simple Forward Projection
    #[test]
    fn test_forward_projection_simple() {
        let target_range = (100, 200);
        let record = (100, 200, 100, 200, Strand::Forward);
        let cigar_ops = vec![CigarOp::new(100, '=').unwrap()];
        let (start, end) = project_target_range_through_alignment(target_range, record, &cigar_ops);
        assert_eq!((start, end), (100, 200));
    }

    // 2. Simple Reverse Projection
    #[test]
    fn test_reverse_projection_simple() {
        let target_range = (100, 200);
        let record = (100, 200, 100, 200, Strand::Reverse);
        let cigar_ops = vec![CigarOp::new(100, '=').unwrap()];
        let (start, end) = project_target_range_through_alignment(target_range, record, &cigar_ops);
        assert_eq!((start, end), (100, 200)); // Adjust for reverse calculation
    }

    // 3. Forward Projection with Insertions
    #[test]
    fn test_forward_projection_with_insertions() {
        let target_range = (50, 150);
        let record = (50, 150, 50, 160, Strand::Forward);
        let cigar_ops = vec![
            CigarOp::new(50, '='), // Match
            CigarOp::new(10, 'I'), // Insertion
            CigarOp::new(50, '='), // Match
        ];
        let cigar_ops: Vec<_> = cigar_ops.into_iter().map(|op| op.unwrap()).collect();
        let (start, end) = project_target_range_through_alignment(target_range, record, &cigar_ops);
        assert_eq!((start, end), (50, 160));
    }

    // 4. Forward Projection with Deletions
    #[test]
    fn test_forward_projection_with_deletions() {
        let target_range = (50, 150);
        let record = (50, 150, 50, 140, Strand::Forward);
        let cigar_ops = vec![
            CigarOp::new(50, '='), // Match
            CigarOp::new(10, 'D'), // Deletion
            CigarOp::new(40, '='), // Match
        ];
        let cigar_ops: Vec<_> = cigar_ops.into_iter().map(|op| op.unwrap()).collect();
        let (start, end) = project_target_range_through_alignment(target_range, record, &cigar_ops);
        assert_eq!((start, end), (50, 140));
    }

    // 5. Reverse Projection with Mixed Operations
    #[test]
    fn test_reverse_projection_with_mixed_operations() {
        let target_range = (150, 250);
        let record = (100, 200, 200, 300, Strand::Reverse);
        let cigar_ops = vec![
            CigarOp::new(50, '='), // 150, 250
            CigarOp::new(10, 'D'), // 150, 260
            CigarOp::new(10, 'I'), // 150, 250
            CigarOp::new(40, '='), // 150, 250
        ];
        let cigar_ops: Vec<_> = cigar_ops.into_iter().map(|op| op.unwrap()).collect();
        let (start, end) = project_target_range_through_alignment(target_range, record, &cigar_ops);
        assert_eq!((start, end), (200, 250));
    }

    // 6. Edge Case Projection
    #[test]
    fn test_edge_case_projection() {
        let target_range = (0, 10);
        let record = (0, 50, 0, 40, Strand::Forward);
        let cigar_ops = vec![
            CigarOp::new(10, '='), // Match
            CigarOp::new(20, 'D'), // Deletion in target
            CigarOp::new(8, '='), // Match
            CigarOp::new(1, 'X'), // Match
            CigarOp::new(1, '='), // Match
            CigarOp::new(10, 'I'), // Insertion in query
            CigarOp::new(10, '='), // Match
        ];
        let cigar_ops: Vec<_> = cigar_ops.into_iter().map(|op| op.unwrap()).collect();
        let (start, end) = project_target_range_through_alignment(target_range, record, &cigar_ops);
        println!("{} {}", start, end);
        assert_eq!((start, end), (0, 10));
    }

    #[test]
    fn test_parse_cigar_to_delta_basic() {
        let cigar = "10=5I5D";
        let cigar_ops = vec![
            CigarOp::new(10, '='),
            CigarOp::new(5, 'I'),
            CigarOp::new(5, 'D'),
        ];
        let cigar_ops: Vec<_> = cigar_ops.into_iter().map(|op| op.unwrap()).collect();
        let ops = parse_cigar_to_delta(cigar).unwrap();
        assert_eq!(ops, cigar_ops);
    }

    #[test]
    fn test_parse_cigar_to_delta_invalid() {
        let cigar = "10=5Q"; // Q is not a valid CIGAR operation
        assert!(parse_cigar_to_delta(cigar).is_err());
    }

    #[test]
    fn test_parse_paf_valid() {
        let paf_data = b"seq1\t100\t10\t20\t+\tt1\t200\t30\t40\t10\t20\t255\tcg:Z:10M\n";
        let reader = BufReader::new(&paf_data[..]);
        let expected_records = vec![
            PafRecord {
                query_name: "seq1".to_string(),
                query_length: 100,
                query_start: 10,
                query_end: 20,
                target_name: "t1".to_string(),
                target_length: 200,
                target_start: 30,
                target_end: 40,
                cigar: Some("10M".to_string()),
                strand: Strand::Forward,
            },
            // Add more test records as needed
        ];
        let records = parse_paf(reader).unwrap();
        assert_eq!(records, expected_records);
    }

}
