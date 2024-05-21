use clap::Parser;
use std::fs::File;
// use std::hash::Hash;
use std::collections::HashMap;
use std::io::{self, BufReader, BufWriter};
use std::num::NonZeroUsize;
use noodles::bgzf;
use impg::impg::{Impg, SerializableImpg, AdjustedInterval, check_intervals};
use coitrees::{IntervalTree, BasicCOITree, Interval};
use impg::paf;
use rayon::ThreadPoolBuilder;
use std::io::BufRead;

use itertools::Itertools;

/// Command-line tool for querying overlaps in PAF files.
#[derive(Parser, Debug)]
#[clap(author, version, about)]
struct Args {
    /// Path to the PAF file. If specified without an index, the tool will look for or generate an associated index file.
    #[clap(short='p', long, value_parser)]
    paf_file: Option<String>,

    /// Force the regeneration of the index, even if it already exists.
    #[clap(short='I', long, action)]
    force_reindex: bool,

    /// Target range in the format `seq_name:start-end`.
    #[clap(short='r', long, value_parser)]
    target_range: Option<String>,

    /// Path to the BED file containing target regions.
    #[clap(short='b', long, value_parser)]
    target_bed: Option<String>,

    /// Window size to create PAF files from.
    #[clap(short='w', long, value_parser)]
    window_size: Option<i32>,    

    /// Enable transitive overlap requests.
    #[clap(short='x', long, action)]
    transitive: bool,

    /// Output results in PAF format.
    #[clap(short='P', long, action)]
    output_paf: bool,
        
    /// Print stats about the index.
    #[clap(short='s', long, action)]
    stats: bool,

    /// Number of threads for parallel processing.
    #[clap(short='t', long, value_parser, default_value_t = NonZeroUsize::new(1).unwrap())]
    num_threads: NonZeroUsize,

    /// Check the projected intervals, reporting the wrong ones (slow, useful for debugging).
    #[clap(short='c', long, action)]
    check_intervals: bool,
}

fn main() -> io::Result<()> {
    let args = Args::parse();

    // Configure the global thread pool to use the specified number of threads
    ThreadPoolBuilder::new().num_threads(args.num_threads.into()).build_global().unwrap();

    let impg = match args {
        Args { paf_file: Some(paf), force_reindex: false, .. } => load_or_generate_index(&paf, args.num_threads)?,
        Args { paf_file: Some(paf), force_reindex: true, .. } => generate_index(&paf, args.num_threads)?,
        _ => return Err(io::Error::new(io::ErrorKind::InvalidInput, "A PAF file must be provided")),
    };

    if args.stats {
        print_stats(&impg);
    }

    if let Some(target_range) = args.target_range {
        let (target_name, target_range) = parse_target_range(&target_range)?;
        let results = perform_query(&impg, &target_name, target_range, args.transitive);
        if args.check_intervals {
            let invalid_cigars = check_intervals(&impg, &results);
            if !invalid_cigars.is_empty() {
                for (row, error_reason) in invalid_cigars {
                    eprintln!("{}; {}", error_reason, row);
                }
                panic!("Invalid intervals encountered.");
            }
        }
        if args.output_paf {
            output_results_paf(&impg, results, &target_name, None);
        } else {
            output_results_bed(&impg, results);
        }
    } else if let Some(target_bed) = args.target_bed {
        let targets = parse_bed_file(&target_bed)?;
        for (target_name, target_range, name) in targets {
            let results = perform_query(&impg, &target_name, target_range, args.transitive);
            if args.check_intervals {
                let invalid_cigars = check_intervals(&impg, &results);
                if !invalid_cigars.is_empty() {
                    for (row, error_reason) in invalid_cigars {
                        eprintln!("{}; {}", error_reason, row);
                    }
                    panic!("Invalid intervals encountered.");
                }
            }
            if args.output_paf {
                output_results_paf(&impg, results, &target_name, name);
            } else {
                output_results_bedpe(&impg, results, &target_name, name);
            }
        }
    } else if let Some(window_size) = args.window_size {
        // println!("{window_size}");
        let mut seen: HashMap<u32, Vec<Interval<()>>> = HashMap::new();
        for key in impg.trees.keys().sorted() {
            println!("key: {}", key);
            impg.trees.get(key);
            let target_length = impg.seq_index.get_len_from_id(*key).expect("Target length not found in index");
            println!("target length: {}", target_length);
            let target_name = impg.seq_index.get_name(*key).unwrap();
            println!("target name: {}", target_name);
            let mut i: i32 = 0;
            // check if key in coitree (once per coitree is sufficient)
            let interval_arr = if let Some(interval_arr) = seen.get_mut(key) {
                interval_arr
            } else { // else we insert a new empty vec
                seen.entry(*key).or_insert(Vec::new())
            };
            while i < target_length.try_into().unwrap() {
                println!("IIIII: {}", i);
                let end;
                if i + window_size < target_length.try_into().unwrap() {
                    end = i + window_size;
                } else {
                    end = target_length.try_into().unwrap();
                }
                // if already in coitree, extract overlapping intervals
                // TODO NOT SURE HERE
                // if not in coitree, add a new key with vec
                // TODO NOT SURE HERE

                // transitive stuff
                let results = impg.query(key.clone(), i, end);
                // add new intervals to coitree
                interval_arr.push(Interval::new(i, end, ()));
                output_results_paf(&impg, results, &target_name, None);
                i = i + window_size;
            }
        }
    }   
    Ok(())
}

fn parse_bed_file(bed_file: &str) -> io::Result<Vec<(String, (i32, i32), Option<String>)>> {
    let file = File::open(bed_file)?;
    let reader = BufReader::new(file);
    let mut ranges = Vec::new();

    for line in reader.lines() {
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() < 3 {
            return Err(io::Error::new(io::ErrorKind::InvalidInput, "Invalid BED file format"));
        }

        let (start, end) = parse_range(&parts[1..=2])?;
        let name = parts.get(3).map(|s| s.to_string());
        ranges.push((parts[0].to_string(), (start, end), name));
    }

    Ok(ranges)
}

fn parse_target_range(target_range: &str) -> io::Result<(String, (i32, i32))> {
    let parts: Vec<&str> = target_range.rsplitn(2, ':').collect();
    if parts.len() != 2 {
        return Err(io::Error::new(io::ErrorKind::InvalidInput, "Target range format should be `seq_name:start-end`"));
    }

    let (start, end) = parse_range(&parts[0].split('-').collect::<Vec<_>>())?;
    Ok((parts[1].to_string(), (start, end)))
}

fn parse_range(range_parts: &[&str]) -> io::Result<(i32, i32)> {
    if range_parts.len() != 2 {
        return Err(io::Error::new(io::ErrorKind::InvalidInput, "Range format should be `start-end`"));
    }

    let start = range_parts[0].parse::<i32>().map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "Invalid start value"))?;
    let end = range_parts[1].parse::<i32>().map_err(|_| io::Error::new(io::ErrorKind::InvalidInput, "Invalid end value"))?;

    if start >= end {
        return Err(io::Error::new(io::ErrorKind::InvalidInput, "Start value must be less than end value"));
    }

    Ok((start, end))
}

fn load_or_generate_index(paf_file: &str, num_threads: NonZeroUsize) -> io::Result<Impg> {
    let index_file = format!("{}.impg", paf_file);
    if std::path::Path::new(&index_file).exists() {
        load_index(paf_file)
    } else {
        generate_index(paf_file, num_threads)
    }
}

fn generate_index(paf_file: &str, num_threads: NonZeroUsize) -> io::Result<Impg> {
    let file = File::open(paf_file)?;
    let reader: Box<dyn io::Read> = if [".gz", ".bgz"].iter().any(|e| paf_file.ends_with(e)) {
        Box::new(bgzf::MultithreadedReader::with_worker_count(num_threads, file))
    } else {
        Box::new(file)
    };
    let reader = BufReader::new(reader);
    let records = paf::parse_paf(reader).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("Failed to parse PAF records: {:?}", e)))?;
    let impg = Impg::from_paf_records(&records, paf_file).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("Failed to create index: {:?}", e)))?;

    let index_file = format!("{}.impg", paf_file);
    let serializable = impg.to_serializable();
    let file = File::create(index_file)?;
    let writer = BufWriter::new(file);
    bincode::serialize_into(writer, &serializable).map_err(|e| io::Error::new(io::ErrorKind::Other, format!("Failed to serialize index: {:?}", e)))?;

    Ok(impg)
}

fn load_index(paf_file: &str) -> io::Result<Impg> {
    let index_file = format!("{}.impg", paf_file);
    
    let paf_file_metadata = std::fs::metadata(paf_file)?;
    let index_file_metadata = std::fs::metadata(index_file.clone())?;
    if let (Ok(paf_file_ts), Ok(index_file_ts)) = (paf_file_metadata.modified(), index_file_metadata.modified()) {
        if paf_file_ts > index_file_ts
        {
            eprintln!("WARNING:\tPAF file has been modified since impg index creation.");
        }
    } else {
        eprintln!("WARNING:\tUnable to compare timestamps of PAF file and impg index file. PAF file may have been modified since impg index creation.");
    }

    let file = File::open(index_file)?;
    let reader = BufReader::new(file);
    let serializable: SerializableImpg = bincode::deserialize_from(reader).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("Failed to deserialize index: {:?}", e)))?;
    Ok(Impg::from_paf_and_serializable(paf_file, serializable))
}

fn perform_query(impg: &Impg, target_name: &str, target_range: (i32, i32), transitive: bool) -> Vec<AdjustedInterval> {
    let (target_start, target_end) = target_range;
    let target_id = impg.seq_index.get_id(target_name).expect("Target name not found in index");
    let target_length = impg.seq_index.get_len_from_id(target_id).expect("Target length not found in index");
    if target_end > target_length as i32 {
        panic!("Target range end ({}) exceeds the target sequence length ({})", target_end, target_length);
    }
    if transitive {
        impg.query_transitive(target_id, target_start, target_end)
    } else {
        impg.query(target_id, target_start, target_end)
    }
}

fn output_results_bed(impg: &Impg, results: Vec<AdjustedInterval>) {
    for (overlap, _, _) in results {
        let overlap_name = impg.seq_index.get_name(overlap.metadata).unwrap();
        let (first, last, strand) = if overlap.first <= overlap.last {
            (overlap.first, overlap.last, '+')
        } else {
            (overlap.last, overlap.first, '-')
        };
        println!("{}\t{}\t{}\t.\t{}", overlap_name, first, last, strand);
    }
}

fn output_results_bedpe(impg: &Impg, results: Vec<AdjustedInterval>, target_name: &str, name: Option<String>) {
    for (overlap_query, _, overlap_target) in results {
        let overlap_name = impg.seq_index.get_name(overlap_query.metadata).unwrap();
        let (first, last, strand) = if overlap_query.first <= overlap_query.last {
            (overlap_query.first, overlap_query.last, '+')
        } else {
            (overlap_query.last, overlap_query.first, '-')
        };
        println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t0\t{}\t+",
                 overlap_name, first, last,
                 target_name, overlap_target.first, overlap_target.last,
                 name.as_deref().unwrap_or("."), strand);
    }
}

fn output_results_paf(impg: &Impg, results: Vec<AdjustedInterval>, target_name: &str, name: Option<String>) { 
    let target_length = impg.seq_index.get_len_from_id(impg.seq_index.get_id(target_name).unwrap()).unwrap();  
    for (overlap_query, cigar, overlap_target) in results {
        let overlap_name = impg.seq_index.get_name(overlap_query.metadata).unwrap();
        let (first, last, strand) = if overlap_query.first <= overlap_query.last {
            (overlap_query.first, overlap_query.last, '+')
        } else {
            (overlap_query.last, overlap_query.first, '-')
        };

        let query_length = impg.seq_index.get_len_from_id(overlap_query.metadata).unwrap();  

        let has_m_operation = cigar.iter().any(|op| op.op() == 'M');
        let (matches, block_len) = if has_m_operation {
            // We overestimate the number of matches by counting all M operations
            cigar.iter().fold((0, 0), |(matches, block_len), op| {
                let len = op.len();
                match op.op() {
                    'M' => (matches + len, block_len + len),
                    'I' | 'D' => (matches, block_len + len),
                    _ => (matches, block_len),
                }
            })
        } else {
            cigar.iter().fold((0, 0), |(matches, block_len), op| {
                let len = op.len();
                match op.op() {
                    '=' => (matches + len, block_len + len),
                    'X' | 'I' | 'D' => (matches, block_len + len),
                    _ => (matches, block_len),
                }
            })
        };
        let cigar_str : String = cigar.iter().map(|op| format!("{}{}", op.len(), op.op())).collect();

        match name {
            Some(ref name) => println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tcg:Z:{}\tan:Z:{}",
                                    overlap_name, query_length, first, last, strand,
                                    target_name, target_length, overlap_target.first, overlap_target.last,
                                    matches, block_len, 255, cigar_str, name),
            None => println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tcg:Z:{}",
                                overlap_name, query_length, first, last, strand,
                                target_name, target_length, overlap_target.first, overlap_target.last,
                                matches, block_len, 255, cigar_str),
        }
    }
}

fn print_stats(impg: &Impg) {
    println!("Number of sequences: {}", impg.seq_index.len());
    println!("Number of overlaps: {}", impg.trees.values().map(|tree| tree.len()).sum::<usize>());
}
