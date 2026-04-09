extern crate bio;
extern crate flate2;
extern crate seq_io;
extern crate ncbitaxonomy;

use std::fs::{File, create_dir};
use std::io;
use std::io::{Read, Write, BufRead, BufWriter, BufReader};
use std::path::Path;
use std::process;
use std::vec::Vec;
use std::fmt;
use std::collections::HashMap;

use flate2::Compression;
use flate2::read::GzDecoder;
use flate2::write::GzEncoder;
use clap::{Arg, ArgGroup, Command};
use seq_io::fastq::Record;
use ncbitaxonomy::{NcbiTaxonomy, NcbiSqliteTaxonomy};

enum FilterTool {
    Centrifuge,
    Kraken2
}

impl fmt::Display for FilterTool {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match *self {
            FilterTool::Centrifuge => write!(f, "FilterTool::Centrifuge"),
            FilterTool::Kraken2 => write!(f, "FilterTool::Kraken2")
        }
    }
}

fn filter_fastq(fastq_filename: &Path, tax_report_filename: &str,
                taxonomy: &dyn NcbiTaxonomy,
                output_dir: &Path, filter_tool: &FilterTool, ancestor_id: i32) {
    if !output_dir.exists() {
        create_dir(output_dir).unwrap_or_else(|_| panic!("Failed to create output dir {}", output_dir.display()));
    }
    let fastq_file = File::open(fastq_filename).unwrap_or_else(|_| panic!("Failed to open input FASTQ file ({})", fastq_filename.display()));
    let fastq_decoder: Box<dyn Read> = if fastq_filename.to_str().unwrap().ends_with(".gz") {
        Box::new(GzDecoder::new(fastq_file))
    } else {
        Box::new(fastq_file)
    };
    let mut fastq_reader = seq_io::fastq::Reader::new(BufReader::new(fastq_decoder));

    let tax_report_file = File::open(tax_report_filename).unwrap_or_else(|_| panic!("Failed to open input Centrifuge file ({})", tax_report_filename));
    let tax_report_reader = io::BufReader::new(tax_report_file);

    let mut read_valid: HashMap<String, i32> = HashMap::new();
    for line in tax_report_reader.lines() {
        let line = line.expect("Unable to read line from centrifuge file");
        let fields = line.split('\t').collect::<Vec<&str>>();

        match filter_tool {
            FilterTool::Centrifuge => {
                if line.starts_with("readID") {
                    // skip the header
                    continue;
                }
                let id = fields[0].to_owned();
                let score = fields[3].parse::<i32>().unwrap_or_else(|_| panic!("Failed to read score from ({})", fields[3]));
                let current_score = match read_valid.get(&id) {
                    Some(value) => *value,
                    None => 0
                };
                if score >= current_score {
                    let taxid = fields[2].parse::<i32>().unwrap();

                    if taxonomy.is_descendant_taxid(taxid, ancestor_id) {
                        read_valid.insert(id, score);
                    } else if score > current_score {
                        // only reset this to zero if this non-descendant taxid is a better fit
                        read_valid.insert(id, 0);
                    }
                }
            },
            FilterTool::Kraken2 => {
                // kraken2 adds an entry each time it sees a read, so for paired end
                // reads there are 2 entries which might not agree with each other
                //
                // this code treats read as invalid if it *ever* shows up as invalid,
                // i.e. if either read in a pair is invalid
                let is_classified = fields[0];
                let id = fields[1].to_owned();
                if is_classified == "U" {
                    // this is an unclassified read
                    read_valid.insert(id, 0);
                } else {
                    // is_classified is C
                    if !read_valid.contains_key(&id) || *(read_valid.get(&id).unwrap()) != 0 {
                        // only insert key if it is either new or was not previously noted
                        // as unclassified or not a descendant

                        let name_or_taxid = fields[2];
                        let taxid = if name_or_taxid.contains("(taxid") {
                            // taxon name output format
                            let name_parts = name_or_taxid.split(' ').collect::<Vec<&str>>();
                            let last_part = name_parts.last().unwrap();
                            let taxid_part: &str = last_part.split(')').collect::<Vec<&str>>()[0];
                            taxid_part.parse::<i32>().unwrap()
                        } else {
                            name_or_taxid.parse::<i32>().unwrap()
                        };

                        if taxonomy.is_descendant_taxid(taxid, ancestor_id) {
                            read_valid.insert(id, 1000);  // make up a score for kraken2
                        } else  {
                            read_valid.insert(id, 0);
                        }

                    }
                }
            },
        };
    }

    let mut valid_records = 0;
    let mut total_records = 0;
    let filename_parts: Vec<&str>= fastq_filename.file_name().and_then(|s| s.to_str()).unwrap().split('.').collect();
    let output_filename = output_dir.to_str().unwrap().to_owned() + "/" + filename_parts[0] + ".filtered." + &filename_parts[1..].join(".");
    let output_file = File::create(&output_filename).unwrap_or_else(|_| panic!("Failed to create output file: {}", output_filename));
    let output_encoder: Box<dyn Write> = if output_filename.ends_with(".gz") {
        Box::new(GzEncoder::new(output_file, Compression::default()))
    } else {
        Box::new(output_file)
    };
    let mut output_writer = BufWriter::new(output_encoder );
    while let Some(result) = fastq_reader.next() {
        let record = result.expect("Error reading record");
        let id = record.id().unwrap();
        total_records += 1;
        if read_valid.contains_key(id) && *read_valid.get(id).unwrap() > 0 {
            record.write_unchanged(&mut output_writer).unwrap_or_else(|_| panic!("Failed to write record to output file"));
            valid_records += 1;
        }
    }
    eprintln!("{} records written out of {} total records", valid_records, total_records);
}

pub fn main() {
    // TODO: write test suite
    let matches = Command::new("taxonomy_filter_fastq")
        .version(ncbitaxonomy::VERSION)
        .author("Peter van Heusden <pvh@sanbi.ac.za>")
        .about("Filter FASTQ files whose reads have been classified by Centrifuge or Kraken2, only retaining reads in taxa descending from given ancestor")
        .arg(Arg::new("TAXDB_URL")
            .short('b')
            .long("db")
            .value_name("URL")
            .help("URL for SQLite taxonomy database"))
        .arg(Arg::new("ANCESTOR_ID")
            .short('A')
            .long("ancestor_taxid")
            .value_name("TAXID")
            .required(true)
            .help("Taxonomy ID of ancestor to use as ancestor filter"))
        .group(ArgGroup::new("filter_tool")
            .arg("centrifuge")
            .arg("kraken2")
            .required(true))
        .arg(Arg::new("centrifuge")
            .short('C')
            .long("centrifuge")
            .action(clap::ArgAction::SetTrue)
            .help("Filter using report from Centrifuge"))
        .arg(Arg::new("kraken2")
            .short('K')
            .long("kraken2")
            .action(clap::ArgAction::SetTrue)
            .help("Filter using report from Kraken2"))
        .arg(Arg::new("OUTPUT_DIR")
            .short('o')
            .long("output_dir")
            .value_name("OUTDIR")
            .help("Directory to deposit filtered output files in"))
        .arg(Arg::new("TAXONOMY_REPORT_FILENAME")
            .short('F')
            .long("tax_report_filename")
            .value_name("FILE")
            .required(true)
            .help("Output from Kraken2 (default) or Centrifuge"))
        .arg(Arg::new("INPUT_FASTQ")
            .value_name("FILE")
            .num_args(1..)
            .required(true)
            .help("FASTQ file(s) with reads to filter"))
        .get_matches();

    let output_dir_str = matches.get_one::<String>("OUTPUT_DIR");
    let output_dir = match output_dir_str {
        Some(path) => Path::new(path.as_str()),
        None => Path::new(".")
    };

    let filter_tool = if matches.get_flag("centrifuge") {
        FilterTool::Centrifuge
    } else {
        FilterTool::Kraken2  // default to kraken2
    };
    eprintln!("filter tool {}", filter_tool);

    let ancestor_id_str = matches.get_one::<String>("ANCESTOR_ID").expect("ANCESTOR_ID is required");
    let ancestor_id = ancestor_id_str.parse::<i32>().unwrap_or_else(|_| panic!("Failed to interpret ({}) as a taxonomy ID", ancestor_id_str));

    let tax_report_filename = matches.get_one::<String>("TAXONOMY_REPORT_FILENAME").expect("TAXONOMY_REPORT_FILENAME is required");

    let taxdb_url = matches.get_one::<String>("TAXDB_URL").map(|s| s.as_str());

    let taxonomy = NcbiSqliteTaxonomy::new(taxdb_url);
    if !taxonomy.contains_id(ancestor_id) {
        eprintln!("Taxonomy does not contain an ancestor with taxid {}", ancestor_id);
        process::exit(1);
    }

    let input_files: Vec<&str> = matches.get_many::<String>("INPUT_FASTQ")
        .expect("INPUT_FASTQ is required")
        .map(|s| s.as_str())
        .collect();
    for input_file in input_files.iter() {
        let input_file_path = Path::new(input_file);
        eprintln!("processing {}", input_file_path.file_name().and_then(|s| s.to_str()).unwrap());
        filter_fastq(input_file_path, tax_report_filename, &taxonomy,
                     output_dir, &filter_tool, ancestor_id);
    }
}
