extern crate bio;
extern crate ncbitaxonomy;

use clap::Command;

use std::cmp;
use std::fs::File;
use std::io;
use std::process;
use std::vec::Vec;

use bio::io::fasta;
use bio::utils::TextSlice;

use ncbitaxonomy::{NcbiTaxonomy, NcbiSqliteTaxonomy};

// wrap a TextSlice (a rust-bio name for a &[u8] i.e. byte array)
// at a certain width (e.g. 80 to look like NCBI RefSeq)
fn wrap(seq: TextSlice, width: usize) -> Vec<u8> {
    let mut wrapped_seq_vec: Vec<u8> = Vec::new();
    let mut eol= [0; 1];
    let seqlen = seq.len();
    '\n'.encode_utf8(&mut eol);
    for start in (0..seq.len()).step_by(width) {
        let end = cmp::min(start + width, seqlen);
        wrapped_seq_vec.extend_from_slice(&seq[start..end]);
        if end != seqlen {
            // insert a '\n' but only if we are not at the last
            // block of sequence data
            wrapped_seq_vec.extend_from_slice(&eol);
        }
    }
    wrapped_seq_vec
}

pub fn main() {
    // TODO: use functions, write testing suite
    let matches = Command::new("taxonomy_filter_refseq")
        .version(ncbitaxonomy::VERSION)
        .author("Peter van Heusden <pvh@sanbi.ac.za>")
        .about("Filter NCBI RefSeq FASTA files by taxonomic lineage")
        .arg(clap::Arg::new("TAXDB_URL")
            .short('d')
            .long("db")
            .value_name("URL")
            .help("URL for SQLite taxonomy database"))
        .arg(clap::Arg::new("NO_PREDICTED")
            .long("no_predicted")
            .action(clap::ArgAction::SetTrue)
            .help("Don't accept computationally predicted RNAs and proteins (XM_, XR_ and XP_ accessions)"))
        .arg(clap::Arg::new("NO_CURATED")
            .long("no_curated")
            .action(clap::ArgAction::SetTrue)
            .help("Don't accept curated RNAs and proteins (NM_, NR_ and NP_ accessions)"))
        .arg(clap::Arg::new("INPUT_FASTA")
            .required(true)
            .help("FASTA file with RefSeq sequences"))
        .arg(clap::Arg::new("ANCESTOR_NAME")
            .required(true)
            .help("Name of ancestor to use as ancestor filter"))
        .arg(clap::Arg::new("OUTPUT_FASTA")
            .help("Output FASTA filename (or stdout if omitted)"))
        .get_matches();

    let no_predicted = matches.get_flag("NO_PREDICTED");

    let no_curated = matches.get_flag("NO_CURATED");

    let input_fasta_filename = matches.get_one::<String>("INPUT_FASTA").expect("INPUT_FASTA is required").as_str();
    let input_fasta = File::open(input_fasta_filename).unwrap_or_else(|_| panic!("Failed to open input FASTA file ({})", input_fasta_filename));
    let input_fasta_reader = fasta::Reader::new(input_fasta);

    let taxdb_url = matches.get_one::<String>("TAXDB_URL").map(|s| s.as_str());
    let taxonomy = NcbiSqliteTaxonomy::new(taxdb_url);

    // the use of Box here is inspired by:
    // https://stackoverflow.com/questions/26378842/how-do-i-overcome-match-arms-with-incompatible-types-for-structs-implementing-sa
    // in short, it is means to present each match 'arm' as returning the same (Box<io::Write>) type
    let output_file = match matches.get_one::<String>("OUTPUT_FASTA") {
        Some(name) => Box::new(File::create(name).unwrap_or_else(|_| panic!("Failed to open output file ({})", name))) as Box<dyn io::Write>,
        None => Box::new(io::stdout()) as Box<dyn io::Write>,
    };

    let mut output_fasta = fasta::Writer::new(output_file);

    let ancestor_name = matches.get_one::<String>("ANCESTOR_NAME").expect("ANCESTOR_NAME is required").as_str();

    if !taxonomy.contains_name(ancestor_name) {
        eprintln!("Taxonomy does not contain an ancestor named {}", ancestor_name);
        process::exit(1);
    }

    for record in input_fasta_reader.records() {
        let record = record.unwrap();
        let description = record.desc().unwrap_or("unknown");
        let division = record.id().as_bytes()[0];
        let species_start = description.find('[').unwrap_or_else(|| panic!("[ missing in description ({})", description));
        let species_end = description.rfind(']').unwrap_or_else(|| panic!("] missing in description ({})", description));
        let species_name = &description[(species_start+1)..species_end];
        if !(no_predicted && (division == b'X' || division == b'Y')) && !(no_curated && (division == b'N' || division == b'A' || division == b'W')) && taxonomy.contains_name(species_name) && taxonomy.is_descendant(species_name, ancestor_name) {
            output_fasta.write(record.id(), record.desc(), wrap(record.seq(), 80).as_slice()).unwrap();
        }
    }
}