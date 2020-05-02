#[macro_use]
extern crate clap;
extern crate bio;
extern crate ncbitaxonomy;

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
    let matches = clap_app!(taxonomy_filter_refseq =>
        (version: ncbitaxonomy::VERSION)
        (author: "Peter van Heusden <pvh@sanbi.axc.za>")
        (about: "Filter NCBI RefSeq FASTA files by taxonomic lineage")
        (@arg TAXDB_URL: -d --db +takes_value "URL for SQLite taxonomy database")
        (@arg NO_PREDICTED: --no_predicted "Don't accept computationally predicted RNAs and proteins (XM_, XR_ and XP_ accessions)")
        (@arg NO_CURATED: --no_curated "Don't accept curated RNAs and proteins (NM_, NR_ and NP_ accessions)")
        (@arg INPUT_FASTA: +required "FASTA file with RefSeq sequences")
        (@arg ANCESTOR_NAME: +required "Name of ancestor to use as ancestor filter")
        (@arg OUTPUT_FASTA: "Output FASTA filename (or stdout if omitted)")
        ).get_matches();

    let no_predicted = match matches.occurrences_of("NO_PREDICTED") {
        0 => false,
        _ => true
    };

    let no_curated = match matches.occurrences_of("NO_CURATED") {
        0 => false,
        _ => true
    };

    let input_fasta_filename = matches.value_of("INPUT_FASTA").unwrap();
    let input_fasta = File::open(input_fasta_filename).unwrap_or_else(|_| panic!("Failed to open input FASTA file ({})", input_fasta_filename));
    let input_fasta_reader = fasta::Reader::new(input_fasta);

    let taxdb_url = if matches.is_present("TAXDB_URL") { Some(matches.value_of("TAXDB_URL").unwrap()) } else { None };
    let taxonomy = NcbiSqliteTaxonomy::new(taxdb_url);

    // the use of Box here is inspired by:
    // https://stackoverflow.com/questions/26378842/how-do-i-overcome-match-arms-with-incompatible-types-for-structs-implementing-sa
    // in short, it is means to present each match 'arm' as returning the same (Box<io::Write>) type
    let output_file = match matches.value_of("OUTPUT_FASTA") {
        Some(name) => Box::new(File::create(name).unwrap_or_else(|_| panic!("Failed to open output file ({})", name))) as Box<dyn io::Write>,
        None => Box::new(io::stdout()) as Box<dyn io::Write>,
    };

    let mut output_fasta = fasta::Writer::new(output_file);

    let ancestor_name = matches.value_of("ANCESTOR_NAME").unwrap();

    if !taxonomy.contains_name(ancestor_name) {
        eprintln!("Taxonomy does not contain an ancestor named {}", ancestor_name);
        process::exit(1);
    }

    for record in input_fasta_reader.records() {
        let record = record.unwrap();
        let description = match record.desc() {
            Some(desc) => desc,
            None => "unknown"
        };
        let division = record.id().as_bytes()[0];
        let species_start = description.find('[').unwrap_or_else(|| panic!("[ missing in description ({})", description));
        let species_end = description.rfind(']').unwrap_or_else(|| panic!("] missing in description ({})", description));
        let species_name = &description[(species_start+1)..species_end];
        if !(no_predicted && (division == b'X' || division == b'Y')) && !(no_curated && (division == b'N' || division == b'A' || division == b'W')) && taxonomy.contains_name(species_name) && taxonomy.is_descendant(species_name, ancestor_name) {
            output_fasta.write(record.id(), record.desc(), wrap(record.seq(), 80).as_slice()).unwrap();
        }
    }
}