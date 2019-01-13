#[macro_use]
extern crate clap;

extern crate bio;

use std::cmp;
use std::fs::File;
use std::io;
use std::path::Path;
use std::process;
use std::vec::Vec;

use bio::io::fasta;
use bio::utils::TextSlice;

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
    let matches = clap_app!(taxonomy_filter_refseq =>
        (version: "0.1.2")
        (author: "Peter van Heusden <pvh@sanbi.axc.za>")
        (about: "Filter NCBI RefSeq FASTA files by taxonomic lineage")
        (@arg TAXONOMY_FILENAME_PREFIX: -t --tax_prefix +takes_value "String to prepend to names of nodes.dmp and names.dmp")
        (@arg INPUT_FASTA: +required "FASTA file with RefSeq sequences")
        (@arg TAXONOMY_DIR: +required "Directory containing the NCBI taxonomy nodes.dmp and names.dmp files")
        (@arg ANCESTOR_NAME: +required "Name of ancestor to use as ancestor filter")
        (@arg OUTPUT_FASTA: "Output FASTA filename (or stdout if omitted)")
        ).get_matches();

    let input_fasta_filename = matches.value_of("INPUT_FASTA").unwrap();
    let input_fasta = File::open(input_fasta_filename).unwrap_or_else(|_| panic!("Failed to open input FASTA file ({})", input_fasta_filename));
    let input_fasta_reader = fasta::Reader::new(input_fasta);

    let ncbi_taxonomy_path = Path::new(matches.value_of("TAXONOMY_DIR").unwrap());

    let tax_prefix = match matches.value_of("TAXONOMY_FILENAME_PREFIX") {
        Some(name) => name,
        None => ""
    }.to_string();

    let nodes_path = ncbi_taxonomy_path.join(tax_prefix.clone() + "nodes.dmp");
    if ! nodes_path.exists() {
        eprintln!("NCBI Taxonomy {}nodes.dmp file not found in {}", tax_prefix, ncbi_taxonomy_path.to_str().unwrap());
        process::exit(1);
    }

    let names_path = ncbi_taxonomy_path.join(tax_prefix.clone() + "names.dmp");
    if ! names_path.exists() {
        eprintln!("NCBI Taxonomy {}names.dmp file not found in {}", tax_prefix, ncbi_taxonomy_path.to_str().unwrap());
        process::exit(1);
    }

    // the use of Box here is inspired by:
    // https://stackoverflow.com/questions/26378842/how-do-i-overcome-match-arms-with-incompatible-types-for-structs-implementing-sa
    // in short, it is means to present each match 'arm' as returning the same (Box<io::Write>) type
    let output_file = match matches.value_of("OUTPUT_FASTA") {
        Some(name) => Box::new(File::create(name).unwrap_or_else(|_| panic!("Failed to open output file ({})", name))) as Box<io::Write>,
        None => Box::new(io::stdout()) as Box<io::Write>,
    };

    let mut output_fasta = fasta::Writer::new(output_file);

    let ancestor_name = matches.value_of("ANCESTOR_NAME").unwrap();

    let taxonomy = ncbitaxonomy::NcbiTaxonomy::from_ncbi_files(
        nodes_path.as_path().to_str().unwrap(),
        names_path.as_path().to_str().unwrap()).expect("Failed to load NCBI Taxonomy");

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
        let species_start = description.find('[').unwrap_or_else(|| panic!("[ missing in description ({})", description));
        let species_end = description.rfind(']').unwrap_or_else(|| panic!("] missing in description ({})", description));
        let species_name = &description[(species_start+1)..species_end];
        if taxonomy.contains_name(species_name) && taxonomy.is_descendant(species_name, ancestor_name) {
            output_fasta.write(record.id(), record.desc(), wrap(record.seq(), 80).as_slice()).unwrap();
        }
    }
}