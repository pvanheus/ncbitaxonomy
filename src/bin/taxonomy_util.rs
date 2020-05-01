#[macro_use]
extern crate clap;
extern crate ncbitaxonomy;

use std::path::Path;
use std::process;
use ncbitaxonomy::{NcbiTaxonomy, NcbiSqliteTaxonomy};

fn common_ancestor_distance(taxonomy: &dyn NcbiTaxonomy, name1: &str, name2: &str, only_canonical: bool) {
    match taxonomy.get_distance_to_common_ancestor(name1, name2, only_canonical) {
        Some(distance) => {
            println!("{}", distance);
        },
        None => {
            eprintln!("no common ancestor found");
            process::exit(1)
        }
    }
}

pub fn main() {
    let app_m = clap_app!(taxonomy_filter_refseq =>
        (version: ncbitaxonomy::VERSION)
        (author: "Peter van Heusden <pvh@sanbi.axc.za>")
        (about: "Utilities for working with the NCBI taxonomy database")
        (@arg TAXONOMY_FILENAME_PREFIX: -t --tax_prefix +takes_value "String to prepend to names of nodes.dmp and names.dmp")
        (@arg TAXONOMY_DIR: +required "Directory containing the NCBI taxonomy nodes.dmp and names.dmp files")
        (@subcommand common_ancestor_distance =>
            (about: "find the tree distance to te common ancestor between two taxa")
            (@arg CANONICAL: --only_canonical "Only consider canonical taxonomic ranks")
            (@arg NAME1: +required "Name of first taxon")
            (@arg NAME2: +required "Name of second taxon")
        )
        (@subcommand find_id =>
            (about: "find taxonomy ID for name")
            (@arg NAME: +required "Name of taxon")
        )
        (@subcommand find_name =>
            (about: "find name for taxonomy ID")
            (@arg ID: +required "Taxonomy ID to look up")
        )
        (@subcommand sqlite =>
            (about: "sqlite testing")
        )
    ).get_matches();

    let ncbi_taxonomy_path = Path::new(app_m.value_of("TAXONOMY_DIR").unwrap());

    let tax_prefix = match app_m.value_of("TAXONOMY_FILENAME_PREFIX") {
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

    eprintln!("loading taxonomy");
    let taxonomy = ncbitaxonomy::NcbiFileTaxonomy::from_ncbi_files(
        nodes_path.as_path().to_str().unwrap(),
        names_path.as_path().to_str().unwrap()).expect("Failed to load NCBI Taxonomy");
    eprintln!("taxonomy loaded");

    match app_m.subcommand() {
        ("common_ancestor_distance", Some(sub_m)) => {
            let only_canonical = sub_m.is_present("CANONICAL");
            let name1 = sub_m.value_of("NAME1").unwrap();
            let name2 = sub_m.value_of("NAME2").unwrap();
            common_ancestor_distance(&taxonomy, name1, name2, only_canonical);
        },
        ("find_id", Some(sub_m)) => {
            let taxid = (sub_m.value_of("NAME").unwrap()).parse::<i32>().unwrap();
            let taxonomy = NcbiSqliteTaxonomy::new(None);
            println!("{}", taxonomy.contains_id(taxid))

        },
        ("find_name", Some(sub_m)) => {

        },
        ("sqlite", Some(sub_m)) => {
            taxonomy.save_to_sqlite().expect("failed to save taxonomy database to SQLite");
        },
        _ => {
            eprintln!("Unknown subcommand");
            eprintln!("{}", app_m.usage());
            process::exit(1);
        }
    }
}