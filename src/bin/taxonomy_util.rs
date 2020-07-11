#[macro_use]
extern crate clap;
extern crate ncbitaxonomy;

use std::path::Path;
use std::process;
use ncbitaxonomy::{NcbiTaxonomy, NcbiSqliteTaxonomy};

fn common_ancestor_distance(taxonomy: &dyn NcbiTaxonomy, name1: &str, name2: &str, only_canonical: bool) {
    match taxonomy.get_distance_to_common_ancestor(name1, name2, only_canonical) {
        Some((distance, common_ancestor_name)) => {
            println!("{}\t{}", distance, common_ancestor_name);
        },
        None => {
            eprintln!("no common ancestor found");
            process::exit(1)
        }
    }
}

pub fn main() {
    // TODO:
    // * write get_lineage - print lineage of taxon
    let app_m = clap_app!(taxonomy_util =>
        (version: ncbitaxonomy::VERSION)
        (author: "Peter van Heusden <pvh@sanbi.axc.za>")
        (about: "Utilities for working with the NCBI taxonomy database")
        (@arg TAXDB_URL: -d --db +takes_value "URL for SQLite taxonomy database")
        (@subcommand common_ancestor_distance =>
            (about: "find the tree distance to te common ancestor between two taxa")
            (@arg CANONICAL: --only_canonical "Only consider canonical taxonomic ranks")
            (@arg NAME1: +required "Name of first taxon")
            (@arg NAME2: +required "Name of second taxon")
        )
        (@subcommand get_id =>
            (about: "find taxonomy ID for name")
            (@arg NAME: +required "Name of taxon")
        )
        (@subcommand get_name =>
            (about: "find name for taxonomy ID")
            (@arg ID: +required "Taxonomy ID to look up")
        )
        (@subcommand get_lineage =>
            (about: "get lineage for name [unimplemented]")
            (@arg SHOW_NAMES: --show_names -S "Show taxon names, not just IDs")
            (@arg DELIMITER: --delimiter -D +takes_value "Delimiter for lineage string")
            (@arg NAME: +required "Name of taxon")
        )
        (@subcommand to_sqlite =>
            (about: "save taxonomy database loaded from files to SQLite database file")
            (@arg TAXONOMY_FILENAME_PREFIX: -t --tax_prefix +takes_value "String to prepend to names of nodes.dmp and names.dmp")
            (@arg TAXONOMY_DIR: +required "Directory containing the NCBI taxonomy nodes.dmp and names.dmp files")
        )
    ).get_matches();

    let taxdb_url = if app_m.is_present("TAXDB_URL") { Some(app_m.value_of("TAXDB_URL").unwrap()) } else { None };

    let taxonomy = NcbiSqliteTaxonomy::new(taxdb_url);

    match app_m.subcommand() {
        ("common_ancestor_distance", Some(sub_m)) => {
            let only_canonical = sub_m.is_present("CANONICAL");
            let name1 = sub_m.value_of("NAME1").unwrap();
            let name2 = sub_m.value_of("NAME2").unwrap();
            common_ancestor_distance(&taxonomy, name1, name2, only_canonical);
        },
        ("get_id", Some(sub_m)) => {
            let name = sub_m.value_of("NAME").unwrap();
            match taxonomy.get_id_by_name(name) {
                Some(val) => println!("{}", val),
                None => eprintln!("name {} not found in taxomomy", name)
            }
        },
        ("get_name", Some(sub_m)) => {
            let taxid = (sub_m.value_of("ID").unwrap()).parse::<i32>().unwrap();
            match taxonomy.get_name_by_id(taxid) {
                Some(val) => println!("{}", val),
                None => eprintln!("id {} not found in taxonomy", taxid)
            }
        },
        ("get_lineage", Some(sub_m)) => {
            let show_names = sub_m.is_present("SHOW_NAMES");
            let delimiter = sub_m.value_of("DELIMITER").unwrap_or(";");
            let name = sub_m.value_of("NAME").unwrap();

            match taxonomy.get_lineage(name) {
                None => eprintln!("{} not found in taxonomy", name),
                Some(lineage) => {
                    let output_list: Vec<String> = lineage.iter().map(|id| {
                        match show_names {
                            false => id.to_string(),
                            true => taxonomy.get_name_by_id(*id).unwrap() + " (" + id.to_string().as_str() + ")"
                        }
                    }).collect();
                    println!("{}", output_list.join(delimiter));
                }
            }
        }
        ("to_sqlite", Some(sub_m)) => {
            let ncbi_taxonomy_path = Path::new(sub_m.value_of("TAXONOMY_DIR").unwrap());

            let tax_prefix = match sub_m.value_of("TAXONOMY_FILENAME_PREFIX") {
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

            taxonomy.save_to_sqlite(taxdb_url).expect("failed to save taxonomy database to SQLite");
        },
        _ => {
            eprintln!("Unknown subcommand");
            eprintln!("{}", app_m.usage());
            process::exit(1);
        }
    }
}