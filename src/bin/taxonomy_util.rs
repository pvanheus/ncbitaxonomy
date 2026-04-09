extern crate ncbitaxonomy;

use clap::{Command, Arg};

use std::path::Path;
use std::process;
use ncbitaxonomy::{NcbiTaxonomy, NcbiSqliteTaxonomy};
use std::process::exit;

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
    let app = Command::new("taxonomy_util")
        .version(ncbitaxonomy::VERSION)
        .author("Peter van Heusden <pvh@sanbi.ac.za>")
        .about("Utilities for working with the NCBI taxonomy database")
        .arg(Arg::new("TAXDB_URL")
            .short('d')
            .long("db")
            .value_name("URL")
            .required(true)
            .help("URL for SQLite taxonomy database"))
        .subcommand(
            Command::new("common_ancestor_distance")
                .about("find the tree distance to the common ancestor between two taxa")
                .arg(Arg::new("CANONICAL")
                    .long("only_canonical")
                    .action(clap::ArgAction::SetTrue)
                    .help("Only consider canonical taxonomic ranks"))
                .arg(Arg::new("NAME1")
                    .required(true)
                    .help("Name of first taxon"))
                .arg(Arg::new("NAME2")
                    .required(true)
                    .help("Name of second taxon"))
        )
        .subcommand(
            Command::new("get_id")
                .about("find taxonomy ID for name")
                .arg(Arg::new("NAME")
                    .required(true)
                    .help("Name of taxon"))
        )
        .subcommand(
            Command::new("get_name")
                .about("find name for taxonomy ID")
                .arg(Arg::new("ID")
                    .required(true)
                    .help("Taxonomy ID to look up"))
        )
        .subcommand(
            Command::new("get_lineage")
                .about("get lineage for name")
                .arg(Arg::new("SHOW_NAMES")
                    .long("show_names")
                    .short('S')
                    .action(clap::ArgAction::SetTrue)
                    .help("Show taxon names, not just IDs"))
                .arg(Arg::new("DELIMITER")
                    .long("delimiter")
                    .short('D')
                    .value_name("DELIM")
                    .help("Delimiter for lineage string"))
                .arg(Arg::new("NAME")
                    .required(true)
                    .help("Name of taxon"))
        )
        .subcommand(
            Command::new("to_sqlite")
                .about("save taxonomy database loaded from files to SQLite database file")
                .arg(Arg::new("TAXONOMY_FILENAME_PREFIX")
                    .short('t')
                    .long("tax_prefix")
                    .value_name("PREFIX")
                    .help("String to prepend to names of nodes.dmp and names.dmp"))
                .arg(Arg::new("TAXONOMY_DIR")
                    .required(true)
                    .help("Directory containing the NCBI taxonomy nodes.dmp and names.dmp files"))
        );
    
    let matches = app.get_matches();

    let taxdb_url = matches.get_one::<String>("TAXDB_URL").expect("TAXDB_URL is required").as_str();
    // sqlite URLs are filename paths (the :memory: URL does not make sense for this application)
    match matches.subcommand_name() {
        Some("to_sqlite") => {}, // valid to have a nonexistent database here, we will write to it
        _ => {
            let db_path = Path::new(taxdb_url);
            if ! db_path.exists() {
                eprintln!("No database found at {}", taxdb_url);
                exit(1);
            }
        }
    }

    let taxonomy = NcbiSqliteTaxonomy::new(Some(taxdb_url));

    match matches.subcommand() {
        Some(("common_ancestor_distance", sub_m)) => {
            let only_canonical = sub_m.get_flag("CANONICAL");
            let name1 = sub_m.get_one::<String>("NAME1").expect("NAME1 is required").as_str();
            let name2 = sub_m.get_one::<String>("NAME2").expect("NAME2 is required").as_str();
            common_ancestor_distance(&taxonomy, name1, name2, only_canonical);
        },
        Some(("get_id", sub_m)) => {
            let name = sub_m.get_one::<String>("NAME").expect("NAME is required").as_str();
            match taxonomy.get_id_by_name(name) {
                Some(val) => println!("{}", val),
                None => eprintln!("name {} not found in taxonomy", name)
            }
        },
        Some(("get_name", sub_m)) => {
            let taxid = (sub_m.get_one::<String>("ID").expect("ID is required").as_str()).parse::<i32>().unwrap();
            match taxonomy.get_name_by_id(taxid) {
                Some(val) => println!("{}", val),
                None => eprintln!("id {} not found in taxonomy", taxid)
            }
        },
        Some(("get_lineage", sub_m)) => {
            let show_names = sub_m.get_flag("SHOW_NAMES");
            let delimiter = sub_m.get_one::<String>("DELIMITER").map(|s| s.as_str()).unwrap_or(";");
            let name = sub_m.get_one::<String>("NAME").expect("NAME is required").as_str();

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
        Some(("to_sqlite", sub_m)) => {
            let ncbi_taxonomy_path = Path::new(sub_m.get_one::<String>("TAXONOMY_DIR").expect("TAXONOMY_DIR is required").as_str());

            let tax_prefix = match sub_m.get_one::<String>("TAXONOMY_FILENAME_PREFIX") {
                Some(name) => name.to_string(),
                None => "".to_string()
            };

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

            taxonomy.save_to_sqlite(Some(taxdb_url)).expect("failed to save taxonomy database to SQLite");
        },
        _ => {
            eprintln!("Unknown subcommand");
            process::exit(1);
        }
    }
}