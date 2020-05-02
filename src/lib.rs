#![recursion_limit = "1024"]
#[macro_use]
extern crate diesel;
extern crate dotenv;
extern crate thiserror;
extern crate indextree;
extern crate core;
extern crate seq_io;
extern crate clap;

/// ncbitaxonomy: a module for working with a local copy of the NCBI taxonomy database

use thiserror::Error;
use std::io;

#[derive(Error, Debug)]
pub enum NcbiTaxonomyError {
    #[error("I/O error opening or reading file")]
    Io(#[from] io::Error),
    #[error("format error in nodes.dmp in line {0}")]
    NodeFileFormatError(String),
    #[error("failed to parse integer from string {0}")]
    ParseIntError(#[from] ::std::num::ParseIntError)
}

use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufReader,BufRead};
use indextree::{Arena, NodeId, Traverse};
pub use indextree::NodeEdge;
use std::iter::FromIterator;

pub mod models;
pub mod schema;

use diesel::prelude::*;
use diesel::result::Error as DieselError;
use diesel::sqlite::SqliteConnection;
use dotenv::dotenv;
use std::env;

use self::models::*;
use diesel::expression::dsl::count;

pub const VERSION: &str = env!("CARGO_PKG_VERSION");

fn establish_connection() -> SqliteConnection {
    dotenv().ok();

    let database_url = env::var("DATABASE_URL").expect("DATABASE_URL must be set");
    SqliteConnection::establish(&database_url)
        .unwrap_or_else(|_| panic!("Error connecting to {}", database_url))
}

fn get_canonical_ranks() -> HashSet<String> {
    // canonical ranks (+ superkingdom) as they appear in the NCBI taxonomy database
    HashSet::from_iter(vec!["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"].iter().map(|x| (*x).to_string()))
}

pub trait NcbiTaxonomy {
    fn contains_id(&self, taxid: i32) -> bool;
    fn contains_name(&self, name: &str) -> bool;
    fn is_descendant(&self, name: &str, ancestor_name: &str) -> bool;
    fn is_descendant_taxid(&self, taxid: i32, ancestor_taxid: i32) -> bool;
    fn get_name_by_id(&self, taxid: i32) -> Option<String>;
    fn get_id_by_name(&self, name: &str) -> Option<i32>;
    fn get_distance_to_common_ancestor_taxid(&self, taxid1: i32, taxid2: i32, only_canonical: bool) -> Option<i32>;
    fn get_distance_to_common_ancestor(&self, name1: &str, name2: &str, only_canonical: bool) -> Option<i32>;
}

#[derive(Debug)]
pub struct NcbiFileTaxonomy {
    arena: Arena<i32>,
    name_to_node: HashMap<String, NodeId>,
    id_to_node: HashMap<i32, NodeId>,
    id_to_name: HashMap<i32, String>,
    id_to_rank: HashMap<i32, String>
}

impl NcbiFileTaxonomy {

    /// from_ncbi_files
    ///
    /// Reads the `nodes.dmp` file and `names.dmp` file from the NCBI Taxonomy database to
    /// generate a NcbiTaxonomy structure
    ///
    /// # Examples
    ///
    /// ```
    /// use ncbitaxonomy::*;
    ///
    /// let taxonomy = NcbiFileTaxonomy::from_ncbi_files("data/nodes.dmp", "data/names.dmp");
    /// ```
    pub fn from_ncbi_files(nodes_filename: &str, names_filename: &str) -> Result<NcbiFileTaxonomy, NcbiTaxonomyError> {
        let mut child_ids_by_parent_id: HashMap<i32, Vec<i32>> = HashMap::new();
        let mut id_to_rank = HashMap::new();
        let nodes_file = File::open(nodes_filename)?;
        for line_maybe in BufReader::new(nodes_file).lines() {
            let line = line_maybe?;
            let mut fields = line.split("\t|\t");
            let id_str = fields.next().ok_or_else(|| NcbiTaxonomyError::NodeFileFormatError(line.clone()))?;
            let parent_id_str = fields.next().ok_or_else(|| NcbiTaxonomyError::NodeFileFormatError(line.clone()))?;
            let rank = fields.next().ok_or_else(|| NcbiTaxonomyError::NodeFileFormatError(line.clone()))?.to_string();
            let id = id_str.parse::<i32>().or_else(|e| Err(NcbiTaxonomyError::ParseIntError(e)))?;
            let parent_id = parent_id_str.parse::<i32>().or_else(|e| Err(NcbiTaxonomyError::ParseIntError(e)))?;
            id_to_rank.insert(id, rank);
            if parent_id != id {  // this happens for the root node
                // thanks to https://stackoverflow.com/questions/33243784/append-to-vector-as-value-of-hashmap/33243862
                // for this way to get the existing entry or insert an empty list.
                child_ids_by_parent_id.entry(parent_id).or_insert_with(Vec::new).push(id);
            }
        }

        let mut keys = child_ids_by_parent_id.keys().collect::<Vec<&i32>>();
        keys.sort_unstable_by(|a, b| a.cmp(b));

        let mut arena: Arena<i32> = Arena::new();

        let mut id_to_node: HashMap<i32, NodeId> = HashMap::new();
        for id in keys {
            let node_id = match id_to_node.get(id) {
                Some(node_id) => *node_id,
                None => arena.new_node(*id),
            };
            id_to_node.insert(*id, node_id);
            for child in child_ids_by_parent_id.get(&id).expect("ID not found in child_ids_by_parent_id") {
                let child_node_id = match id_to_node.get(child) {
                    Some(child_node_id) => *child_node_id,
                    None => arena.new_node(*child),
                };
                id_to_node.insert(*child, child_node_id);
                assert_ne!(node_id, child_node_id, "child node id same as node: {} (for {} {})", node_id, *id, *child);
                node_id.append(child_node_id, &mut arena).unwrap();  // might return Failure, in which case we panic!
            }
        }

        // now its time to read the names_filename that maps names to IDs
        let mut name_to_node = HashMap::new();
        let mut id_to_name = HashMap::new();
        let name_file = File::open(names_filename)?;
        for line_maybe in BufReader::new(name_file).lines() {
            let line = line_maybe?;
            let fields = line.split("\t|\t").collect::<Vec<&str>>();
            if fields[3].starts_with("scientific name") {
                let id_str = fields[0];
                let id = id_str.parse::<i32>().or_else(|e| Err(NcbiTaxonomyError::ParseIntError(e)))?;
                let name = if fields[2] != "" { fields[2].to_string() } else { fields[1].to_string() };
                let node_id = id_to_node.get(&id).expect("ID not found in id_to_node");
                id_to_name.insert(id, name.clone());
                name_to_node.insert(name, *node_id);
            }
        }

        let tree = NcbiFileTaxonomy { arena, name_to_node, id_to_node, id_to_name, id_to_rank };
        Ok(tree)
    }

    pub fn save_to_sqlite(&self) -> Result<(), DieselError> {
        // design of storing a tree in a relational DB inspired by:
        // https://makandracards.com/makandra/45275-storing-trees-in-databases
        use schema::taxonomy;
        let connection = establish_connection();

        connection.transaction::<_, DieselError, _>(|| {
            for (id, nodeid) in self.id_to_node.iter() {
                let mut ancestors_vec = nodeid.ancestors(&self.arena).map(|nodeid| self.get_id_by_node(nodeid).unwrap().to_string()).collect::<Vec<String>>();
                ancestors_vec.reverse();
                let ancestors_string = ancestors_vec.join("/");
                let name = self.id_to_name.get(id).unwrap();
                let taxon_record = NewTaxon {
                    id,
                    ancestry: match ancestors_string  {
                        v if v == "1" => None,
                        _ => Some(&ancestors_string[..])
                    },
                    name,
                    rank: match self.id_to_rank.get(id) {
                        Some(v) => Some(&v[..]),
                        None => None
                    }

                };
                diesel::insert_into(taxonomy::table)
                    .values(&taxon_record   )
                    .execute(&connection)?;
            }
            Ok(())
        })?;
        Ok(())
    }

    /// get_node_by_id
    ///
    /// get a NodeId from a numeric NCBI Taxonomy ID
    pub fn get_node_by_id(&self, id: i32) -> Option<&NodeId> {
        self.id_to_node.get(&id)
    }

    /// traversal
    ///
    /// traverse the tree nodes (in depth first order) from the node with a given NCBI Taxonomy ID
    pub fn traversal(&self, from: i32) -> Option<Traverse<i32>> {
        match self.get_node_by_id(from) {
            Some(node_id) => Some(node_id.traverse(&self.arena)),
            None => None
        }
    }

    /// get_id_by_node
    ///
    /// get the NCBI Taxonomy ID held by the node with a given NodeId
    pub fn get_id_by_node(&self, node_id: NodeId) -> Option<i32> {
        match self.arena.get(node_id)  {
            Some(node) => Some(node.data),
            None => None
        }
    }

    // TODO write tests for get_distance_to_common_ancestor and get_distance_to_common_ancestor_id
}

impl NcbiTaxonomy for NcbiFileTaxonomy {
    /// contains_id
    ///
    /// check whether the taxonomy contains a (number) ID
    fn contains_id(&self, id: i32) -> bool {
        self.id_to_node.contains_key(&id)
    }

    /// contains_name
    ///
    /// check whether the taxonomy contains a node with the specified name
    ///
    /// **note:** the name used is what is reported as a the 'scientific name' in the NCBI Taxonomy database.
    /// synonyms are currently not supported
    fn contains_name(&self, name: &str) -> bool {
        self.name_to_node.contains_key(name)
    }

    /// is_descendant
    ///
    /// check if a certain named node is a descendant of another named named
    fn is_descendant(&self, name: &str, ancestor_name: &str) -> bool {
        let id = match self.get_id_by_name(name) {
            Some(id) => id,
            None => return false
        };
        let ancestor_id = match self.get_id_by_name(ancestor_name) {
            Some(id) => id,
            None => return false
        };
        self.is_descendant_taxid(id, ancestor_id)
    }

    /// is_descendant_taxid
    ///
    /// check if a certain node with taxid is a descendant of another taxid
    fn is_descendant_taxid(&self, taxid: i32, ancestor_taxid: i32) -> bool {
        let id = match self.id_to_node.get(&taxid) {
            Some(id) => id,
            None => return false
        };
        let ancestor_id = match self.id_to_node.get(&ancestor_taxid) {
            Some(id) => id,
            None => return false
        };
        for node in id.ancestors(&self.arena) {
            if node == *ancestor_id {
                return true
            }
        }
        false
    }

    /// get_name_by_id
    ///
    /// get the scientific name associated with a given NCBI Taxonomy ID
    fn get_name_by_id(&self, id: i32) -> Option<String> {
        self.id_to_name.get(&id).cloned()
    }

    fn get_id_by_name(&self, name: &str) -> Option<i32> {
        match self.name_to_node.get(name) {
            Some(nodeid) => self.get_id_by_node(*nodeid),
            None => None
        }
    }

    /// get_distance_to_common_ancestor_id
    ///
    /// get the distance (in steps in the tree) between taxid1 and the common ancestor with taxid2
    fn get_distance_to_common_ancestor_taxid(&self, taxid1: i32, taxid2: i32, only_canonical: bool) -> Option<i32> {
        let canonical_ranks = get_canonical_ranks();
        if taxid1 == taxid2 {
            return Some(0)
        }

        let taxon1 = self.id_to_node.get(&taxid1)?;
        let taxon2 = self.id_to_node.get(&taxid2)?;

        let mut ancestors_distance1 = HashMap::new();
        let mut current_distance = 0;
        let taxid1_rank = self.id_to_rank.get(&taxid1)?;
        if !only_canonical || canonical_ranks.contains(taxid1_rank) {
            // we know that taxid1 != taxid2, so either taxid2 is an ancestor of
            // taxid1 or there is a common ancestor further back. in the first case,
            // taxid2 will be found in the ancestors of taxid1, so thifs distance will
            // not be used. in the second case, taxid1 might be found in the ancestors
            // of taxid2, so this distance will still not be used
            ancestors_distance1.insert(taxid1, 0);
        }
        for node in taxon1.ancestors(&self.arena) {
            let nodeid = self.get_id_by_node(node)?;
            let rank = self.id_to_rank.get(&nodeid)?;
            if !only_canonical || canonical_ranks.contains(rank) {
                current_distance += 1;
                if nodeid == taxid2 {
                    return Some(current_distance)
                }
                ancestors_distance1.insert(nodeid, current_distance);
            }
        }

        // if we got here then we did not see taxid2 yet, i.e. taxid2 is not a
        // direct ancestor of taxid1
        current_distance = 0;
        for node in taxon2.ancestors(&self.arena) {
            let nodeid = self.get_id_by_node(node).unwrap();
            let rank = self.id_to_rank.get(&nodeid)?;
            if !only_canonical || canonical_ranks.contains(rank) {
                current_distance += 1;
                if ancestors_distance1.contains_key(&nodeid) {
                    // the distance to te common ancestor is the distance from taxon2
                    // to an ancestor that is also an ancestor to taxon1
                    return Some(current_distance)
                }
            }
        }
        None
    }

    /// get_distance_to_common_ancestor
    ///
    /// find the distance in the tree between name1 and name2
    fn get_distance_to_common_ancestor(&self, name1: &str, name2: &str, only_canonical: bool) -> Option<i32> {
        let taxon1 = self.name_to_node.get(name1)?;

        let taxon2 = self.name_to_node.get(name2)?;

        self.get_distance_to_common_ancestor_taxid(self.get_id_by_node(*taxon1).unwrap(),
                                                   self.get_id_by_node(*taxon2).unwrap(), only_canonical)
    }
}

pub struct NcbiSqliteTaxonomy {
    connection: SqliteConnection
}

impl NcbiSqliteTaxonomy {
    pub fn new(db_url: Option<&str>) -> Self {
        dotenv().ok();

        let database_url = match db_url {
            Some(database_url) => database_url.to_string(),
            None => env::var("DATABASE_URL").expect("DATABASE_URL must be set")
        };
        NcbiSqliteTaxonomy {
            connection: SqliteConnection::establish(&database_url).unwrap_or_else(|_| panic!("Error connecting to {}", database_url))
        }
    }

    fn get_ancestry_for_taxid(&self, taxid: i32) -> Option<String> {
        use schema::taxonomy::dsl::*;

        let results: Vec<Option<String>> = taxonomy.filter(id.eq(taxid))
            .select(ancestry)
            .load(&self.connection)
            .expect("Error loading taxonomy");

        match results.len() {
            1 => results[0].clone(),
            _ => panic!("taxid {} not found in taxonomy", taxid)
        }
    }

    fn get_ancestors(&self, taxid: i32) -> Vec<i32> {
        let ancestry_string = self.get_ancestry_for_taxid(taxid);
        match ancestry_string {
            None => vec![], // the root taxon has no ancestry
            Some(val) => {
                let mut ancestors: Vec<i32> = val.split('/').map(|id_str| id_str.parse::<i32>().unwrap()).collect();
                ancestors.reverse();
                ancestors
            }
        }
    }

    fn get_rank(&self, taxid: i32) -> Option<String> {
        use schema::taxonomy::dsl::*;

        let results: Vec<Option<String>> = taxonomy.filter(id.eq(taxid))
            .select(rank)
            .load(&self.connection)
            .expect("Error loading taxonomy");

        match results.len() {
            1 => results[0].clone(),
            _ => panic!("taxid {} not found in taxonomy", taxid)
        }
    }
}

impl NcbiTaxonomy for NcbiSqliteTaxonomy {

    fn contains_id(&self, taxid: i32) -> bool {
        use schema::taxonomy::dsl::*;

        let results: Vec<i64> = taxonomy.filter(id.eq(taxid))
            .select(count(id))
            .load(&self.connection)
            .expect("Error loading taxonomy");

        results[0] == 1
    }

    fn contains_name(&self, name_str: &str) -> bool {
        use schema::taxonomy::dsl::*;

        let results: Vec<i64> = taxonomy.filter(name.eq(name_str))
            .select(count(id))
            .load(&self.connection)
            .expect("Error loading taxonomy");

        results[0] == 1
    }

    fn is_descendant(&self, name_str: &str, ancestor: &str) -> bool {
        let taxid = match self.get_id_by_name(name_str) {
            Some(val) => val,
            None => return false
        };

        let ancestor_taxid = match self.get_id_by_name(ancestor) {
            Some(val) => val,
            None => return false
        };

        self.is_descendant_taxid(taxid, ancestor_taxid)
    }

    fn is_descendant_taxid(&self, taxid: i32, ancestor_taxid: i32) -> bool {
        use schema::taxonomy::dsl::*;

        // ancestor pattern is id/id/id so if ancestor_taxid is an ancestor
        // of taxid, LIKE 'ancestor_taxid/%' OR LIKE '%/ancestor_taxid/%' OR LIKE '%/ancestor_taxid'
        // will be true
        let pattern1 = format!("{}/%", ancestor_taxid);
        let pattern2 = format!("%/{}/%", ancestor_taxid);
        let pattern3 = format!("%/{}", ancestor_taxid);

        let results: Vec<i64> = taxonomy.filter(
                id.eq(taxid).and(
                    ancestry.like(pattern1)
                        .or(ancestry.like(pattern2))
                        .or(ancestry.like(pattern3))
                ))
            .select(count(id))
            .load(&self.connection)
            .expect("Error loading taxonomy");

        match results.len() {
            1 => true,
            _ => false
        }
    }

    fn get_name_by_id(&self, taxid: i32) -> Option<String> {
        use schema::taxonomy::dsl::*;

        let results: Vec<String> = taxonomy.filter(id.eq(taxid))
            .select(name)
            .load(&self.connection)
            .expect("Error loading taxonomy");

        match results.len() {
            1 => Some(results[0].clone()),
            _ => None
        }
    }

    fn get_id_by_name(&self, name_str: &str) -> Option<i32> {
        use schema::taxonomy::dsl::*;

        let results: Vec<i32> = taxonomy.filter(name.eq(name_str))
            .select(id)
            .load(&self.connection)
            .expect("Error loading taxonomy");

        match results.len() {
            1 => Some(results[0]),
            _ => None
        }
    }

    fn get_distance_to_common_ancestor_taxid(&self, taxid1: i32, taxid2: i32, only_canonical: bool) -> Option<i32> {
        // canonical ranks (+ superkingdom) as they appear in the NCBI taxonomy database
        let canonical_ranks = get_canonical_ranks();

        if taxid1 == taxid2 {
            return Some(0)
        }

        let mut ancestors_distance1 = HashMap::new();
        let mut current_distance = 0;
        // TODO: make rank a NON NULL column
        let taxid1_rank = self.get_rank(taxid1)?;
        if !only_canonical || canonical_ranks.contains(&taxid1_rank) {
            // see comment above for why distance is 0
            ancestors_distance1.insert(taxid1, 0);
        }
        for taxid in self.get_ancestors(taxid1) {
            let current_rank = self.get_rank(taxid)?;
            if taxid == taxid2 {
                return Some(current_distance)
            }
            if only_canonical || canonical_ranks.contains(&current_rank) {
                current_distance += 1;
                ancestors_distance1.insert(taxid, current_distance);
            }
        }

        current_distance = 0;
        for taxid in self.get_ancestors(taxid2) {
            let current_rank = self.get_rank(taxid)?;
            eprintln!("{}", self.get_name_by_id(taxid).unwrap());
            if !only_canonical || canonical_ranks.contains(&current_rank) {
                current_distance += 1;
                if ancestors_distance1.contains_key(&taxid) {
                    return Some(current_distance)
                }
            }
        }
        None
    }

    fn get_distance_to_common_ancestor(&self, name1: &str, name2: &str, only_canonical: bool) -> Option<i32> {
        let taxid1 = match self.get_id_by_name(name1) {
            Some(val) => val,
            None => return None
        };

        let taxid2 = match self.get_id_by_name(name2) {
            Some(val) => val,
            None => return None
        };

        self.get_distance_to_common_ancestor_taxid(taxid1, taxid2, only_canonical)
    }
}

#[cfg(test)]
mod tests {
    use super::{NcbiFileTaxonomy, NcbiSqliteTaxonomy, NcbiTaxonomy, NodeEdge};

    pub struct NcbiFileTaxonomyFixture {
        pub taxonomy: NcbiFileTaxonomy,
    }

    impl Default for NcbiFileTaxonomyFixture {
        fn default() -> Self {
            let tree = NcbiFileTaxonomy::from_ncbi_files("data/sample_tree_nodes.dmp", "data/sample_tree_names.dmp").unwrap();
            Self { taxonomy: tree }
        }
    }

    pub struct NcbiSqliteTaxonomyFixture {
        pub taxonomy: NcbiSqliteTaxonomy,
    }


    // TODO: built this in memory on the fly
    impl Default for NcbiSqliteTaxonomyFixture {
        fn default() -> Self {
            let tree = NcbiSqliteTaxonomy::new(Some("data/ncbi_taxonomy.sqlite"));
            Self { taxonomy: tree }
        }
    }

    #[test]
    fn contains_id() {
        let fixture = NcbiFileTaxonomyFixture::default();
        assert!(fixture.taxonomy.contains_id(504556));
    }

    #[test]
    fn sqlite_contains_id() {
        let fixture = NcbiSqliteTaxonomyFixture::default();
        assert!(fixture.taxonomy.contains_id(504556));
    }

    #[test]
    fn contains_name() {
        let fixture = NcbiFileTaxonomyFixture::default();
        assert!(fixture.taxonomy.contains_name("Propionibacterium phage PAS7"));
        assert!(fixture.taxonomy.contains_name("Viruses"));
        assert!(fixture.taxonomy.contains_name("environmental samples <bacteriophages>"));
    }

    #[test]
    fn sqlite_contains_name() {
        let fixture = NcbiSqliteTaxonomyFixture::default();
        assert!(fixture.taxonomy.contains_name("Propionibacterium phage PAS7"));
        assert!(fixture.taxonomy.contains_name("Viruses"));
        assert!(fixture.taxonomy.contains_name("environmental samples <bacteriophages>"));
    }

    #[test]
    fn get_node_by_id() {
        let fixture = NcbiFileTaxonomyFixture::default();
        assert_eq!(fixture.taxonomy.get_node_by_id(999999999), None);
        assert!(match fixture.taxonomy.get_node_by_id(504556) { Some(_) => true, None => false })
    }

    #[test]
    fn get_id_by_node() {
        let fixture = NcbiFileTaxonomyFixture::default();
        let node_id = fixture.taxonomy.get_node_by_id(504556).unwrap();
        assert_eq!(fixture.taxonomy.get_id_by_node(*node_id).unwrap(), 504556);
    }

    #[test]
    fn get_name_by_id() {
        let fixture = NcbiFileTaxonomyFixture::default();
        assert_eq!(fixture.taxonomy.get_name_by_id(370556).unwrap(), "Streptococcus phage 9429.1");
    }

    #[test]
    fn sqlite_get_name_by_id() {
        let fixture = NcbiSqliteTaxonomyFixture::default();
        assert_eq!(fixture.taxonomy.get_name_by_id(370556).unwrap(), "Streptococcus phage 9429.1");
    }

    #[test]
    fn traversal() {
        let fixture = NcbiFileTaxonomyFixture::default();
        let traversal = fixture.taxonomy.traversal(12333);
        match traversal {
            Some(traversal) => {
                let mut counter = 0;
                for node_edge in traversal {
                    match node_edge {
                        NodeEdge::Start(_) => counter += 1,
                        _ => ()
                    }
                }
                assert_eq!(counter, 500)
            }
            None => assert_eq!(true, false, "Failed to load traversal from 12333")
        }
    }

    #[test]
    fn descendants() {
        let fixture = NcbiFileTaxonomyFixture::default();
        assert!(fixture.taxonomy.is_descendant("Propionibacterium phage PAS7", "unclassified bacterial viruses"));
    }

    #[test]
    fn sqlite_descendants() {
        let fixture = NcbiSqliteTaxonomyFixture::default();
        assert!(fixture.taxonomy.is_descendant("Propionibacterium phage PAS7", "unclassified bacterial viruses"));
    }

    #[test]
    fn taxid_descendants() {
        let fixture = NcbiFileTaxonomyFixture::default();
        assert!(fixture.taxonomy.is_descendant_taxid(504556, 12333));
    }

    #[test]
    fn sqlite_taxid_descendants() {
        let fixture = NcbiSqliteTaxonomyFixture::default();
        assert!(fixture.taxonomy.is_descendant_taxid(504556, 12333));
    }
}
