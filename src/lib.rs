#![recursion_limit = "1024"]

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

pub const VERSION: &str = env!("CARGO_PKG_VERSION");

#[derive(Debug)]
pub struct NcbiTaxonomy {
    arena: Arena<u32>,
    name_to_node: HashMap<String, NodeId>,
    id_to_node: HashMap<u32, NodeId>,
    id_to_name: HashMap<u32, String>,
    id_to_rank: HashMap<u32, String>
}

impl NcbiTaxonomy {

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
    /// let taxonomy = NcbiTaxonomy::from_ncbi_files("data/nodes.dmp", "data/names.dmp");
    /// ```
    pub fn from_ncbi_files(nodes_filename: &str, names_filename: &str) -> Result<NcbiTaxonomy, NcbiTaxonomyError> {
        let mut child_ids_by_parent_id: HashMap<u32, Vec<u32>> = HashMap::new();
        let mut id_to_rank = HashMap::new();
        let nodes_file = File::open(nodes_filename)?;
        for line_maybe in BufReader::new(nodes_file).lines() {
            let line = line_maybe?;
            let mut fields = line.split("\t|\t");
            let id_str = fields.next().ok_or_else(|| NcbiTaxonomyError::NodeFileFormatError(line.clone()))?;
            let parent_id_str = fields.next().ok_or_else(|| NcbiTaxonomyError::NodeFileFormatError(line.clone()))?;
            let rank = fields.next().ok_or_else(|| NcbiTaxonomyError::NodeFileFormatError(line.clone()))?.to_string();
            let id = id_str.parse::<u32>().or_else(|e| Err(NcbiTaxonomyError::ParseIntError(e)))?;
            let parent_id = parent_id_str.parse::<u32>().or_else(|e| Err(NcbiTaxonomyError::ParseIntError(e)))?;
            id_to_rank.insert(id, rank);
            if parent_id != id {  // this happens for the root node
                // thanks to https://stackoverflow.com/questions/33243784/append-to-vector-as-value-of-hashmap/33243862
                // for this way to get the existing entry or insert an empty list.
                child_ids_by_parent_id.entry(parent_id).or_insert_with(Vec::new).push(id);
            }
        }

        let mut keys = child_ids_by_parent_id.keys().collect::<Vec<&u32>>();
        keys.sort_unstable_by(|a, b| a.cmp(b));

        let mut arena: Arena<u32> = Arena::new();

        let mut id_to_node: HashMap<u32, NodeId> = HashMap::new();
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
                let id = id_str.parse::<u32>().or_else(|e| Err(NcbiTaxonomyError::ParseIntError(e)))?;
                let name = fields[1].to_string();
                let node_id = id_to_node.get(&id).expect("ID not found in id_to_node");
                id_to_name.insert(id, name.clone());
                name_to_node.insert(name, *node_id);
            }
        }

        let tree = NcbiTaxonomy { arena, name_to_node, id_to_node, id_to_name, id_to_rank };
        Ok(tree)
    }

    /// contains_id
    ///
    /// check whether the taxonomy contains a (number) ID
    pub fn contains_id(&self, id: u32) -> bool {
        self.id_to_node.contains_key(&id)
    }

    /// contains_name
    ///
    /// check whether the taxonomy contains a node with the specified name
    ///
    /// **note:** the name used is what is reported as a the 'scientific name' in the NCBI Taxonomy database.
    /// synonyms are currently not supported
    pub fn contains_name(&self, name: &str) -> bool {
        self.name_to_node.contains_key(name)
    }

    /// is_descendant
    ///
    /// check if a certain named node is a descendant of another named named
    pub fn is_descendant(&self, name: &str, ancestor_name: &str) -> bool {
        let id = match self.name_to_node.get(name) {
            Some(id) => id,
            None => return false
        };
        let ancestor_id = match self.name_to_node.get(ancestor_name) {
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

    /// is_descendant
    ///
    /// check if a certain node with taxid is a descendant of another taxid
    pub fn is_descendant_taxid(&self, taxid: u32, ancestor_taxid: u32) -> bool {
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

    /// get_node_by_id
    ///
    /// get a NodeId from a numeric NCBI Taxonomy ID
    pub fn get_node_by_id(&self, id: u32) -> Option<&NodeId> {
        self.id_to_node.get(&id)
    }

    /// traversal
    ///
    /// traverse the tree nodes (in depth first order) from the node with a given NCBI Taxonomy ID
    pub fn traversal(&self, from: u32) -> Option<Traverse<u32>> {
        match self.get_node_by_id(from) {
            Some(node_id) => Some(node_id.traverse(&self.arena)),
            None => None
        }
    }

    /// get_id_by_node
    ///
    /// get the NCBI Taxonomy ID held by the node with a given NodeId
    pub fn get_id_by_node(&self, node_id: NodeId) -> Option<u32> {
        match self.arena.get(node_id)  {
            Some(node) => Some(node.data),
            None => None
        }
    }

    /// get_name_by_id
    ///
    /// get the scientific name associated with a given NCBI Taxonomy ID
    pub fn get_name_by_id(&self, id: u32) -> Option<&String> {
        self.id_to_name.get(&id)
    }

    /// get_distance_to_common_ancestor_id
    ///
    /// get the distance (in steps in the tree) between taxid1 and the common ancestor with taxid2
    pub fn get_distance_to_common_ancestor_id(&self, taxid1: u32, taxid2: u32, only_canonical: bool) -> Option<u32> {
        // canonical ranks (+ superkingdom) as they appear in the NCBI taxonomy database
        let canonical_ranks: HashSet<String>  = HashSet::from_iter(vec!["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"].iter().map(|x| x.to_string()));
        if taxid1 == taxid2 {
            return Some(0)
        }

        let taxon1 = self.id_to_node.get(&taxid1)?;
        let taxon2 = self.id_to_node.get(&taxid2)?;

        let mut ancestors_distance1 = HashMap::new();
        let mut current_distance = 0;
        ancestors_distance1.insert(taxid1, current_distance);
        for node in taxon1.ancestors(&self.arena) {
            let nodeid = self.get_id_by_node(node)?;
            let rank = self.id_to_rank.get(&nodeid)?;
            if !only_canonical || canonical_ranks.contains(rank) {
                current_distance += 1;
                ancestors_distance1.insert(self.get_id_by_node(node).unwrap(), current_distance);
            }
        }

        // if the ancestors_distance1 map contains taxid2
        // then taxon2 is an ancestor of taxon1 and we return
        // the distnace between taxon1 and this ancestor
        if ancestors_distance1.contains_key(&taxid2) {
            return Some(*ancestors_distance1.get(&taxid2).unwrap());
        }

        current_distance = 0;
        for node in taxon2.ancestors(&self.arena) {
            let nodeid = self.get_id_by_node(node).unwrap();
            let rank = self.id_to_rank.get(&nodeid)?;
            if !only_canonical || canonical_ranks.contains(rank) {
                current_distance += 1;
                if ancestors_distance1.contains_key(&nodeid) {
                    // the distance to te common ancestor is the distance from taxon2
                    // to an ancestor that is also an ancestor to taxon1
                    // eprintln!("common ancestor {} {} {} {} {}", self.id_to_name.get(&nodeid)?, nodeid, rank, canonical_ranks.contains(rank), current_distance);
                    return Some(current_distance)
                }
            }
        }
        None
    }

    /// get_distance_to_common_ancestor
    ///
    /// find the distance in the tree between name1 and name2
    pub fn get_distance_to_common_ancestor(&self, name1: &str, name2: &str, only_canonical: bool) -> Option<u32> {
        let taxon1 = self.name_to_node.get(name1)?;

        let taxon2 = self.name_to_node.get(name2)?;

        self.get_distance_to_common_ancestor_id(self.get_id_by_node(*taxon1).unwrap(),
                                                self.get_id_by_node(*taxon2).unwrap(), only_canonical)
    }

    // TODO write tests for get_distance_to_common_ancestor and get_distance_to_common_ancestor_id
}

#[cfg(test)]
mod tests {
    use super::{NcbiTaxonomy, NodeEdge};

    pub struct NcbiTaxonomyFixture {
        pub taxonomy: NcbiTaxonomy,
    }

    impl Default for NcbiTaxonomyFixture {
        fn default() -> Self {
            let tree = NcbiTaxonomy::from_ncbi_files("data/sample_tree_nodes.dmp", "data/sample_tree_names.dmp").unwrap();
            Self { taxonomy: tree }
        }
    }

    #[test]
    fn contains_id() {
        let fixture = NcbiTaxonomyFixture::default();
        assert!(fixture.taxonomy.contains_id(504556));
    }

    #[test]
    fn contains_name() {
        let fixture = NcbiTaxonomyFixture::default();
        assert!(fixture.taxonomy.contains_name("Propionibacterium phage PAS7"));
        assert!(fixture.taxonomy.contains_name("Viruses"));
        assert!(fixture.taxonomy.contains_name("environmental samples"));
    }


    #[test]
    fn get_node_by_id() {
        let fixture = NcbiTaxonomyFixture::default();
        assert_eq!(fixture.taxonomy.get_node_by_id(999999999), None);
        assert!(match fixture.taxonomy.get_node_by_id(504556) { Some(_) => true, None => false })
    }

    #[test]
    fn get_id_by_node() {
        let fixture = NcbiTaxonomyFixture::default();
        let node_id = fixture.taxonomy.get_node_by_id(504556).unwrap();
        assert_eq!(fixture.taxonomy.get_id_by_node(*node_id).unwrap(), 504556);
    }

    #[test]
    fn get_name_by_id() {
        let fixture = NcbiTaxonomyFixture::default();
        assert_eq!(fixture.taxonomy.get_name_by_id(370556).unwrap(), "Streptococcus phage 9429.1");
    }

    #[test]
    fn traversal() {
        let fixture = NcbiTaxonomyFixture::default();
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
        let fixture = NcbiTaxonomyFixture::default();
        assert!(fixture.taxonomy.is_descendant("Propionibacterium phage PAS7", "unclassified bacterial viruses"));
    }

    #[test]
    fn taxid_descendants() {
        let fixture = NcbiTaxonomyFixture::default();
        assert!(fixture.taxonomy.is_descendant_taxid(504556, 12333));
    }
}
