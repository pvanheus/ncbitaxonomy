#![recursion_limit = "1024"]
#[macro_use]
extern crate error_chain;

extern crate indextree;

pub mod ncbitaxonomy {

    mod errors {
        error_chain! {
            foreign_links {
                Io(::std::io::Error) #[doc = "Io"];
                ParseIntError(::std::num::ParseIntError);
            }
            errors {
                NodeFileFormatError(l: String) {
                    description("format error in nodes.dmp file")
                    display("format error in nodes.dmp in line {}", l)
                }
            }
        }
    }

    use self::errors::*;

    use std::collections::HashMap;
    use std::fs::File;
    use std::io::{BufReader,BufRead};
    use indextree::{Arena, NodeId, Traverse};
    pub use indextree::NodeEdge;

    #[derive(Debug)]
    pub struct NcbiTaxonomy {
        arena: Arena<u32>,
        name_to_node: HashMap<String, NodeId>,
        id_to_node: HashMap<u32, NodeId>,
        id_to_name: HashMap<u32, String>,
    }

    impl NcbiTaxonomy {

        pub fn read_ncbi_files(nodes_filename: &str, names_filename: &str) -> Result<NcbiTaxonomy> {
            let mut child_ids_by_parent_id: HashMap<u32, Vec<u32>> = HashMap::new();
            let nodes_file = File::open(nodes_filename)?;
            for line_maybe in BufReader::new(nodes_file).lines() {
                match line_maybe {
                    Ok(line) => {
                        let mut fields = line.split("\t|\t");
                        let id_str = fields.next().chain_err(|| ErrorKind::NodeFileFormatError(line.clone()))?;
                        let parent_id_str = fields.next().chain_err(|| ErrorKind::NodeFileFormatError(line.clone()))?;
                        let id = id_str.parse::<u32>().chain_err(|| format!("failed to parse id_str as u32: {}", id_str))?;
                        let parent_id = parent_id_str.parse::<u32>().chain_err(|| format!("failed to parse parent_id_str as u32: {}", parent_id_str))?;
                        if parent_id != id {  // this happens for the root node
                            let mut insert_new_entry = false;
                            match child_ids_by_parent_id.get_mut(&parent_id) {
                                Some(entry) => { entry.push(id)},
                                None => { insert_new_entry = true },
                            }
                            if insert_new_entry {
                                child_ids_by_parent_id.insert(parent_id, vec![id]);
                            }
                        }
                    },
                    Err(e) => return Err(ErrorKind::Io(e).into())
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
                let line = match line_maybe {
                    Ok(line) => line,
                    Err(e) => return Err(ErrorKind::Io(e).into())
                };
                let mut fields = line.split("\t|\t").collect::<Vec<&str>>();
                if fields[3].starts_with("scientific name") {
                    let id_str = fields[0];
                    let id = id_str.parse::<u32>().chain_err(|| format!("failed to parse id_str as u32: {}", id_str))?;
                    let name = fields[1].to_string();
                    let node_id = id_to_node.get(&id).expect("ID not found in id_to_node");
                    id_to_name.insert(id, name.clone());
                    name_to_node.insert(name, *node_id);
                }
            }

            let tree = NcbiTaxonomy { arena, name_to_node, id_to_node, id_to_name };
            return Ok(tree);
        }

        pub fn contains_id(&self, id: u32) -> bool {
            self.id_to_node.contains_key(&id)
        }

        pub fn contains_name(&self, name: &str) -> bool {
            self.name_to_node.contains_key(name)
        }

        pub fn is_descendant(&self, name: &str, ancestor_name: &str) -> bool {
            let id = match self.name_to_node.get(name) {
                Some(id) => id,
                None => return false
            };
            let ancestor_id = match self.name_to_node.get(ancestor_name) {
                Some(id) => id,
                None => return false
            };
            for node in id.ancestors(&self.arena).into_iter() {
                if node == *ancestor_id {
                    return true
                }
            }
            return false
        }

        pub fn get_node_by_id(&self, id: u32) -> Option<&NodeId> {
            self.id_to_node.get(&id)
        }

        pub fn traversal(&self, from: u32) -> Option<Traverse<u32>> {
            match self.get_node_by_id(from) {
                Some(node_id) => Some(node_id.traverse(&self.arena)),
                None => None
            }
        }

        pub fn get_id_by_node(&self, node_id: NodeId) -> Option<u32> {
            match self.arena.get(node_id)  {
                Some(node) => Some(node.data),
                None => None
            }
        }

        pub fn get_name_by_id(&self, id: u32) -> Option<&String> {
            self.id_to_name.get(&id)
        }
    }

}


#[cfg(test)]
mod tests {
    use ncbitaxonomy::{NcbiTaxonomy, NodeEdge};

    pub struct NcbiTreeFixture {
        pub tree: NcbiTaxonomy,
    }

    impl Default for NcbiTreeFixture {
        fn default() -> Self {
            let tree = NcbiTaxonomy::read_ncbi_files("data/sample_tree_nodes.dmp", "data/sample_tree_names.dmp").unwrap();
            Self { tree }
        }
    }

    #[test]
    fn contains_id() {
        let fixture = NcbiTreeFixture::default();
        assert!(fixture.tree.contains_id(504556));
    }

    #[test]
    fn contains_name() {
        let fixture = NcbiTreeFixture::default();
        assert!(fixture.tree.contains_name("Propionibacterium phage PAS7"));
    }


    #[test]
    fn get_node_by_id() {
        let fixture = NcbiTreeFixture::default();
        assert_eq!(fixture.tree.get_node_by_id(999999999), None);
        assert!(match fixture.tree.get_node_by_id(504556) { Some(_) => true, None => false })
    }

    #[test]
    fn get_id_by_node() {
        let fixture = NcbiTreeFixture::default();
        let node_id = fixture.tree.get_node_by_id(504556).unwrap();
        assert_eq!(fixture.tree.get_id_by_node(*node_id).unwrap(), 504556);
    }

    #[test]
    fn get_name_by_id() {
        let fixture = NcbiTreeFixture::default();
        assert_eq!(fixture.tree.get_name_by_id(370556).unwrap(), "Streptococcus phage 9429.1");
    }

    #[test]
    fn traversal() {
        let fixture = NcbiTreeFixture::default();
        let traversal = fixture.tree.traversal(12333);
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
        let fixture = NcbiTreeFixture::default();
        assert!(fixture.tree.is_descendant("Propionibacterium phage PAS7", "unclassified bacterial viruses"));
    }
}
