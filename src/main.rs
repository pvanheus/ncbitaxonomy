extern crate tree;
use tree::ncbitaxonomy;

pub fn main() {
    let tree = ncbitaxonomy::NcbiTaxonomy::read_ncbi_files("/home/pvh/Documents/code/Masters/taxonomy_filter/ncbi_data/nodes.dmp", "/home/pvh/Documents/code/Masters/taxonomy_filter/ncbi_data/names.dmp").unwrap();
    println!("finished reading tree");
    let traversal = tree.traversal(1);
    if let Some(traversal) = traversal {
        let mut counter = 0;
        for node_id in traversal {
            if counter > 1000 {
                break;
            }
            counter += 1;
            match node_id {
                tree::ncbitaxonomy::NodeEdge::Start(node_id) => {
                    let id = tree.get_id_by_node(node_id).unwrap();
                    println!("start: {}", id)
                },
                tree::ncbitaxonomy::NodeEdge::End(node_id) => {
                    let id = tree.get_id_by_node(node_id).unwrap();
                    println!("end: {}", id)
                }
            }
        }
    }
}