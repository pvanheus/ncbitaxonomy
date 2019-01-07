/// currently just testing code - please ignore
pub fn main() {
    let tree = ncbitaxonomy::NcbiTaxonomy::from_ncbi_files("/home/pvh/Documents/code/Masters/taxonomy_filter/ncbi_data/nodes.dmp", "/home/pvh/Documents/code/Masters/taxonomy_filter/ncbi_data/names.dmp").unwrap();
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
                ncbitaxonomy::NodeEdge::Start(node_id) => {
                    let id = tree.get_id_by_node(node_id).unwrap();
                    println!("start: {}", id)
                },
                ncbitaxonomy::NodeEdge::End(node_id) => {
                    let id = tree.get_id_by_node(node_id).unwrap();
                    println!("end: {}", id)
                }
            }
        }
    }
}