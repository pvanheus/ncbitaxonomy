### ncbitaxonomy

This is a Rust crate (i.e. library) for working with a local copy of the 
[NCBI Taxonomy database](https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/).
The database can be downloaded (either `taxdump.zip` or `taxdump.tar.gz`) from the
[NCBI Taxonomy FTP site](https://ftp.ncbi.nih.gov/pub/taxonomy/).

Documentation for version 0.1.0 is available at [crates.io](https://docs.rs/ncbitaxonomy/0.1.0/ncbitaxonomy/struct.NcbiTaxonomy.html).

### taxonomy_filter_refseq

(new in 0.1.1)

A tool to filter a NCBI RefSeq FASTA file so that only the ancestors of a given taxon
are retained.

```bash
$ taxonomy_filter_refseq
error: The following required arguments were not provided:
    <INPUT_FASTA>
    <TAXONOMY_DIR>
    <ANCESTOR_NAME>

USAGE:
    taxonomy_filter_refseq [OPTIONS] <INPUT_FASTA> <TAXONOMY_DIR> <ANCESTOR_NAME> [OUTPUT_FASTA]

```

### TODO

* Clean up non-idiomatic code (e.g. the use of the insert_new_entry bool)
* Add testing via CI
* Refactor taxonomy_filter_refseq: move most code to library, add tests