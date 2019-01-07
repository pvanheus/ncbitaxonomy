[![CircleCI](https://circleci.com/gh/pvanheus/ncbitaxonomy.svg?style=svg)](https://circleci.com/gh/pvanheus/ncbitaxonomy)

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
$ taxonomy_filter_refseq --help
taxonomy_filter_refseq 0.1.2
Peter van Heusden <pvh@sanbi.axc.za>
Filter NCBI RefSeq FASTA files by taxonomic lineage

USAGE:
    taxonomy_filter_refseq [OPTIONS] <INPUT_FASTA> <TAXONOMY_DIR> <ANCESTOR_NAME> [OUTPUT_FASTA]

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -t, --tax_prefix <TAXONOMY_FILENAME_PREFIX>    String to prepend to names of nodes.dmp and names.dmp

ARGS:
    <INPUT_FASTA>      FASTA file with RefSeq sequences
    <TAXONOMY_DIR>     Directory containing the NCBI taxonomy nodes.dmp and names.dmp files
    <ANCESTOR_NAME>    Name of ancestor to use as ancestor filter
    <OUTPUT_FASTA>     Output FASTA filename (or stdout if omitted)

```

### TODO

* Clean up non-idiomatic code (e.g. the use of the insert_new_entry bool)
* Refactor taxonomy_filter_refseq: move most code to library, add tests