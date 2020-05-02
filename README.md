[![CircleCI](https://circleci.com/gh/pvanheus/ncbitaxonomy.svg?style=svg)](https://circleci.com/gh/pvanheus/ncbitaxonomy)

### ncbitaxonomy

This is a Rust crate (i.e. library) for working with a local copy of the 
[NCBI Taxonomy database](https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/).
The database can be downloaded (either `taxdump.zip` or `taxdump.tar.gz`) from the
[NCBI Taxonomy FTP site](https://ftp.ncbi.nih.gov/pub/taxonomy/) and reformatted into a SQLite database
using the `taxonomy_util` utility's `to_sqlite` subcommand.

Documentation is available at [crates.io](https://crates.io/crates/ncbitaxonomy).

### taxonomy_filter_refseq

(new in 0.1.1)

A tool to filter a NCBI RefSeq FASTA file so that only the ancestors of a given taxon
are retained.

```bash
$ taxonomy_filter_refseq --help
taxonomy_filter_refseq 1.0.0
Peter van Heusden <pvh@sanbi.axc.za>
Filter NCBI RefSeq FASTA files by taxonomic lineage

USAGE:
    taxonomy_filter_refseq [FLAGS] [OPTIONS] <INPUT_FASTA> <ANCESTOR_NAME> [OUTPUT_FASTA]

FLAGS:
        --no_curated      Don't accept curated RNAs and proteins (NM_, NR_ and NP_ accessions)
        --no_predicted    Don't accept computationally predicted RNAs and proteins (XM_, XR_ and XP_ accessions)
    -h, --help            Prints help information
    -V, --version         Prints version information

OPTIONS:
    -d, --db <TAXDB_URL>    URL for SQLite taxonomy database

ARGS:
    <INPUT_FASTA>      FASTA file with RefSeq sequences
    <ANCESTOR_NAME>    Name of ancestor to use as ancestor filter
    <OUTPUT_FASTA>     Output FASTA filename (or stdout if omitted)
```

### taxonomy_filter_fastq

(new in version 0.2.0)


```bash
$ taxonomy_filter_fastq --help
taxonomy_filter_fastq 1.0.0
Peter van Heusden <pvh@sanbi.axc.za>
Filter FASTQ files whose reads have been classified by Centrifuge or Kraken2, only retaining reads in taxa descending
from given ancestor

USAGE:
    taxonomy_filter_fastq [FLAGS] [OPTIONS] <INPUT_FASTQ>... --ancestor_taxid <ANCESTOR_ID> --tax_report_filename <TAXONOMY_REPORT_FILENAME> <--centrifuge|--kraken2>

FLAGS:
    -d, --output_dir    Directory to deposited filtered output files in
    -C, --centrifuge    Filter using report from Centrifuge
    -h, --help          Prints help information
    -K, --kraken2       Filter using report from Kraken2
    -V, --version       Prints version information

OPTIONS:
    -A, --ancestor_taxid <ANCESTOR_ID>                      Name of ancestor to use as ancestor filter
    -d, --db <TAXDB_URL>                                    URL for SQLite taxonomy database
    -F, --tax_report_filename <TAXONOMY_REPORT_FILENAME>    Output from Kraken2 (default) or Centrifuge

ARGS:
    <INPUT_FASTQ>...    FASTA file with RefSeq sequences
```

### taxonomy_util

(new in 1.0.0)

Utilities to convert NCBI taxonomy database files into SQLite database (the input format used in other tools).

```bash
taxonomy_util 1.0.0
Peter van Heusden <pvh@sanbi.axc.za>
Utilities for working with the NCBI taxonomy database

USAGE:
    taxonomy_util [OPTIONS] [SUBCOMMAND]

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -d, --db <TAXDB_URL>    URL for SQLite taxonomy database

SUBCOMMANDS:
    common_ancestor_distance    find the tree distance to te common ancestor between two taxa
    get_id                      find taxonomy ID for name
    get_lineage                 get lineage for name [unimplemented]
    get_name                    find name for taxonomy ID
    help                        Prints this message or the help of the given subcommand(s)
    to_sqlite                   save taxonomy database loaded from files to SQLite database file
```