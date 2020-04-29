[![CircleCI](https://circleci.com/gh/pvanheus/ncbitaxonomy.svg?style=svg)](https://circleci.com/gh/pvanheus/ncbitaxonomy)

### ncbitaxonomy

This is a Rust crate (i.e. library) for working with a local copy of the 
[NCBI Taxonomy database](https://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/).
The database can be downloaded (either `taxdump.zip` or `taxdump.tar.gz`) from the
[NCBI Taxonomy FTP site](https://ftp.ncbi.nih.gov/pub/taxonomy/).

Documentation for version 0.1.3 is available at [crates.io](https://docs.rs/ncbitaxonomy/0.1.3/ncbitaxonomy/struct.NcbiTaxonomy.html).

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
        --no_curated      Don't accept curated RNAs and proteins (NM_, NR_ and NP_ accessions)
        --no_predicted    Don't accept computationally predicted RNAs and proteins (XM_, XR_ and XP_ accessions)
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

### taxonomy_filter_fastq

(new in version 0.2.0)

```bash
$ taxonomy_filter_fastq --help
taxonomy_filter_refseq 0.1.2
Peter van Heusden <pvh@sanbi.axc.za>
Filter NCBI RefSeq FASTA files by taxonomic lineage

USAGE:
    taxonomy_filter_fastq [FLAGS] [OPTIONS] <INPUT_FASTQ> --ancestor_taxid <ANCESTOR_ID> --taxdir <TAXONOMY_DIR> --tax_report_filename <TAXONOMY_REPORT_FILENAME> <--centrifuge|--kraken2>

FLAGS:
    -d, --output_dir    Directory to deposited filtered output files in
    -C, --centrifuge    Filter using report from Centrifuge
    -h, --help          Prints help information
    -K, --kraken2       Filter using report from Kraken2
    -V, --version       Prints version information

OPTIONS:
    -A, --ancestor_taxid <ANCESTOR_ID>                      Name of ancestor to use as ancestor filter
    -T, --taxdir <TAXONOMY_DIR>
            Directory containing the NCBI taxonomy nodes.dmp and names.dmp files

    -t, --tax_prefix <TAXONOMY_FILENAME_PREFIX>             String to prepend to names of nodes.dmp and names.dmp
    -F, --tax_report_filename <TAXONOMY_REPORT_FILENAME>    Output from Kraken2 (default) or Centrifuge

ARGS:
    <INPUT_FASTQ>    FASTA file with RefSeq sequences
```
