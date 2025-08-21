HumGut2 - a human gut genome collection
================
Lars Snipen

- [Introduction](#introduction)
- [Download HumGut](#download-humgut)
  - [FASTA files](#fasta-files)
  - [Metadata tables](#metadata-tables)
- [Taxonomy](#taxonomy)
  - [Taxonomy tree](#taxonomy-tree)
- [Building a kraken2 database](#building-a-kraken2-database)
  - [The taxonomy](#the-taxonomy)
  - [Other genomes](#other-genomes)
  - [The HumGut library](#the-humgut-library)
  - [Building](#building)
- [Building a bracken database](#building-a-bracken-database)

<style>
.r {background-color:  #CCCCCC; }
.r {color: #000044;}
</style>

# Introduction

The human gut metagenome is the focus of a lot of research in our time.
From this site you can download the *HumGut genome collection* and
accompanying metadata. As described in our paper, genomes encountered in
healthy human guts worldwide were ranked by prevalence and clustered by
whole genome identity (97.5%). Genomes representing the clusters were
retained as the HumGut genome collection.

The HumGut2 is simply a small extension of the original genome
collection, including many more genomes in the scan for relevant
members. The HumGut2 has 31,225 genomes. This is only around 1000
genomes more than the original, indicating that even if we increase the
number of genomes to search through we do not expand the collection of
relevant genomes very much anymore at this resolution.

If you use this resource, we would ask you to cite our paper:

[HumGut: A comprehensive human gut prokaryotic genomes collection
filtered by metagenome data](https://pubmed.ncbi.nlm.nih.gov/34330336/)

<br> <br>

# Download HumGut

<br>

## FASTA files

The HumGut2 collection contains 31,225 genomes, which is a carefully
extracted representative subset from a much larger pool of more than
400,000 genomes to which reads from thousands of human gut metagenomes
map. In this compressed archive:

- **[HumGut2.tar](https://arken.nmbu.no/~larssn/humgut/HumGut2.tar)**

you find the compressed FASTA files, one for each genome. All files have
Header-lines equipped with the proper text for building a kraken2
database, see sections below for more details on this. This compressed
archive is roughly 18GB. <br>

## Metadata tables

The tab-separated text file

- **[HumGut2.tsv](https://arken.nmbu.no/~larssn/humgut/HumGut2.tsv)**

lists metadata about each HumGut genome (31,225 rows). Below you find a
description of all columns.

| Column | Description |
|----|----|
| `HumGut_name` | Unique HumGut name for each genome. |
| `cluster975` | The highest resolution cluster this genome belongs to (97.5% sequence identity). |
| `cluster95` | The coarser resolution cluster this genome belongs to (95% sequence identity). |
| `gtdbtk_tax_id` | GTDB-tk genome taxonomy ID. Note that the GTDB database (<https://gtdb.ecogenomic.org/>) has no such integer identifiers, and we have just artificially created some here. This is required for building kraken2/bracken/krakenUniq databases using this taxonomy. These integers are from 4 000 000 and up, a choice made to not interfere with the NCBI Taxonomy. |
| `gtdbtk_name` | Genome name as given by GTDBTk. |
| `gtdbtk_taxonomy` | The full GTDBTk taxonomy, from domain and down. |
| `prevalence_score` | The average sequence identity across 3,534 healthy human gut metagenomes. |
| `metagenomes_present` | The number of metagenomes where the genome was found present, using ≥ 95% sequence identity as a threshold. |
| `completeness` | The estimated completeness (%) of the genome. |
| `contamination` | The estimated contamination (%) of the genome. |
| `GC` | Genome GC content. |
| `genome_size` | Number of base pairs in genome. |
| `source` | Either RefSeq (<https://ftp.ncbi.nlm.nih.gov/genomes/refseq/>) or UHGG (<https://www.ebi.ac.uk/metagenomics/>) |
| `genome_type` | The completion level as listed in RefSeq, or MAG (all UHGG genomes). |
| `cluster975_size` | Number of genomes in the same cluster of the highest resolution (97.5%). |
| `cluster95_size` | Number of genomes in the same coarse cluster (95%). |
| `genome_file` | The name of the FASTA file in the archive **HumGut2.tar.gz** from above. |

<br> <br>

# Taxonomy

One obvious use of the HumGut collection is to assign some taxonomy to
the reads you have after sequencing a human gut metagenome. We have
assigned all HumGut genomes to the
<a href="https://gtdb.ecogenomic.org/" target="blank">GTDB database</a>
taxonomy (version 2.26) using the
<a href="https://github.com/Ecogenomics/GTDBTk" target="blank">GTDBTk
software</a>.

In the table **HumGut2.tsv** mentioned above, each HumGut genome has a
taxonomy identifier (`gtdb_tax_id`) and taxonomy (`gtdb_taxonomy`). The
GTDB database does not assign numeric identifiers to their taxa, and the
`gtdb_tax_id` in the HumGut metadata are just integer identifiers we
have ‘created’. They are all at least 4 million in value in order not to
overlap with the corresponding NCBI taxonomy identifiers.

Some HumGut genomes may lack taxonomy if they are simply too different
from any taxon listed by GTDB. <br>

## Taxonomy tree

In order to describe the taxonomy tree, we have chosen to use the data
structures used by NCBI Taxonomy
(<https://ftp.ncbi.nih.gov/pub/taxonomy/>). This is also what tools like
kraken2 uses (see below).

Here you can download the two files needed for using the GTDB taxonomy:

- **[GTDB_names.dmp](https://arken.nmbu.no/~larssn/humgut/GTDB_names.dmp)**
- **[GTDB_nodes.dmp](https://arken.nmbu.no/~larssn/humgut/GTDB_nodes.dmp)**

These files contains the GTDB taxonomy names and the `gtdb_tax_id`
mentioned above for all taxa in the clades below domain `Archaea` and
`Bacteria`. For all other clades in the tree of life these files use the
NCBI Taxonomy. Thus, you may use these files when creating a kraken2
database where you combine our HumGut genomes with other genomes like
human, fungi or viruses. The latter will then be classified by the NCBI
taxonomy while all prokaryotes have the GTDB taxonomy, <br> <br>

# Building a kraken2 database

The kraken2 software (<https://github.com/DerrickWood/kraken2>) is a
popular tool for making taxonomic classification of metagenome reads.
Building a kraken2 database from the HumGut genomes gives you an
excellent tool for taxonomic profiling of data from the human gut. Below
we describe a procedure for building such a kraken2 database.

Make a folder in which you want to build the kraken2 database. We refer
to this as `$KRAKEN2` here. <br>

## The taxonomy

Download the files `GTDB_names.dmp` and `GTDB_nodes.dmp` mentioned
above.

Make the subfolder `$KRAKEN2/taxonomy`. Copy `GTDB_names.dmp` and
`GTDB_nodes.dmp` into this:

``` bash
cp GTDB_names.dmp $KRAKEN2/taxonomy/names.dmp
cp GTDB_nodes.dmp $KRAKEN2/taxonomy/nodes.dmp
```

Note that the names of the copied files must be `names.dmp` and
`nodes.dmp` for kraken2 to recognize them. <br>

## Other genomes

We strongly recommend you also include the human genome in the database
as long as your data are from the human gut. Here is the code for
including the human genome:

``` bash
kraken2-build --download-library human --db $KRAKEN2
```

The folder `$KRAKEN2/library` should appear, and inside it, the `human/`
subfolder with some data. There should at least be a file named
`library.fna`, around 3.1GB, with the human genome sequences. Note that
this download may fail, and you may need to repeat this step until it
has downloaded completely. It usually takes a few minutes.

You may also need to use the `--use-ftp` option if the default RSYNC way
of downloading is unavailable.

You may also include other kraken2 libraries in the same way,
e.g. fungi, virus etc. <br>

## The HumGut library

We assume you have extracted the `HumGut.tar.gz` from above into the
subfolder `fna`. If you include all HumGut genomes, simply write all
FASTA-files from this folder into a single uncompressed FASTA file. The
latter because kraken2 cannot build from compressed FASTA files. It can
be done like this

``` bash
zcat fna/*.fna.gz > HumGut975_library.fna
```

The file `HumGut975_library.fna` should be close to 60GB.

Then you add this to the kraken2 database with

``` bash
kraken2-build --add-to-library HumGut975_library.fna --db $KRAKEN2
```

A subfolder `$KRAKEN2/library/added` should now appear. Inside it a
fasta-file should appear, containing the library sequences, i.e. should
be of the same size as the `HumGut975_library.fna` or whatever file you
used.

This step will take some time, depending on the size of the library. You
could speed this up by using multiple threads. After this step, you may
delete the huge FASTA file you created above (`HumGut975_library.fna`).

In the example above, we included all HumGut genomes in the database. If
this requires too much memory, or you simply just want a lower
resolution, you may only use the 95% clustered genomes in the database
instead of the full resolution. Here is some R code for creating a
corresponding library file:

``` r
library(tidyverse)
humgut950.tbl <- read_delim("HumGut2.tsv", delim = "\t") %>% 
  distinct(cluster95, .keep_all = T)
ok <- file.append("HumGut95_library.fna.gz",
                  file.path("fna", humgut950.tbl$genome_file))
```

It basically filters the table in `HumGut2.tsv` to keep only the rows
where you find the *first* occurrence of the unique values in the
`cluster95` column. These are the genomes representing the 95% clusters.
In `HumGut2.tsv` the rows are sorted in descending order by the
`prevalence_score` and hence, the first row for each cluster is the
genome to keep. Note that in the above code we assume the file
`HumGut2.tsv` and the folder `fna` (with all fasta files) is in the
current working directory, please use correct paths if they are
elsewhere. This code produces a compressed FASTA file, and you need to
uncompress it before you call `kraken2-build`. This file should be
around 11GB when uncompressed.

You may of course also select all kinds of other subsets of the HumGut
genomes to include in your database, using a similar approach. <br>

## Building

The last step is just to build the database

``` bash
kraken2-build --build --threads 20 --db $KRAKEN2
```

Here we used 20 threads. This step is the most time-consuming.

The files hash.k2d, opts.k2d, taxo.k2d and seqid2taxid.map should appear
in the `$KRAKEN2` folder when the build is complete.

NOTE: You may run `kraken2-build --clean` to delete files no longer
needed in the `$KRAKEN` folder, but *do not* do this yet if you intend
to also build a bracken database (see below), as some of these files
will still be needed. <br> <br>

# Building a bracken database

Once you have the kraken2 database, it is straightforward to extend this
by building a bracken (<https://github.com/jenniferlu717/Bracken>)
database as well:

``` bash
bracken-build -d $KRAKEN2 -t 10 -k 35 -l 100
```

Here we used `10` threads, k-mers of length `35` (same as kraken2) and
read length `100`. This should add three files to the `$KRAKEN2` folder:
database.kraken, database.100mers.kraken, database.100mers.kmer_distrib.
See the bracken software GitHub site for more details:
<https://github.com/jenniferlu717/Bracken>.

<!-- # Building a krakenUniq database -->

<!-- As an alternative to the kraken2/bracken approach, you may also be interested in using the software krakenUniq (https://github.com/fbreitwieser/krakenuniq). -->

<!-- Once you have the kraken2 database, it is straightforward to build a krakenUniq database. Note that this requires much more memory than the kraken2/bracken. For the full HumGut collection I allocate 500GB of memory for this job, and you should have at least 600GB of free disk space. -->

<!-- First, make a folder for you krakenUniq database, we denote this `$KRAKENUNIQ`. -->

<!-- Next, copy both the `taxonomy/` and the `library/` folders from the kraken2 database to the new folder: -->

<!-- ```{bash, eval=FALSE} -->

<!-- cp -r $KRAKEN2/taxonomy $KRAKENUNIQ -->

<!-- cp -r $KRAKEN2/library $KRAKENUNIQ -->

<!-- ``` -->

<!-- Finally, build the krakenUniq database: -->

<!-- ```{.bash} -->

<!-- krakenuniq-build --threads 20 --db $KRAKENUNIQ -->

<!-- ``` -->

<!-- You may free a lot of disk space by running -->

<!-- ```{bash, eval=FALSE} -->

<!-- krakenUniq-build --clean -db $KRAKENUNIQ -->

<!-- ``` -->

<!-- In my case the full HumGut database was reduced to occupy around 350GB disk space after this step. A similar cleaning may also be done for the `KRAKEN2` database now. -->
