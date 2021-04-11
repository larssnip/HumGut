HumGut - a human gut genome collection
================
Lars Snipen

-   [Introduction](#introduction)
-   [Download HumGut](#download-humgut)
    -   [FASTA files](#fasta-files)
    -   [Metadata tables](#metadata-tables)
-   [Taxonomy](#taxonomy)
    -   [Taxonomy tree](#taxonomy-tree)
-   [Building a kraken2 database](#building-a-kraken2-database)
    -   [The taxonomy](#the-taxonomy)
    -   [The human genome](#the-human-genome)
    -   [The HumGut library](#the-humgut-library)
    -   [Building](#building)
-   [Building a bracken database](#building-a-bracken-database)
-   [Building a krakenUniq database](#building-a-krakenuniq-database)

# Introduction

The human gut metagenome is the focus of a lot of research in our time.
From this site you can download the *HumGut genome collection* and
accompanying metadata. As described in our paper, genomes encountered in
healthy human guts worldwide were ranked by prevalence and clustered by
whole genome identity (97.5%). Genomes representing the clusters, 30 691
in total, were retained as HumGut.

If you use this resource, we would ask you to cite our paper:

[HumGut: A comprehensive human gut prokaryotic genomes collection
filtered by metagenome
data](https://www.biorxiv.org/content/10.1101/2020.03.25.007666v2.full)

as well as the underlying data repositories at EMBL-EBI
(<https://www.ebi.ac.uk/metagenomics/>) and NCBI
(<https://www.ncbi.nlm.nih.gov/genome>). <br> <br>

# Download HumGut

<br>

## FASTA files

The HumGut collection contains 30,691 genomes. In this compressed
archive:

-   **[HumGut.tar.gz](http://arken.nmbu.no/~larssn/humgut/HumGut.tar.gz)**

you find the compressed FASTA files, one for each genome. All files have
Header-lines equipped with the proper text for building a kraken2 or
krakenUniq database, see sections below for more details on this. This
archive is roughly 18GB. <br>

## Metadata tables

The tab-separated text file

-   **[HumGut.tsv](http://arken.nmbu.no/~larssn/humgut/HumGut.tsv)**

lists metadata about each HumGut genome (30,691 rows). Below you find a
description of all columns.

| Column                 | Description                                                                                                                                                                                                                                                                                                                                                                                                |
|------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `HumGut_name`          | Unique HumGut name for each genome.                                                                                                                                                                                                                                                                                                                                                                        |
| `HumGut_tax_id`        | Unique HumGut identifier for each genome. These are integers from 3 000 000 and up, a choice made to not interfere with the NCBI Taxonomy database integers.                                                                                                                                                                                                                                               |
| `cluster975`           | The highest resolution cluster this genome belongs to (97.5% sequence identity).                                                                                                                                                                                                                                                                                                                           |
| `cluster95`            | The coarser resolution cluster this genome belongs to (95% sequence identity).                                                                                                                                                                                                                                                                                                                             |
| `gtdbtk_tax_id`        | GTDB-tk genome taxonomy ID. Note that the GTDB database (<https://gtdb.ecogenomic.org/>) has no such integer identifiers, and we have just artificially created some here. This is required for building kraken2/bracken/krakenUniq databases using this taxonomy. These integers are from 4 000 000 and up, a choice made to not interfere with the NCBI Taxonomy or the HumGut\_tax\_id mentioned above. |
| `gtdbtk_organism_name` | Genome name as given by GTDB-tk.                                                                                                                                                                                                                                                                                                                                                                           |
| `gtdbtk_taxonomy`      | The full GTDB-tk taxonomy, from domain and down.                                                                                                                                                                                                                                                                                                                                                           |
| `ncbi_tax_id`          | The taxonomy identifier from the NCBI Taxonomy database (<https://www.ncbi.nlm.nih.gov/taxonomy/>).                                                                                                                                                                                                                                                                                                        |
| `ncbi_organism_name`   | Genome name at the NCBI Taxonomy database.                                                                                                                                                                                                                                                                                                                                                                 |
| `ncbi_rank`            | The rank at the NCBI Taxonomy database.                                                                                                                                                                                                                                                                                                                                                                    |
| `prevalence_score`     | The average sequence identity across 3,534 healthy human gut metagenomes.                                                                                                                                                                                                                                                                                                                                  |
| `metagenomes_present`  | The number of metagenomes where the genome was found present, using ≥ 95% sequence identity as a threshold.                                                                                                                                                                                                                                                                                                |
| `completeness`         | The estimated completeness (%) of the genome.                                                                                                                                                                                                                                                                                                                                                              |
| `contamination`        | The estimated contamination (%) of the genome.                                                                                                                                                                                                                                                                                                                                                             |
| `GC`                   | Genome GC content.                                                                                                                                                                                                                                                                                                                                                                                         |
| `genome_size`          | Number of base pairs in genome.                                                                                                                                                                                                                                                                                                                                                                            |
| `source`               | Either RefSeq (<https://ftp.ncbi.nlm.nih.gov/genomes/refseq/>) or UHGG (<https://www.ebi.ac.uk/metagenomics/>)                                                                                                                                                                                                                                                                                             |
| `genome_type`          | The completion level as listed in RefSeq, or MAG (all UHGG genomes).                                                                                                                                                                                                                                                                                                                                       |
| `cluster975_size`      | Number of genomes in the same cluster of the highest resolution (97.5%).                                                                                                                                                                                                                                                                                                                                   |
| `cluster95_size`       | Number of genomes in the same coarse cluster (95%).                                                                                                                                                                                                                                                                                                                                                        |
| `genome_file`          | The name of the FASTA file in the archive **HumGut.tar.gz** from above.                                                                                                                                                                                                                                                                                                                                    |
| `ftp_download`         | The ftp address from where we downloaded the genome.                                                                                                                                                                                                                                                                                                                                                       |

The tab-separated text file

-   **[All\_genomes.tsv](http://arken.nmbu.no/~larssn/humgut/All_genomes.tsv)**

lists metadata about all the 381,779 genomes used for obtaining the
HumGut clusters (381,779 rows). The columns are a subset of those in the
above table, see the above description of column names. Note that we do
not provide the FASTA files for all these genomes at this website, since
they are publicly available elsewhere. The FTP addresses in the column
`ftp_download` shows where each genome is found. <br> <br>

# Taxonomy

One obvious use of the HumGut collection is to assign some taxonomy to
the reads you have after sequencing a human gut metagenome. In the table
**HumGut.tsv** mentioned above, each HumGut genomes has a taxonomy
identifier (`HumGut_tax_id`). You also find the columns `gtdbtk_tax_id`
and `ncbi_tax_id` in the same table. These are the *parents* of the
HumGut genome in the GTDB or the NCBI taxonomy, respectively. Be aware
that even if each HumGut genome is clustered at a sub-species identity
threshold, its parent may not be a species, but sometimes a genus or
even higher rank. This may be because some branches in the taxonomy do
not contain all ranks, or the HumGut genome is simply too different from
any species listed by GTDB or NCBI, and has been assigned directly under
some higher rank. <br>

## Taxonomy tree

In order to describe the taxonomy tree, we have chosen to use the data
structures used by NCBI Taxonomy
(<https://ftp.ncbi.nih.gov/pub/taxonomy/>). This is also what tools like
kraken2 and krakenUniq uses (see below).

Here you can download the two files needed for using the GTDB taxonomy:

-   **[gtdb\_names.dmp](http://arken.nmbu.no/~larssn/humgut/gtdb_names.dmp)**
-   **[gtdb\_nodes.dmp](http://arken.nmbu.no/~larssn/humgut/gtdb_nodes.dmp)**

Here you can download the files needed for using the NCBI taxonomy:

-   **[ncbi\_names.dmp](http://arken.nmbu.no/~larssn/humgut/ncbi_names.dmp)**
-   **[ncbi\_nodes.dmp](http://arken.nmbu.no/~larssn/humgut/ncbi_nodes.dmp)**

Both these sets of files contain the `HumGut_tax_id` for all HumGut
genomes. Their parents are the corresponding `gtdbtk_tax_id` or
`ncbi_tax_id`. The files then contain the branches leading down to these
taxa. Note that the files *do not* contain the full taxonomy for all
taxa in either database. They have been pruned to only contain the
relevant branches leading to some HumGut genome.

In addition, both pairs of files include the human genome branch, with
its NCBI taxonomy. This has been included since we believe the human
genome should always be included in a reference database for reads from
the human gut (as possible contaminations).

Both taxonomies above are from January 2021. This changes slowly over
time. We will make efforts to update this at regular intervals. <br>
<br>

# Building a kraken2 database

The kraken2 software (<https://github.com/DerrickWood/kraken2>) is a
popular tool for making taxonomic classification of metagenome reads.
Building a kraken2 database from the HumGut genomes gives you an
excellent tool for taxonomic profiling of data from the human gut. Below
we describe a procedure for building such a kraken2 database.

Make a folder in which you want to build the kraken2 database. We refer
to this as `$KRAKEN2` here. <br>

## The taxonomy

Download the files `gtdb_names.dmp` and `gtdb_nodes.dmp` mentioned
above, or the corresponding `ncbi_` files if you want to use the NCBI
taxonomy. The procedure is exactly the same in both cases, and we use
the `gtdb_`files in the code examples.

Make the subfolder `$KRAKEN2/taxonomy`. Copy `gtdb_names.dmp` and
`gtdb_nodes.dmp` into this:

``` bash
cp gtdb_names.dmp $KRAKEN2/taxonomy/names.dmp
cp gtdb_nodes.dmp $KRAKEN2/taxonomy/nodes.dmp
```

Note that the names of the copied files must be `names.dmp` and
`nodes.dmp` for kraken2 to recognize them. <br>

## The human genome

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
of downloading is unavailable. <br>

## The HumGut library

We assume you have extracted the `HumGut.tar.gz` from above into the
folder `$FNA`. If you include all HumGut genomes, simply write all
FASTA-files from this folder into a single uncompressed FASTA file. The
latter because kraken2 cannot build from compressed FASTA files. It can
be done like this

``` bash
zcat $FNA/*.fna.gz > HumGut975_library.fna
```

The file `HumGut975_library.fna` should be close to 60GB.

Then you add this to the kraken2 database with

``` bash
kraken2-build --add-to-library HumGut975_library.fna --db $KRAKEN2
```

A subfolder `$KRAKEN2/library/added` should now appear. Inside it a
fasta-file (`.fna`) should appear, containing the library sequences,
i.e. should be of the same size as the `HumGut975_library.fna` or
whatever file you used.

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
read_delim("HumGut.tsv", delim = "\t") %>% 
  distinct(cluster95, .keep_all = T) -> humgut950.tbl
ok <- file.append("HumGut95_library.fna.gz",
                  file.path("FNA", humgut950.tbl$genome_file))
```

It basically filters the table in `HumGut.tsv` to keep only the rows
where you find the *first* occurrence of the unique values in the
`cluster95` column. These are the genomes representing the 95% clusters.
In `HumGut.tsv` the rows are sorted in descending order by the
`prevalence_score` and hence, the first row for each cluster is the
genome to keep.

Note that in the above code we assume the file `HumGut.tsv` and the
folder `FNA` (with all fasta files) is in the current working directory,
please use correct paths if they are elsewhere.

This code above produces a compressed FASTA file, and you need to
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
to also build a bracken and/or krakenUniq database (see below), as some
of these files will still be needed. <br> <br>

# Building a bracken database

Once you have the kraken2 database, it is straightforward to extend this
by building a bracken (<https://github.com/jenniferlu717/Bracken>)
database as well:

``` bash
bracken-build -d $KRAKEN2 -t 10 -k 35 -l 100
```

Here we used `10` threads, k-mers of length `35` (same as kraken2) and
read length `100`. This should add three files to the `$KRAKEN2` folder:
database.kraken, database.100mers.kraken,
database.100mers.kmer\_distrib. See the bracken software GitHub site for
more details: <https://github.com/jenniferlu717/Bracken>. <br> <br>

# Building a krakenUniq database

As an alternative to the kraken2/bracken approach, you may also be
interested in using the software krakenUniq
(<https://github.com/fbreitwieser/krakenuniq>).

Once you have the kraken2 database, it is straightforward to build a
krakenUniq database. Note that this requires much more memory than the
kraken2/bracken. For the full HumGut collection I allocate 500GB of
memory for this job, and you should have at least 600GB of free disk
space.

First, make a folder for you krakenUniq database, we denote this
`$KRAKENUNIQ`.

Next, copy both the `taxonomy/` and the `library/` folders from the
kraken2 database to the new folder:

``` bash
cp -r $KRAKEN2/taxonomy $KRAKENUNIQ
cp -r $KRAKEN2/library $KRAKENUNIQ
```

Finally, build the krakenUniq database:

``` bash
krakenuniq-build --threads 20 --db $KRAKENUNIQ
```

You may free a lot of disk space by running

``` bash
krakenUniq-build --clean -db $KRAKENUNIQ
```

In my case the full HumGut database was reduced to occupy around 350GB
disk space after this step. A similar cleaning may also be done for the
`KRAKEN2` database now.
