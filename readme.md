# Simple python package for storing and querying large Multi-Omics

This package aims to aggregates multiple modalities of data (currently only RNA-Seq and small variants (SNPs and inde;s)
are supported but more to come). The inspirations was the constant frustration of parsing files many time just to switch the 
analysis from one modality to another. I decided to create a database schema that would support large (-ish) studides
with multi-omic datasets. The package does not perfrom any kind of analysis (yet). It is only there for data retrival. 

## Create a project database
`create_project.py` creates a sqlite database that contains several linked tables. The currently supported tables are

+ Samples (this is the sample metadata)
+ gene expression (currently RSEM output will make it more flexible later, see `vcf.yaml`)
+ isoform expression (same as above)
+ Junctions (this is STAR's SJ.out.tab file again same with expression will include more fiexiblity)
+ Filtered junctions (same as junctions but these are the "interesting" ones)
+ RNA-Variants (if you have performed varriant calling from your RNA-Seq data)
+ Filtered RNA-Variant (same as above)
+ Variants
+ Filtered Variants

For junctions, variants and their filtered counterparts you will need one file per sample (it's perfecly ok if a sample is 
missing one or more of the files)

To create the database you some additional files in addition to the processed data above. 

1. A samples csv with all the metadata you need (see `config.yaml`)
2. files file for each data modality where processed files are mapped to sample ids
3. A `config.yaml` (see below)
4. A `vcf.yaml` that describes the vcf fields you are intereted in my go-to is usually ensembl VEP (see below). 

### config.yaml

This file describes what kind of files and data modalities are present in the study. Please take a look at the example file
provided. It should be fairly self explanatory

### vcf.yaml

Same as above, this files describes the vcf files and the fields in the INFO field that you are interested in. Currently I'm only
supporting one but supporting multiple is definitely possible by re-running the same functions for different fields.

You can create the databse as such:

```bash
python3 create_project.py -y config.yaml
```

This script will create the project database, if you have set the create flag in the `config.yaml` then samples will be added to the database
instead of creating a new one. Keep in mind that if you have provided a samples files this needs to have the NEW SAMPLES ONLY since trying to 
add existing samples will break the primary key constraint on the samples table. 

## Querying the database

First you need to create a genome class from [pytxdb](https://github.com/celalp/pytxdb) and use that in your project class

```python
from sqlalchemy import create_engine
import pybiomart as biomart
from pytxdb import Genome
from clinpy import Project, Junction

db=create_engine("sqlite:///genome.db")
mart=dataset = biomart.Dataset(name='hsapiens_gene_ensembl',
                              host='http://grch37.ensembl.org')

genome=Genome(db, "hs37d5_spikein_fixed.fa", mart)

project_db=create_engine("sqlite:///project.db")
project=Project(project_db, genome)
```

The `Project` class has several methods:

+ `Project.samples()` is for searching samples by cohort or samplename returns a dataframe
+ `Project.expression()` is for searching expression level by gene or transcript, it can return long or wide format dataframe
+ `Project.junctions()` is for searching junctions (filtered only at the moment) it can return a simple dataframe or a list of 
`Junction()` instances (see below)

## Junction class

While dealing with expression levels is somewhat straighforward dealing with splice junctions in short reads is not. Currently, 
there are no packages that manages to do this elegantly (not that mine is any better at the moment). 

The junction class has many methods that I think will make it easier to deal with some of the quirks that I face. 

You can create a `Junction` class manually or return a list of junctions from the `Project` instance. 

A manual creation will look like this:
```python
junc=Junction(chrom="1", start=1581906, end=1645137, strand="-", uniq_map=11, 
              multi_map=1)
```

As you can see we are not adding the genome and project information in the junction class. I am not sure if that is a good 
idea but doing so creates a lot of duplication of the Juction and Genome classes.

### Using the `Junction` class:

For each junction instance there are many methods one can use some are rather self explanatory:

+ `Junction.genes()` search for genes that span the junction multiple genes are supported
+ `Junction.transcripts()` same as above not all transcript of the genes might be fetched because some isoforms are 
larger than others. 
+ `Junction.features()` instead of returning the whole span of the junction this one only returns the exon/intron of the
transcript(s) of interest for the start/end locations
+ `Junction.samples()` return a list of samplename with that have the same junction, you can add some tolerance either in 5' or 3' 
or specify % overlap (reciprocal or not)
+ `Junction.new_transcript()` you can get the sequence of the new transcript with that contains the junciton. It takes the output
of features. 
+ I am currently working on adding a `Junctions.filter()` method where you can specify a function to filter the junctions file
Then the method will go to the junctions table and iterate over samples one by one and create the filtered junctions tables. (or not 
there will be several options for different results)

## Variants class (under development)

This will have methods to query and filter variants wheter they are called form RNA-Seq data or not.  It will have methods to 
display available info and qual/gt columns and an additiona function to filter variants based on arbitrary constraints on the former. 

I am working on these methods as we speak and will try to make them as flexible as possible. 

## Variant class

This class represent a single variant. You will have methods to search for samples that contain the variant or allele freq or allele
count of the variant. You can specify a list of samples or a cohort so you can limit to search to a subset of samples. 

## in the future
This is nowhere near complete or tested in some order

1. Create a more flexible way to add expression and splicing data
2. add support for structural variants
3. copy number variation
4. repeat expansions
5. allele specific expression
6. Ability to give user to design their own assay and how that relates to each sample