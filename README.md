# metagenomeannotation

nextflow pipeline for metagenome annotation in chunks with Bakta or Prokka

The assembly data will be split in chunks, annotated, simplified by gffred, and merged.

The pipeline is specifically made to annotate large assembles of [nf-core/mag](https://nf-co.re/mag) [Krakau et al. 2022](https://doi.org/10.1093/nargab/lqac007)
(from the [nf-core workflow collection](https://dx.doi.org/10.1038/s41587-020-0439-x)) 
and serve as input to for metatranscriptome analysis with [nf-core/rnaseq](https://nf-co.re/rnaseq).

## Software requirement

- Java
- nextflow
- singularity / apptainer

## Usage

### Input data

Input data is a fasta file, uncompressed or compressed (`*.gz`). In the current implementation, the input file name may not contain `.part_`.

### Run the pipeline

#### With Prokka

```bash
NXF_VER=24.04.4 nextflow run d4straub/metagenomeannotation -r main \
    --input MEGAHIT-group-all.contigs.head100.fa.gz \
    --nchunks 999 \
    --tool "prokka" \
    --outdir results
```

#### With Bakta

First prepare the database with

```bash
singularity pull  --name bakta-1.10.4--pyhdfd78af_0.img https://depot.galaxyproject.org/singularity/bakta:1.10.4--pyhdfd78af_0
singularity exec bakta-1.10.4--pyhdfd78af_0.img bakta_db download --type full --output bakta_db
```

and run the pipeline with

```bash
NXF_VER=24.04.4 nextflow run d4straub/metagenomeannotation -r main \
    --input MEGAHIT-group-all.contigs.head100.fa.gz \
    --nchunks 999 \
    --tool "bakta" \
    --baktadb bakta_db/db \
    --outdir results \
    -resume
```

### Tips & warnings

- The input file name may not contain `.part_`.

- Use a config to specify your cluster specifications, e.g. the here available `cfc_resources.config` by appending to the above command `-c cfc_resources.config`.

- The number of chunks can be as high as possible, but it was here tested only with 999, i.e. `--nchunks 999`, higher numbers might upset script input limits.

- Alternatively, the number of chunks can be omitted (i.e. not using `--nchunks`) but rather the size of chunks can be defined with `--chunksize`.

- For maximum execution speed make sure that `executor.queueSize` is at least as large as the number of chunks with `--nchunks`. The `cfc_resources.config` uses `executor.queueSize = 1000`.

## Output

- `cat/`
  - `<input>.<tool>.gtf`: annotation in gtf format
  - `<input>.<tool>.faa`: protein sequences
  - `<input>.<tool>.tsv`: tab-separated file with annotation details

Where `<input>` is the base name of the file supplied with `--input` and `<tool>` is "prokka" or "bakta".

## Credits

This pipeline was originally written by Daniel Straub ([@d4straub](https://github.com/d4straub)) for use at the [Quantitative Biology Center (QBiC)](http://www.qbic.life).

Code snippets (modules) were based on [nf-core modules](https://nf-co.re/modules/) of the [nf-core](https://nf-co.re) collection of workflows ([Ewels et al., 2020](https://dx.doi.org/10.1038/s41587-020-0439-x)).

## Citations

Please cite all employed tools, such as

- [Nextflow](https://pubmed.ncbi.nlm.nih.gov/28398311/)

  > Di Tommaso P, Chatzou M, Floden EW, Barja PP, Palumbo E, Notredame C. Nextflow enables reproducible computational workflows. Nat Biotechnol. 2017 Apr 11;35(4):316-319. doi: 10.1038/nbt.3820. PubMed PMID: 28398311.

- [Singularity](https://pubmed.ncbi.nlm.nih.gov/28494014/)

  > Kurtzer GM, Sochat V, Bauer MW. Singularity: Scientific containers for mobility of compute. PLoS One. 2017 May 11;12(5):e0177459. doi: 10.1371/journal.pone.0177459. eCollection 2017. PubMed PMID: 28494014; PubMed Central PMCID: PMC5426675.

- [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

- [seqkit](https://bioinf.shenwei.me/seqkit/)

  > Shen, W., Le, S., Li, Y., & Hu, F. (2016). SeqKit: A Cross-Platform and Ultrafast Toolkit for FASTA/Q File Manipulation. In Q. Zou (Ed.), PLOS ONE (Vol. 11, Issue 10, p. e0163962). Public Library of Science (PLoS). doi:10.1371/journal.pone.0163962

- [Prokka](https://pubmed.ncbi.nlm.nih.gov/24642063/)

  > Seemann T. Prokka: rapid prokaryotic genome annotation. Bioinformatics. 2014 Jul 15;30(14):2068-9. doi: 10.1093/bioinformatics/btu153. Epub 2014 Mar 18. PMID: 24642063.

- [Bakta](https://doi.org/10.1099/mgen.0.000685)

  > Schwengers O., Jelonek L., Dieckmann M. A., Beyvers S., Blom J., Goesmann A. (2021). Bakta: rapid and standardized annotation of bacterial genomes via alignment-free sequence identification. Microbial Genomics, 7(11). https://doi.org/10.1099/mgen.0.000685

- [GFFREAD](https://f1000research.com/articles/9-304)

  > Pertea G and Pertea M. GFF Utilities: GffRead and GffCompare [version 2; peer review: 3 approved]. F1000Research 2020, 9:304 (https://doi.org/10.12688/f1000research.23297.2) 

