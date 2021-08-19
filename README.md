# cre
Excel variant report generator and scripts to process WES data (cram/bam/fastq -> variant calls -> annotated variant calls -> prioritized variants -> excel report). Uses results from [bcbio variant2](https://bcbio-nextgen.readthedocs.io/en/latest/contents/pipelines.html#germline-variant-calling) germline variant calling pipeline.

# 1. Installation

1.1 Install **bcbio-nextgen**.
*  [Installation manual](https://bcbio-nextgen.readthedocs.io/en/latest/contents/installation.html)
*  [Installation and update examples](https://github.com/naumenko-sa/cre/blob/master/cre.bcbio.upgrade.sh)

Use HPC or server. Add bcbio to PATH and PYTHONPATH in .bash_profile:
```
export PATH=[installation_path]/tools/bcbio/bin:[installation_path]/tools/bcbio/anaconda/bin:$PATH
export PYTHONPATH=[installation_path]/tools/bcbio/anaconda/lib/python2.7:$PYTHONPATH
```
bcbio installs many other useful tools (including java and R) and datasets through bioconda and cloudbiolinux.
  
1.2 Clone **cre** to ~/cre and add it to PATH.
```
cd
git clone https://github.com/naumenko-sa/cre
```
.bash_profile: `export PATH=~/cre:$PATH`

1.3 (Optional) Install/update OMIM.

By default CRE uses [cre/data/omim.txt](../master/data/omim.txt) and [cre/data/omim.inheritance.csv](../master/data/omim.inheritance.csv).
* Goto https://omim.org/downloads/ and request the latest database.
* In a couple of days you will receive: genemap2.txt, genemap.txt, mim2gene.txt, mimTitles.percent.txt, mimTitles.txt, morbidmap.txt. 
* Preprocess OMIM with [cre.omim.sh](../master/cre.omim.sh): `cd OMIM_DIR;~/cre/cre.omim.sh`
  
Result - omim.txt with omim description of diseases related to ~ 4000 genes
We use the improved OMIM inheritance table from [https://www.cs.toronto.edu/~buske/cheo/](https://www.cs.toronto.edu/~buske/cheo/).Download the second file with inheritance mappings. It references genes by gene name (symbol) rather than by Ensembl_id which is a requirement for CRE. Most gene names (symbols) could be mapped automatically with Ensembl biomart [genes.R](https://github.com/naumenko-sa/bioscripts/blob/master/genes.R), but some genes (not many) might need manual curation to assign the correct ENSEMBL_ID.

1.4 (Optional) Install/update Orphanet: `cd ~/cre/data; ~/cre/cre.orphanet.sh`
Orphanet provides descriptions for ~3600 genes:. By default CRE uses [orphanet.txt](../master/data/orphanet.txt)
1.5 (Optional) Update Gnomad gene contraint scores: `Rscript ~/cre/cre.gnomad_scores.R`
By default using [~/cre/data/gnomad_scores.csv](../master/data/gnomad_scores.csv)

1.6 (Optional) Imprinted genes.
By default using [~/cre/data/imprinting.txt](../master/data/imprinting.txt).

1.7 (Optional) Install HGMD pro database
Install HGMD pro and dump information with [~/cre/cre.hgmd2csv.sql](../master/cre.hgmd2csv.sql).

# 2. (optional) Alignment to grch37 with decoy

By default, bcbio does not have decoy in grch37 reference, decoy is supported only in grch38. Using decoy improves FDR by ~0.5%. Two step approach could be applied to use decoy in bcbio https://github.com/bcbio/bcbio-nextgen/issues/2489:
- install custom grch37d5 reference with decoy:  [cre.bcbio.custom_genome.sh](../master/cre.bcbio.custom_genome.sh)
- run alignment step vs grch37d5 reference: `cre.prepare_bcbio_run.sh <project> align_decoy`
- keep bam file aligned vs grch37d5 for storage
- run variant calling with noalt_calling and bam_clean: remove_extracontigs (SV calling in WGS required additional processing of decoy aligned bam file, see crg).
- when saving project to the storage, use bam files aligned to decoy

# 3. Set up bcbio project for alignment, variant caling and annotation

* Prepare input files: family_sample_1.fq.gz, family_sample_2.fq.gz, or family_sample.bam and place them into family/input folder.
* There might be many samples in a family(project).
* TEST: NIST Ashkenazim trio: ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio (download OsloUniversityHospital exomes).
* run [cre.prepare_bcbio_run.sh](../master/cre.prepare_bcbio_run.sh) [family].

By default it uses [bcbio.templates.wes.yaml](../master/bcbio.templates.wes.yaml) config with following features:
  * 4 callers, the order is important, because variant metrics in ensemble calling (like AD) are picked up from the first caller in the list (gatk)
  * ensemble calling
  * realignment and recalibration. There is no much gain in precision/sensitivity with RR, but to make bam files consistent with other projects it is on here. Actually, realignment helps samtools to call indels better.
  * no bed file. Let callers call every variant which has coverage, we will filter poorly covered variants later. Modern exome capture kits are so perfect, that
we can discover a useful non-coding variant. No sense to filter them out during that stage.
  * effects: VEP. There is a holywar VEP/snpEff/Annovar. My choice is VEP. You will see later both Ensembl and Refseq in the report, so no reason for using Annovar.
  * effects_transcripts: all. We want all effects of a variant on all transcripts to be reported.
  * aligner: bwa. Even staring with bam files, bwa is used. Sometimes input bam files aligned against older reference, or different (chr) naming scheme. It is better to have a bam file consistent with calls made.
  * custom annotation [cre.vcfanno.conf](../master/cre.vcfanno.conf) using data sources installed in bcbio.
  * creates gemini database with vcf2db

## 3a. Input files are in Illumina basespace.
* use basespace-cli to dump bcl files to HPC, then do 1b.

## 3b. Input is Illumina run (bcl files).
* create a sample sheet and run [bcl2fq.sh](../master/bcl2fq.sh).

## 3c. Input is cram file.
* Run [cram2fastq_samtools.pbs](../master/cram2fastq_samtools.pbs) for each cram file.
* Usage: run `samtools view -H <cram_file>` to find the reference path as `old_ref`, then do `qsub cram2fastq_samtools.pbs -v cram=/path/to/cram,old_ref=reference path in cram file header,sample=familyid_sampleid,dir=output directory`. This script switches the reference used in cram to using our local reference, and then it generates fastq files.
* For eg. `qsub cram2fastq_samtools.pbs -v cram=/hpf/largeprojects/ccm_dccforge/dccforge/uploads/CHEO/2248_CH2188/2211891.cram,old_ref=UR:/mnt/hnas/reference/hg19/hg19.fa,sample=2248_CH2188,dir=/hpf/largeprojects/ccmbio/ccmmarvin_shared/exomes/in_progress/2248/input`  
* I would suggest to avoid crams when possible. A damaged bam file could be recovered with [cre.bam_recovery.sh](../master/cre.bam_recovery.sh), but nothing could be done for cram.

# 4. Run bcbio

* Single project: `qsub ~/cre/bcbio.pbs -v project=[project_name]`
Project should have a folder project_name in the current directory.

* Multiple projects: `qsub -t 1-N ~/cre/bcbio.array.pbs`
Current directory should have a list of projects in projects.txt.

# 5. Clean result dir and create project.csv report: 

`qsub ~/cre/cre.sh -v family=[family],cleanup=1`
 * moves project results and sample bam files to family dir
 * removes work and final dirs from bcbio project
 * removes gemini databases for individual callers (we need only ensemble gemini database)

  During the report generation step:
  * dumps variants from gemini database to tab text file
  * dumps variant impacts from gemini database to tab text file
  * annotates variants with refseq in addition to ensembl
  * gets coverage from GATK Haplotype calls, freebayes, and platypus
  * build excel report based on gemini variants table, variant impacts, coverage information and some other fields.

# 6. Step 5 in detail

6.1 [Report description](https://sickkidsca-my.sharepoint.com/:w:/r/personal/arun_ramani_sickkids_ca/_layouts/15/Doc.aspx?sourcedoc=%7B24220052-BF91-4B96-989E-E3143F964A39%7D&file=SNV_FilterColumnDescription.May.2021.docx&action=default&mobileredirect=true).\
6.2 [Report example for Ashkenazim trio from NIST](https://drive.google.com/open?id=0B_bLL10GwDnsN29vY3RRdGlXMWM). \
6.3 [gemini.gemini2txt.sh](../master/gemini.gemini2txt.sh) [project-ensembl.db] - dumps a gemini database into text file. \
6.4 [gemini.variant_impacts.sh](../master/gemini.variant_impacts.sh) [project-ensembl.db] dumps variant impacts from gemini. \
6.5 creates a vcf file with rare and potentially deleterious variants, the same set of variants is in the excel report.
```
cat ${family}-ensemble.db.txt | cut -f 23,24  | sed 1d | sed s/chr// > ${family}-ensemble.db.txt.positions
    bcftools view -R ${family}-ensemble.db.txt.positions -o ${family}.vcf.gz -O z ${family}-ensemble-annotated-decomposed.vcf.gz

```
6.6 coverage from VCFs produced by GATK, platypus, and freebayes - requires gatk wrapper from bcbio.
```
vcf.freebayes.getAO.sh ${family}-freebayes-annotated-decomposed.vcf.gz
vcf.gatk.get_depth.sh ${family}-gatk-haplotype-annotated-decomposed.vcf.gz
vcf.platypus.getNV.sh ${family}-platypus-annotated-decomposed.vcf.gz

```
6.7 Rscript ~/cre/[cre.R](../master/cre.R) [family] - creates report family.csv.

# 7. How to create a database of variants

- cre.database.sh [input_dir] [output_dir] - creates sample-wise and variant-wise reports, which are necessary for annotation with cre.R.
- cre.database.pull_gene.sh [database_prefix] [gene_name] - pull a gene report from the database.

# 8. Coverage plots

- ~/bioscripts/genes.R - pull a bed file from Ensembl
- ~/bioscripts/bam.coverage.bamstats05.sh - calculate coverage
- cheo.R - plot coverage pictures

# 9. List of all scripts

* bcbio.array.pbs
* bcbio.pbs
* bcbio.prepare_families.sh
* bcbio.rename_old_names.sh
* bcl2fastq.sh
* cheo.R - mostly Venn diagrams to compare pipeline validations + some coverage analysis
* cre.prepare_bcbio_run.sh
* cre.bam.validate.sh
* cre.bam.remove_decoy_reads.sh - removes decoy reads for grch37d5 with VariantBam and grep.
* cre.bcbio.upgrade.sh - examples of bcbio installation and upgrade
* cre.coverage.bamstats05.sh - calculate coverage
* cre.fixit.sh - fixes sample names
* cre.gemini_load.sh loads vep-annotated vcf to gemini db.
* cre.gemini.get_variants4gene.sh - pull all varaints for a specific gene.
* cre.gnomad_scores.R - download and parse gnomad scores.
* cre.immunopanels.R - annotates CRE report with 6 immunopanels.
* cre.kinship.R - to plot relatedness (kinship) diagram for a group of samples. Sometimes helps to detect and solve mislabelling.
* cre.package.sh
* cre.rtg.validate.sh - validates NA12878 calls vs genome in a bottle callset with RTG and a bed file
* cre.sh - master script to produce variant reports from bcbio output
* cre.topmed.R - pull variant frequency from TopMed having rs_id
* cre.roh.h3m2.sh: a robust method of ROH/LOH analysis with h3m2, calls variants, accounts for exonic regions, LD, plots picture.
* cre.roh.naive.sh: retrieves MAF<5% high quality variants from gemini.db and reports stretches of >9 HOM variants (start, end, length, genes).
* cre.vcf.has2dp.sh fixes input vcf file from HAS pipeline (Illumina, TCAG) filling DP field
* omim.sh
* orphanet.sh
* vcf.freebayes.getAO.sh
* vcf.gatk.get_depth.sh
* vcf.platypus.getNV.sh
* vcf.samtools.get_depth.sh
* vcf.split_multi.sh
* vep4seqr_hg38.sh
* vep4seqr.sh

# 10. Credits

This work was inspired by 
* [bcbio](https://github.com/chapmanb/bcbio-nextgen/) and [gemini](https://github.com/arq5x/gemini) teams. Thank you all!
* Kristin Kernohan from Children's Hospital of Eastern Ontario (CHEO), who generated most ideas about the report contents. Thank you, Kristin, for all of the discussions!

Thank you colleagues at [CCM](https://ccm.sickkids.ca/), for seminars and personal discussions.
