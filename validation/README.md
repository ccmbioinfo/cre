# CRE Validation 
This folder contains scripts for validating updates to cre report generation. These changes may be any of the following:
1. Annotation updates
2. Tool version updates
3. Changes to filtering thresholds

We will make updates to cre every six months, and we will perform validation at both the report (csv) level, and vcf level.

## Report validation
1. Compare report generated prior to updates (first csv) to report generated post-updates (second csv): 
2. `qsub -F ~/cre/validation/compare_reports.pbs "<first_csv> <second_csv> <first_db> <second_db> </path/to/tables/first report> </path/to/tables/second report>"`
   
   where `<first_db>` refers to the gemini database associated with the first report, and `<second_db>` refers to the gemini database associated with the second report. `</path/to/tables/[first|second] report>` refers to the family/project directory within which the report was generated, where the *table files corresponding to GATK, freebayes, and platypus are found.
   To retrieve the updated/newly added column names from the db for comparison, you must manually add these column names as last arguments in the lines calling "gemini_compare.sh" inside the `compare_reports.pbs`
   
   Outputs:
    *  summary.txt: summarizes the number of variants shared between reports, and variants unique to each report. Also lists the column headers that differ between reports.
        
        ```
        Comparing reports: 428.wes.regular.2020-04-03.csv, 428.wes.regular.2020-10-30.csv
        Columns unique to 428.wes.regular.2020-04-03.csv: [u'omim_gene_description']
        Columns unique to 428.wes.regular.2020-10-30.csv: [u'omim_phenotype']
        Number of variants unique to 428.wes.regular.2020-04-03.csv: 198
        Number of variants unique to 428.wes.regular.2020-10-30.csv: 67
        Number of variants shared: 790
        ``` 
    *  <family_id>.different_annotations.csv: 


           annotations that differ as a set for variants shared between reports. Annotation columns for the first and second reports will be adjacent in the csv file, e.g. sift_score_new,sift_score_old,splicing_new,splicing_old,variation_new,variation_old. Check this file to make sure there are no unexpected changes to annotations, for instance the splicing column being filled with NA values as opposed to the expected numerical values.

    *  <first_csv>.uniq.pos, <second_csv>.uniq.pos:

           positions of variants unique to first and second reports respectively

    *  <first_csv>.uniq.db.txt, <second_csv>.uniq.db.txt:

           gemini database entries associated with both reports for variants unique to first report, and vice-versa

  1. Determine explanations for inclusion/exclusion of variants in first and second reports. Because the the updates to report generation will be different in each validation (for instance, we may update ClinVar annotations and VEP version in one validation, and then six months later update OMIM annotations and alt depth threshold), the validation script will change. For the October 2021 validation, run the validation_2021-10.py script as follows using python3 (first load python3 with `module load python/3.7.1`):
     ```python3
     python3 ~/cre/validation/validation_2021-10.py \
        -db_output1 <first_csv>.uniq.db.txt \
        -db_output2 <second_csv>.uniq.db.txt \
        -prefix1 <first_csv> \
        -prefix2 <second_csv> \
        -dir1 </path/to/tables/first report> \
        -dir2 /path/to/tables/first report
     ```

     Outputs:
      * validation_summary_unique_in_<first_csv>.<date_generated>.csv, validation_summary_unique_in_<second_csv>.<date_generated>.csv:
          
          These are tables that look like this:
          | Variant                                  | Explanation                                  |
          | ---------------------------------------- | -------------------------------------------- |
          | 15:38641773:38641774:A:T                 | Change in clinvar_sig from uncertain to None |
          | 19:633528:633529:G:GGCGCCGCCGCCGCCGCCGCC | Change in impact_severity from HIGH to LOW   |

   1. Some variants unique to one report may have the following explanation for their exclusion from the other report: "Alt depth less than 3 in other report tables". These variants have an apparent depth of >=3, and yet are not present in a report. The explanation is convoluted. When cre.sh is run, it creates [tables of the alt depths](https://github.com/ccmbioinfo/cre/blob/master/cre.sh#L209) for each variant from GATK, freebayes, and platypus. It doesn't create a table for samtools, because samtools doesn't always include the alt depth for a variant. When we query the gemini db, we select variants with alt depth >= 3 OR missing. The purpose of the alt depth tables seems to be to populate the alt depths in the report for variants where it is missing (-1) in the gemini database, i.e. for variants called by samtools and freebayes/platypus, so that the alt depth filter can be applied at the end of report generation to catch variants for which the AD was missing in the gemini db, but is actually <3. So, a variant may have AD missing in the gemini db if called by samtools and platypus, but we then [pull the AD from the platypus table](https://github.com/ccmbioinfo/cre/blob/master/cre.vcf2db.R#L455). If that AD is >=3, the variant will be included in the report. Note, however, that some samtools variants do have alt depth listed, which I believe is where we run into issues when doing the report validations. A variant may have an AD of >=3 in the gemini db which has been pulled from the samtools vcf, but is then filtered out in the report generation step if the AD is <3 in the other caller (e.g. platypus or freebayes) because cre does not generate a samtools table. So, we incorporated these tables into the validation_2021-10.py script to explain why a small number of variants with an apparent depth of >=3 in the gemini db from a run are not included in that report, as the AD was <3 in another caller.  

## VCF validation 

Validate pipeline performance using GIAB benchmark datasets. 
1. Test set: Download and run benchmark dataset with your pipeline 
   1. WGS: ftp:/­/­ftp.­sra.­ebi.­ac.­uk/­vol1/­run/­ERR323/­ERR3239334/­NA12878.­final.­cram
   2. WES: https://www.internationalgenome.org/data-portal/sample/NA12878 
2. Truth set: Download the latest variant calls (VCF) and high-confidence regions (BED) from [here](http://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/release/NA12878_HG001/latest/GRCh37/)
3. Use the following script to compare truth vs test. "bed" and "restrict" are optional, and their equivalent argument in vcfeval are,
   * bed: `-e, --evaluation-regions` 
   * restrict: `--bed-regions`
   ```bash
   qsub rtg-validations.pbs -v truth="<path to truth vcf>",test="<path to test vcf>",out="output folder name"[,bed="<path to high-conf bed>",restrict="<path to exomes.bed>"(optional)]
   ```

4. The truth set calls are from Whole-genome, and so there will be higher number of false negatives when benchmarking WES calls. You can do one of the below steps to restrict regions to exomes,
   1. intersect the giab high-confidence bed file with capture kit regions, and callable BED file from bcbio run. Pass this resulting file as `bed`
   to the script (see [here](https://github.com/bcbio/bcbio-nextgen/blob/747045809e493a5cca0dad0ec4ff053afafd6708/config/examples/NA12878.validate.sh))

### Outputs

* summary.txt: TP, FP, FN, Precision, Specificity, F-measure (use the 'None' row)
* fn.vcf.gz: false-negative variants
* fp.vcf.gz: false-positive variants
* tp.vcf.gz: true-positive variants in test set
* tp-baseline.vcf.gz: true-positive variants in truth set
* <prefix>_roc.png: roc created from weigthed_roc.tsv.gz
* <prefix>_pr.png: precision-recall curve created from weigthed_roc.tsv.gz
* <prefix>_snv-indel_pr.png: precision-recall curve created for snp/indel split using non_snp_roc.tsv.gz and snp_roc.tsv.gz
* <prefix>_snv-indel_roc.png: roc created for snp/indel split using non_snp_roc.tsv.gz and snp_roc.tsv.gz
  
You can also create a combined ROC/Precision-recall curve in single plot using any number of validation results like below,

```bash
rtg rocplot validation1/weigthed_roc.tsv.gz validation2/weigthed_roc.tsv.gz --png=validation12_roc.png
```
