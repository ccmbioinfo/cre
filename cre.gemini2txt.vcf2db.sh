#!/bin/bash
#  exports gemini.db database to gemini.db.txt file
#  database schema: https://gemini.readthedocs.io/en/latest/content/database_schema.html#the-variants-table
#  when using v.chr = g.chr AND v.gene = g.gene it becomes very slow
#  by default bcbio writes PASS only variants to the database

#  example call: cre.gemini2txt.sh S28-ensemble.db 5 ALL
#  when using vcfanno/vcfdb loader some fields are different
#  for some reason \n in the query string does not work here

#PBS -l walltime=1:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=10g,mem=10g

if [ -z $file ]
then
    file=$1
fi

#10 reads for WES 5 reads for RNA-seq
depth_threshold=$2

severity_threshold=$3
#echo $severity_threshold

max_af=$4

gemini query -q "select name from samples order by name" $file > samples.txt

sQuery="select \
        chrom as Chrom,\
        start+1 as Pos,\
        variant_id as Variant_id,\
        ref as Ref,\
        alt as Alt,\
        impact as Variation,\
        dp as Depth,\
        qual as Quality,\
        gene as Gene,\
        clinvar_pathogenic as Clinvar,\
        ensembl_gene_id as Ensembl_gene_id,\
        transcript as Ensembl_transcript_id,\
        aa_length as AA_position,\
        exon as Exon,\
        domains as Protein_domains,\
        rs_ids as rsIDs,\
        gnomad_af as Gnomad_af,\
        gnomad_af_popmax as Gnomad_af_popmax,\
        gnomad_ac as Gnomad_ac,\
        gnomad_hom as Gnomad_hom,\
        sift_score as Sift_score,\
        polyphen_score as Polyphen_score,\
        cadd_phred as Cadd_score,\
        vest3_score as Vest3_score,\
        revel_score as Revel_score,\
        gerp_score as Gerp_score,\
        aa_change as AA_change,\
        hgvsc as Codon_change,\
        phylop20way_mammalian as Conserved_in_20_mammals,
        gts,"

while read sample
do
	sQuery=$sQuery"gts."$sample","
	sQuery=$sQuery"gt_alt_depths."$sample","
	sQuery=$sQuery"gt_depths."$sample","
done < samples.txt

# gene_detailed may contain 2 records per single transcript - because of synonymous gene names, and some genes may have None in the name,for example TSRM
# https://groups.google.com/forum/#!topic/gemini-variation/U3uEvWCzuQo
# v.depth = 'None' see https://github.com/chapmanb/bcbio-nextgen/issues/1894

if [[ "$severity_threshold" == "ALL" ]]
then
#used for RNA-seq = 20k variants in the report
    severity_filter=""
elif [[ "$severity_threshold" == "wes.synonymous" ]]
then
    severity_filter=" and (is_coding=1 or is_splicing=1)"
#for WES = 1k variants in the report
else
    severity_filter=" and impact_severity<>'LOW' "
fi

#max_aaf_all frequency is from gemini.conf and does not include gnomad WGS frequencing, gnomad WES only
#gnomad_af includes gnomad WGS
sQuery=$sQuery"hgvsc as Nucleotide_change_ensembl,\
		hgvsp as Protein_change_ensembl,\
        old_multiallelic as Old_multiallelic
		from variants \
        where \
	        (dp >= "$depth_threshold" or dp = '' or dp is null) "$severity_filter" and gnomad_af_popmax <= "$max_af

s_gt_filter=''
if [ -n "$denovo" ] && [ "$denovo" == 1 ]
then
    # https://www.biostars.org/p/359117/
    proband=`gemini query -q "select name from samples where phenotype=2" $file`
    mom=`gemini query -q "select name from samples where phenotype=1 and sex=2" $file`
    dad=`gemini query -q "select name from samples where phenotype=1 and sex=1" $file`
    
		s_gt_filter="gt_types."$proband" == 1 and (gt_types."$dad" == 0 or gt_types."$dad" == 3) and (gt_types."$mom" == 0 or gt_types."$mom" == 3)"
    
		# otherwise a lot of trash variants
    sQuery=$sQuery" and qual>=500"
    gemini query -q "$sQuery" --gt-filter "$s_gt_filter" --header $file
else
    gemini query --header -q "$sQuery" $file
fi
