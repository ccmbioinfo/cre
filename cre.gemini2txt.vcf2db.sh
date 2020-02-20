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
alt_depth_3=$5
keep_clinvar=$6

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

if [[ "$severity_threshold" == 'ALL' || "$severity_threshold" == "wes.synonymous" ]]
then
#used for RNA-seq = 20k variants in the report
    severity_filter=""
#for WES = 1k variants in the report
else
    severity_filter=" and impact_severity<>'LOW' "
fi


sQuery=$sQuery"hgvsc as Nucleotide_change_ensembl,\
        hgvsp as Protein_change_ensembl,\
        old_multiallelic as Old_multiallelic from variants"

initialQuery=$sQuery # keep the field selection part for later use

#max_aaf_all frequency is from gemini.conf and does not include gnomad WGS frequencing, gnomad WES only
#gnomad_af includes gnomad WGS
if [ $alt_depth_3 -ne 1 ]
then
    sQuery=$sQuery" where (dp >= "${depth_threshold}" or dp = '' or dp is null) "${severity_filter}" and gnomad_af_popmax <= "${max_af}""
else
    sQuery=$sQuery" where gnomad_af_popmax <= "${max_af}" "${severity_filter}""
fi

s_gt_filter=''
if [ -n "$denovo" ] && [ "$denovo" == 1 ]
then
    # https://www.biostars.org/p/359117/
    proband=`gemini query -q "select name from samples where phenotype=2" $file`
    mom=`gemini query -q "select name from samples where phenotype=1 and sex=2" $file`
    dad=`gemini query -q "select name from samples where phenotype=1 and sex=1" $file`
    
    s_gt_filter="gt_types."$proband" == HET and gt_types."$dad" == HOM_REF and gt_types."$mom" == HOM_REF"
    #(gt_types."$proband" == HOM_ALT and gt_types."$dad" == HOM_REF and gt_types."$mom" == HET)"
    # otherwise a lot of trash variants
    sQuery=$sQuery" and qual>=500"
    gemini query -q "$sQuery" --gt-filter "$s_gt_filter" --header $file
elif [ -n "$alt_depth_3" ] && [ "$alt_depth_3" == 1 ]
then
    # keep variant where the alt depth is >=3 in any one of the samples
    s_gt_filter="(gt_alt_depths).(*).(>=3).(any)"
    gemini query -q "$sQuery" --gt-filter "${s_gt_filter}" --header $file
else
    gemini query --header -q "$sQuery" $file
fi

# if $keep_clinvar is set, complement selection to add clinvar-annotated variants that were excluded
# N.B. this will only work for below cases:
#   type is unset, alt_depth_3=0
#   type is unset, alt_depth_3=1
#   type=wes.synonymous,alt_depth_3=0
#   type=wes.synonymous,alt_depth_3=1
if [ -n "$keep_clinvar" ] && [ "$keep_clinvar" == 1 ]
then
    if [[ "$severity_threshold" == "wes.synonymous" ]]
    then
        cQuery=$initialQuery # grab earlier field selection
        cQuery=$cQuery" where (\
        (gnomad_af_popmax > "${max_af}" or (is_coding=0 and is_splicing=0)) \
        or (gnomad_af_popmax <= "${max_af}" and (is_coding=0 and is_splicing=0)) \
        or (gnomad_af_popmax > "${max_af}" and (is_coding=1 or is_splicing=0))\
        ) and clinvar_sig <> ''"
        if [ -n "$alt_depth_3" ] && [ "$alt_depth_3" == 1 ]
        then
            # here only want variants that were excluded previously (i.e all samples had < 3 alt depth)
            s_gt_filter="(gt_alt_depths).(*).(<3).(all)"
            #TODO: add Clinvar condition here
            gemini query -q "$cQuery" --gt-filter "$c_gt_filter" $file
        else
            # just need to add the clinvar filter
            gemini query -q "$cQuery" $file
        fi
    else
        # reversal of all previous filters (De morgan's laws for !(A && B && C) = (!A|!B|!C))
        cQuery=$initialQuery # grab earlier specs
        cQuery=$cQuery" where ((dp < "${depth_threshold}" or impact_severity == 'LOW' or gnomad_af_popmax > "${max_af}") and clinvar_sig <> '')"
        echo $cQuery
        if [ -n "$alt_depth_3" ] && [ "$alt_depth_3" == 1 ]
            then
                # here only want variants that were excluded previously (i.e all samples had < 3 alt depth)
                s_gt_filter="(gt_alt_depths).(*).(<3).(all)"
                gemini query -q "$cQuery" --gt-filter "$c_gt_filter" $file
        else
            # just need to add the clinvar filter
            gemini query -q "$cQuery" $file
        fi
    fi
fi