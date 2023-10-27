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

max_af=$4

type=$5

alt_depth=3

gemini query -q "select name from samples order by name" $file > samples.txt

#if pipeline is cre, filter out variants only called by one of freebayes, samtools, platypus
#else pipeline is mosaic/crg (uses one caller), do not filter by caller, i.e. no "callers" in the gemini db
callers=`gemini db_info $file | grep -w "variants" | grep -w "callers"` 
if [ ! -z "$callers" ] #variable $callers is not an empty string, i.e. it exists in the gemini db
then
	callers="callers"
	caller_filter="and Callers not in ('freebayes', 'samtools', 'platypus')"
else	
	callers="00"
	caller_filter=""
fi

if [[ "$type" == 'wgs' || "$type" == 'denovo' ]]
then
    noncoding_anno="uce_100bp as UCE_100bp, uce_200bp as UCE_200bp,
            dnasei_hypersensitive_site as DNaseI_hypersensitive_site,
            ctcf_binding_site as CTCF_binding_site, 
            enh_cellline_tissue as ENH_cellline_tissue,
            tf_binding_sites as TF_binding_sites,
            c4r_wgs_counts as C4R_WGS_counts,
            c4r_wgs_samples as C4R_WGS_samples"
else
    noncoding_anno="00 as noncoding"
fi

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
		COALESCE(clinvar_pathogenic, '') || COALESCE( ';' || NULLIF(clinvar_sig,''), '') as Clinvar, \
        clinvar_status as Clinvar_status, \
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
	gnomad_male_ac as Gnomad_male_ac, \
        hprc_af as HPRC_af, \
        hprc_ac as HPRC_ac, \
        hprc_hom as HPRC_hom, \
        cmh_af as CMH_af, \
        cmh_ac as CMH_ac, \
        sift_score as Sift_score,\
        polyphen_score as Polyphen_score,\
        cadd_phred as Cadd_score,\
        vest4_score as Vest4_score,\
        revel_score as Revel_score,\
        gerp_score as Gerp_score,\
        aa_change as AA_change,\
        hgvsc as Codon_change,\
        "$callers" as Callers,\
        phylop30way_mammalian as Conserved_in_30_mammals,\
        COALESCE(spliceai_score, '') as SpliceAI_score, \
        uce_100bp as UCE_100bp, uce_200bp as UCE_200bp, \
        Dark_genes as Dark_genes, \
        $noncoding_anno, \
        gts,"

while read sample
do
	sQuery=$sQuery"gts."$sample","
	sQuery=$sQuery"gt_alt_depths."$sample","
	sQuery=$sQuery"gt_depths."$sample","
    sQuery=$sQuery"gt_quals."$sample","
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
sQuery=$sQuery" where gnomad_af_popmax <= "${max_af}" "$caller_filter""${severity_filter}""

s_gt_filter=''
# denovo 0/1 is exported in cre.sh
if [ -n "$denovo" ] && [ "$denovo" == 1 ]
then
    # https://www.biostars.org/p/359117/
    proband=`gemini query -q "select name from samples where paternal_id != -9 and paternal_id != 0 and maternal_id != -9 and maternal_id != 0" $file`

    # print header
    header=$sQuery" limit 0"
    gemini query -q "$header" --header $file

    # get de novo variants for each affected proband
    for p in $proband
    do
        mom=`gemini query -q "select maternal_id from samples where name == '$p'" $file`
        dad=`gemini query -q "select paternal_id from samples where name == '$p'" $file`
        s_gt_filter="((gt_types."$p" == HET or gt_types."$p" == HOM_ALT) and gt_types."$dad" == HOM_REF and gt_types."$mom" == HOM_REF) \
        and (gt_alt_depths."$p" >="${alt_depth}" or (gt_alt_depths).(*).(==-1).(all)) \
        and ((gt_alt_depths."$dad" < 10 and gt_alt_depths."$mom" < 10)  or (gt_alt_depths).(*).(==-1).(all))"
        # otherwise a lot of trash variants
        sQuery=$sQuery" and qual>=400 and Old_multiallelic is null"
        gemini query -q "$sQuery" --gt-filter "$s_gt_filter" $file
    done
else
    # keep variant where the alt depth is >=3 in any one of the samples or they're all -1 (sometimes happens for freebayes called variants?)
    s_gt_filter="(gt_alt_depths).(*).(>="${alt_depth}").(any) or (gt_alt_depths).(*).(==-1).(all)"
	gemini query -q "$sQuery" --gt-filter "${s_gt_filter}" --header $file

    # also get the clinvar variants (duplicates will be removed later)
    cQuery=$initialQuery
    cQuery=$cQuery" where gnomad_af_popmax <= ${max_af} "$caller_filter" and Clinvar <> ''"
    # only get variants where AD >= 1 (any sample with an alternate read)
    c_gt_filter="(gt_alt_depths).(*).(>=1).(any) or (gt_alt_depths).(*).(==-1).(all)"
    gemini query -q "$cQuery" --gt-filter "$c_gt_filter" $file

    # if allele frequency is > 1% and Clinvar is pathogenic, likely pathogenic or conflicting and any status except for no assertion
    cQuery=$initialQuery
    cQuery=$cQuery" where gnomad_af_popmax > ${max_af} "$caller_filter" and Clinvar_status != 'no_assertion_criteria_provided' and Clinvar in ('Pathogenic', 'Likely_pathogenic', 'Conflicting_interpretations_of_pathogenicity')"
    gemini query -q "$cQuery" $file

fi
