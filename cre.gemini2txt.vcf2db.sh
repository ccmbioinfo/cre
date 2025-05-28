#!/bin/bash
#  exports gemini.db database to gemini.db.txt file
#  database schema: https://gemini.readthedocs.io/en/latest/content/database_schema.html#the-variants-table
#  when using v.chr = g.chr AND v.gene = g.gene it becomes very slow
#  by default bcbio writes PASS only variants to the database

#  example call: cre.gemini2txt.sh S28-ensemble.db 5 ALL
#  when using vcfanno/vcfdb loader some fields are different
#  for some reason \n in the query string does not work here

#SBATCH --time=1:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=10G

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
callers=`gemini db_info $file | grep -w "variants" | grep -w "callers"` 
if [ ! -z "$callers" ]
then
	callers="callers"
	caller_filter="and Callers not in ('freebayes', 'samtools', 'platypus')"
else	
	callers="00"
	caller_filter=""
fi

if [[ "$type" == 'wgs' || "$type" == 'denovo' ]]
then
    noncoding_anno="GreenDB_variant_type as GreenDB_variant_type,
            GreenDB_closest_gene as GreenDB_closest_gene,
            GreenDB_controlled_gene as GreenDB_controlled_gene"
    noncoding_scores="ncER as ncER_score, ReMM as ReMM_score, LinSight_Score as LINSIGHT_score"
else
    noncoding_anno="00 as noncoding"
    noncoding_scores="00 as noncoding_scores"
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
        gnomad_exome_af as gnomad_exome_af, \
        gnomad_exome_nhomalt as gnomad_exome_nhomalt,\
        gnomad_exome_af_popmax as gnomad_exome_af_popmax,\
        gnomad_exome_af_nfe as gnomad_exome_af_nfe,\
        gnomad_exome_af_afr as gnomad_exome_af_afr,\
        gnomad_exome_af_amr as gnomad_exome_af_amr,\
        gnomad_exome_af_eas as gnomad_exome_af_eas,\
        gnomad_exome_af_asj as gnomad_exome_af_asj, \
        gnomad_exome_af_fin as gnomad_exome_af_fin, \
        gnomad_exome_af_oth as gnomad_exome_af_oth, \
        gnomad_genome_af as gnomad_genome_af,\
        gnomad_genome_nhomalt as gnomad_genome_nhomalt,\
        gnomad_genome_af_popmax as gnomad_genome_af_popmax,\
        gnomad_genome_af_nfe as gnomad_genome_af_nfe,\
        gnomad_genome_af_afr as gnomad_genome_af_afr,\
        gnomad_genome_af_amr as gnomad_genome_af_amr,\
        gnomad_genome_af_eas as gnomad_genome_af_eas,\
        gnomad_genome_af_asj as gnomad_genome_af_asj,\
        gnomad_genome_af_fin as gnomad_genome_af_fin,\
        gnomad_genome_af_oth as gnomad_genome_af_oth,\
        af_1kg as af_1kg, \
        af_afr_1kg as af_afr_1kg, \
        af_amr_1kg as af_amr_1kg, \
        af_eas_1kg as af_eas_1kg, \
        af_eur_1kg as af_eur_1kg, \
        af_sas_1kg as af_sas_1kg, \
        sift_score as Sift_score,\
        polyphen_score as Polyphen_score,\
        cadd_phred as Cadd_score,\
        vest4_score as Vest4_score,\
        revel_score as Revel_score,\
        gerp_score as Gerp_score,\
        AlphaMissense as AlphaMissense,\
        $noncoding_scores, \
        aa_change as AA_change,\
        hgvsc as Codon_change,\
        "$callers" as Callers,\
        phylop30way_mammalian as Conserved_in_30_mammals,\
        COALESCE(spliceai_score, '') as SpliceAI_score, \
        $noncoding_anno, \
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

#max_aaf_all frequency is from gemini.conf and does not include gnomad WGS frequency, gnomad WES only
#gnomad_genome_af_popmax includes gnomad WGS. Use genome af_popmax to filter 
sQuery=$sQuery" where gnomad_genome_af_popmax <= "${max_af}" "$caller_filter""${severity_filter}""

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
    cQuery=$cQuery" where gnomad_genome_af_popmax <= ${max_af} "$caller_filter" and Clinvar <> ''"
    # only get variants where AD >= 1 (any sample with an alternate read)
    c_gt_filter="(gt_alt_depths).(*).(>=1).(any) or (gt_alt_depths).(*).(==-1).(all)"
    gemini query -q "$cQuery" --gt-filter "$c_gt_filter" $file

    # if allele frequency is > 1% and Clinvar is pathogenic, likely pathogenic or conflicting and any status except for no assertion
    cQuery=$initialQuery
    cQuery=$cQuery" where gnomad_genome_af_popmax > ${max_af} "$caller_filter" and Clinvar_status != 'no_assertion_criteria_provided' and Clinvar in ('Pathogenic', 'Likely_pathogenic', 'Conflicting_interpretations_of_pathogenicity')"
    gemini query -q "$cQuery" $file

fi
