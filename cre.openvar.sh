#!/bin/bash
####################################################################################################
#   creates csv report for small variants that impact alternative ORFs as called by OpenVar

#   parameters:
#	family = [family_id] (=project_id=case_id=folder_name, main result file should be family/family-ensemble.db)
#	database = path to folder where c4r count files and hgmd.csv are found
#   reference = path to reference fasta; default is /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/seq/GRCh37.fa
#   cre_path = path to cre repo; default is ~/cre
####################################################################################################

#SBATCH --job-name=openvar_report
#SBATCH --time=72:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --mem=4G
#SBATCH --output=%x-%j.out


function f_make_report
{
    cd $family
    $cre/cre.gemini2txt.vcf2db.sh ${family}-ensemble.db $depth_threshold $severity_filter $max_af $type > $family.variants.all.txt
    $cre/cre.gemini.variant_impacts.vcf2db.sh ${family}-ensemble.db $depth_threshold $severity_filter $max_af > $family.variant_impacts.all.txt

    # remove duplicate lines from results of gemini query
    awk '!a[$0]++' $family.variants.all.txt > $family.variants.txt
    awk '!a[$0]++' $family.variant_impacts.all.txt > $family.variant_impacts.txt

    tail -n +2 $family.variants.txt > $family.variants.trim.txt
    mv $family.variants.trim.txt $family.variants.txt

    cd ..


    echo GENERATING REPORT WITH TYPE: "${type}"
    Rscript $cre/cre.vcf2db.R $family "${type}" "${database}" "${cre}/data"
    
}

if [ -z $family ]
then
    family=$1
fi

echo $family


#default database path
if [ -z $database ]
then
    database="/hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/GRCh37/variation"
fi


# set path to cre
if [ -z $cre ]
then
    cre=~/cre
fi

if [ -z $max_af ]
then
    max_af=0.01
fi

export max_af
export depth_threshold=10
export severity_filter=ALL
export type=aORF


f_make_report


