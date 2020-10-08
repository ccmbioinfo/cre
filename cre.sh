#!/bin/bash
####################################################################################################
#   keeps only important files from bcbio run: qc, vcf, gemini, bam
#   creates csv report for small variants

#   parameters:
#	family = [family_id] (=project_id=case_id=folder_name, main result file should be family/family-ensemble.db)
#	cleanup= [0|1] default = 0
#	make_report=[0|1] default = 1, don't make report for WGS analysis first
# 	type = [ wes.regular (default) | wes.synonymous | wes.fast | rnaseq | wgs | annotate (only for cleaning) | 
# 	    denovo (all rare variants in wgs, proband should have phenotype=2, parents=phenotype1 also sex for parents in gemini.db) ]
#	max_af = af filter, default = 0.01
####################################################################################################

#PBS -l walltime=30:00:00,nodes=1:ppn=1
#PBS -joe .
#PBS -d .
#PBS -l vmem=40g,mem=40g

# default settings:
# max af > 0.1, ad > 3, plus and any variants with a clinvar entry and ad > 1

# cleanup is different for wes.fast template - don't remove gatk db
function f_cleanup
{
    # better to look for project-summary than hardcode the year
    # keep bam files for new samples
    if [ -z $family ] 
    then
        echo "Project (family) folder does not exist. Exiting"
        exit 1
    fi
    
    cd $family
    
    project_summary=`find final -name project-summary.yaml`
    
    # if there are multiple folders
    n_summaries=`echo $project_summary | grep -o yaml | wc -l`
    if [ $n_summaries -gt 1 ]
    then
	echo "Multiple project-summary.yaml"
	exit 1
    fi

    echo "Project summary: " $project_summary

    # if summary exists
    if [ -f $project_summary ] 
    then
        result_dir=`echo $project_summary | sed s/"\/project-summary.yaml"//`
        echo "result_dir =" $result_dir
    
        # if result_dir is empty that might cause copying entire /
        if [ -d $result_dir ] && [ -n "$result_dir" ]
        then
            mv $result_dir/* .
            mv final/*/*.bam .
            mv final/*/*.bai .
            # keep validation picture
            mv final/*/*.png .
	
            # keep sv calls
            if [ "$type" == "wgs" ]
            then
                mv final sv
            fi

            rm -rf final/
            rm -rf work/
    
	       #proceed only if there is a result dir
           #don't remove input files for new projects
           #rm -rf input/
           #rename bam files to match sample names
        
            for f in *ready.bam;do mv $f `echo $f | sed s/"-ready"//`;done;
            for f in *ready.bam.bai;do mv $f `echo $f | sed s/"-ready"//`;done;

	       #make bam files read only
            for f in *.bam;do chmod 444 $f;done;

	       #calculate md5 sums
            for f in *.bam;do md5sum $f > $f.md5;done;

	       #validate bam files
            for f in *.bam;do cre.bam.validate.sh $f;done;
    
            if [ "$type" == "wes.fast" ] || [ "$type" == "wgs" ]
            then
	           ln -s ${family}-gatk-haplotype.db ${family}-ensemble.db
	           ln -s ${family}-gatk-haplotype-annotated-decomposed.vcf.gz ${family}-ensemble-annotated-decomposed.vcf.gz
	           ln -s ${family}-gatk-haplotype-annotated-decomposed.vcf.gz.tbi ${family}-ensemble-annotated-decomposed.vcf.gz.tbi
            elif [ "$type" == "annotate" ]
            then
	           ln -s ${family}-precalled.db ${family}-ensemble.db
	           ln -s ${family}-precalled-annotated-decomposed.vcf.gz ${family}-ensemble-annotated-decomposed.vcf.gz
	           ln -s ${family}-precalled-annotated-decomposed.vcf.gz.tbi ${family}-ensemble-annotated-decomposed.vcf.gz.tbi
            else
	           # we don't need gemini databases for individual variant callers
	           rm ${family}-freebayes.db
	           rm ${family}-gatk-haplotype.db
	           rm ${family}-samtools.db
	           rm ${family}-platypus.db
	       fi
        fi
    fi
    cd ..
}

#workaround to fix: https://github.com/quinlan-lab/vcf2db/issues/52
#vcf2db changes - to _ in the sample names
function f_fix_sample_names
{
    bcftools query -l $1.subset.vcf.gz > $1.samples.txt
    if grep -q "-" $1.samples.txt;
    then
        cat $1.samples.txt | sed s/"-"/"_"/g > $1.samples.fixed.txt
        echo "VCF2DB fixed sample names, fixing sample names in " $1 " to match..."
        mv $1.subset.vcf.gz $1.subset.tmp.vcf.gz
        bcftools reheader -s $1.samples.fixed.txt $1.subset.tmp.vcf.gz > $1.subset.vcf.gz
        tabix $1.subset.vcf.gz
        rm $1.subset.tmp.vcf.gz
        rm $1.samples.fixed.txt
    fi
    rm $1.samples.txt
}

function f_make_report
{
    cd $family

    echo "Type="$type

    if [ "$type" == "rnaseq" ]
    then
	   export depth_threshold=5
    else
	   export depth_threshold=10
    fi

    if [ "$type" == "wgs" ] || [ "$type" == "rnaseq" ] || [ "$type" == "denovo" ] || [ "$type" == "wes.all" ]
    then
	   export severity_filter=ALL
    elif [ "$type" == "wes.synonymous" ]
    then
	   export severity_filter=$type
    else
	   export severity_filter=HIGHMED
    fi
    
    if [ "$type" == "denovo" ]
    then
			export denovo=1
    fi

    ~/cre/cre.gemini2txt.vcf2db.sh ${family}-ensemble.db $depth_threshold $severity_filter $max_af > $family.variants.all.txt
    ~/cre/cre.gemini.variant_impacts.vcf2db.sh ${family}-ensemble.db $depth_threshold $severity_filter $max_af > $family.variant_impacts.all.txt

    # remove duplicate lines from results of gemini query
    awk '!a[$0]++' $family.variants.all.txt > $family.variants.txt
    awk '!a[$0]++' $family.variant_impacts.all.txt > $family.variant_impacts.txt

    for f in *.vcf.gz
    do
	if [ ! -f $f.tbi ]
	then 
	   tabix $f
	fi
    done

    # report filtered vcf for import in phenotips
    # note that if there is a multiallelic SNP, with one rare allele and one frequent one, both will be reported in the VCF,
    # and just a rare one in the excel report
    cat $family.variants.txt | cut -f 1,2 | sed 1d | sed s/chr// | sort -k1,1 -k2,2n > ${family}-ensemble.db.txt.positions
    
    # this may produce duplicate records if two positions from positions file overlap with a variant 
    # (there are 2 positions and 2 overlapping variants, first reported twice)
    bcftools view -R ${family}-ensemble.db.txt.positions ${family}-ensemble-annotated-decomposed.vcf.gz | bcftools sort | vt uniq - | vt rminfo -t CSQ,Interpro_domain,MutPred_Top5features,MutationTaster_AAE - -o $family.vcf.gz
    tabix $family.vcf.gz
    
    rm $family.tmp.vcf.gz $family.tmp.vcf.gz.tbi

    #individual vcfs for uploading to phenome central
    vcf.split_multi.sh $family.vcf.gz

    reference=$(readlink -f `which bcbio_nextgen.py`)
    reference=`echo $reference | sed s/"anaconda\/bin\/bcbio_nextgen.py"/"genomes\/Hsapiens\/GRCh37\/seq\/GRCh37.fa"/`
    
    echo $reference

    vcf.ensemble.getCALLERS.sh $family.vcf.gz $reference

    #decompose first for the old version of bcbio!
    #gemini.decompose.sh ${family}-freebayes.vcf.gz
    fprefix=${family}-freebayes-annotated-decomposed
    if [ -f $fprefix.vcf.gz ]
    then
        bcftools view -R ${family}-ensemble.db.txt.positions $fprefix.vcf.gz | bcftools sort | vt decompose -s - | vt uniq - -o $fprefix.subset.vcf.gz
        tabix $fprefix.subset.vcf.gz
	
        f_fix_sample_names $fprefix
        vcf.freebayes.getAO.sh $fprefix.subset.vcf.gz $reference
    fi

    #gemini.decompose.sh ${family}-gatk-haplotype.vcf.gz
    fprefix=${family}-gatk-haplotype-annotated-decomposed
    if [ -f $fprefix.vcf.gz ]
    then
        bcftools view -R ${family}-ensemble.db.txt.positions $fprefix.vcf.gz | bcftools sort | vt decompose -s - | vt uniq - -o $fprefix.subset.vcf.gz
        tabix $fprefix.subset.vcf.gz
	
        f_fix_sample_names $fprefix
        vcf.gatk.get_depth.sh $fprefix.subset.vcf.gz $reference
    fi

    #gemini.decompose.sh ${family}-platypus.vcf.gz
    fprefix=${family}-platypus-annotated-decomposed
    if [ -f $fprefix.vcf.gz ]
    then
        bcftools view -R ${family}-ensemble.db.txt.positions $fprefix.vcf.gz | bcftools sort | vt decompose -s - | vt uniq - -o $fprefix.subset.vcf.gz
        tabix $fprefix.subset.vcf.gz
        f_fix_sample_names $fprefix
        vcf.platypus.getNV.gatk3.sh $fprefix.subset.vcf.gz $reference
    fi

    cd ..

    # using Rscript from bcbio
    if [ "$type" == "wgs" ] || [ "$type" == "rnaseq" ]
    then
        noncoding="noncoding"
    else
        noncoding=""
    fi

    echo GENERATING REPORT WITH TYPE: "${type}"
    Rscript ~/cre/cre.vcf2db.R $family "${type}"
    
    cd $family
    #rm $family.create_report.csv $family.merge_reports.csv
    #for vcaller in {freebayes,gatk-haplotype,platypus}
    #do
#	rm ${family}-${vcaller}-annotated-decomposed.subset.vcf.gz ${family}-${vcaller}-annotated-decomposed.subset.vcf.gz.tbi
#   done
    cd ..
}

if [ -z $family ]
then
    family=$1
fi

echo $family

if [ -z "$rnaseq" ]
then
    rnaseq=0
fi

export depth_threshold=10

if [ "$type" == "rnaseq" ]
then
    export depth_threshold=5
    export severity_filter=ALL
fi

# set wes.regular as default type
if [ -z $type ]
then
    type="wes.regular"
fi

#no cleanup by default
if [ -z $cleanup ]
then
    cleanup=0
fi

if [ $cleanup -eq 1 ]
then
    f_cleanup
fi

if [ -z $max_af ]
then
    max_af=0.01
fi
export max_af

#make report by default
if [ -z $make_report ] || [ $make_report -eq 1 ]
then
    f_make_report
fi

#if cleanup set, also make the synonymous report
if [ $cleanup -eq 1 && "$type" == "wes.regular" ]
then
    type="wes.synonymous"
    f_make_report
fi
