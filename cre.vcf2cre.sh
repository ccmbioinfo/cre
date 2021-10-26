#!/bin/bash
#PBS -l walltime=100:00:00,nodes=1:ppn=14
#PBS -joe .
#PBS -d .
#PBS -l vmem=70g,mem=70g

# 50g is crucial -20,30 crashes sometimes

# annotates vcf with VEP and vcfanno for cre for cre report generation
# parameters:
# original_vcf = file.vcf.gz, no tabix index is needed in the dir
# project = case = family = S11 (example)
# [optional] ped = file.ped

# qsub ~/cre/cre.vcf2cre.sh -v original_vcf=file.vcf,project=412

# in old bcbio vcf files from rna-seq pipeline, i.e. combined-annotated-rnaedit.vcf.gz, bcbio-nextgen1.0.4, some info field formats are wrong:
# AC, AF, MLEAC, MLEAF Number=1 not A. Because of that vt is not properly decomposing multiallelic variants, and vcf2db can't create gemini database
# solution is to put Number=A in the vcf header or rerun variant calling with the latest bcbio

# if input is from TCAG (HAS) it does not have DP INFO field, we need to fake it from FORMAT DP for SNVs and from DPI for indels:
# cre.vcf.has2dp.sh
# gunzip -c 331606_S1.flt.nochr.vcf.gz | grep "^#" > 331606.vcf
# add to header: 
##INFO=<ID=DP,Number=1,Type=Integer,Description="Approximate read depth; some reads may have been filtered">
# gunzip -c 331606_S1.flt.nochr.vcf.gz | grep -v "^#"  | grep PASS | sed s/":DPI:"/":DP:"awk -F ':' '{print $0"\tDP="$9}' | awk -F "\t" '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$11";"$8"\t"$9"\t"$10}' >> 331606.vcf

bname=`basename $original_vcf .vcf.gz`

echo "###############################################"
echo `date` ": Removing annotation of " $original_vcf
echo "###############################################"

cre.annotation.strip.sh $original_vcf

gunzip -c $bname.no_anno.vcf.gz | grep "^#"  > $project.vcf
gunzip -c $bname.no_anno.vcf.gz | grep -v "^#" | grep PASS | grep -v possible_rnaedit | egrep -v "^GL000" >> $project.vcf
bgzip $project.vcf
tabix $project.vcf.gz

bcftools sort -Oz $project.vcf.gz > $project.sorted.vcf.gz
tabix $project.sorted.vcf.gz

cre.vt.decompose.sh $project.sorted.vcf.gz

echo "################################################"
echo `date` "Annotating with VEP..."
echo "################################################"

cre.vep.sh $project.sorted.decomposed.vcf.gz

echo "################################################"
echo `date` "Annotating with vcfanno ..."
echo "################################################"
cre.vcfanno.sh $project.sorted.decomposed.vepeffects.vcf.gz

if [ -z $ped ]
then
    echo "no ped file, generating ..."
    bcftools query -l $original_vcf > samples.txt
    > $project.ped
    for sample in `cat samples.txt`
    do
    	echo -e "1\t"$sample"\t0\t0\t0\t0" >> $project.ped
    done
    ped=$project.ped
fi

echo "#################################################"
echo `date` "Generating gemini database"
echo "#################################################"

mkdir $project
vcf2db.py $project.sorted.decomposed.vepeffects.annotated.vcf.gz 	$ped	${project}/${project}-ensemble.db

mv $project.sorted.decomposed.vepeffects.annotated.vcf.gz ${project}/${project}-ensemble-annotated-decomposed.vcf.gz
mv $project.sorted.decomposed.vepeffects.annotated.vcf.gz.tbi ${project}/${project}-ensemble-annotated-decomposed.vcf.gz.tbi

cd $project
ln -s ${project}-ensemble-annotated-decomposed.vcf.gz ${project}-gatk-haplotype-annotated-decomposed.vcf.gz
ln -s ${project}-ensemble-annotated-decomposed.vcf.gz.tbi ${project}-gatk-haplotype-annotated-decomposed.vcf.gz.tbi

cd ..

rm $bname.no_anno.vcf.gz
rm $bname.no_anno.vcf.gz.tbi

rm $project.vcf.gz
rm $project.vcf.gz.tbi

rm $project.sorted*.vcf.gz
rm $project.sorted*.vcf.gz.tbi

echo "#####################################################"
echo `date` " DONE"
echo "#####################################################"
