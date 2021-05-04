#!/bin/bash
# generate both the regular and synonymous reports for a family
# pass in the family id as the first positional parameter, use -email to recieve an email when both jobs are done
# usage: generate_reports.sh <family id> <report_type> [optional -email] 

family=$1
report_type=$2
ped=$3
curr_date=$(date +"%Y-%m-%d")

rerun_folder="${family}_${curr_date}"


if [ "$report_type" == "wes" ] || [ "$report_type" == "wes.both" ]; then
	eval script="~/cre/cre.vcf2cre.sh"
elif [ "$report_type" == "wgs" ] || [ "$report_type" == "denovo" ]; then
	eval script="~/crg/crg.vcf2cre.sh"
else
	echo "Please enter a report type (wes, wes.both, or wgs)"
	exit
fi

if [ -d "$family" ]; then
	cd $family
	family_vcf="${family}-ensemble-annotated-decomposed.vcf.gz"
else
	echo "Family folder ${family} not found"
	exit
fi

if [ -f $family_vcf ]; then
	mkdir $rerun_folder
	cd $rerun_folder
	# could link instead of copying but not sure whether there
	# would be side effects on the linked file (#TODO: test)
	cp ../${family_vcf} .
	# for exome reports, need variant tables from the four variant callers
	cp ../*table . 
	if [ ! z "$ped" ]
		then
			vcf2cre_job="$(qsub "${script}" -v original_vcf="${family_vcf}",project=${family},ped=${ped})"
		else
			vcf2cre_job="$(qsub "${script}" -v original_vcf="${family_vcf}",project=${family})"
	fi
else
	echo "${family_vcf} not present, exiting."
	cd ../..
	exit
fi

# generate wes.regular and wes.synonymous reports only for exomes
if [ "$report_type" = "wes.both" ]; then
	first_report_job="$(qsub ~/cre/cre.sh -W depend=afterany:"${vcf2cre_job}" -v family=${family},type=wes.regular)"
	echo "Regular WES Report Job ID: ${first_report_job}"
	report_job="$(qsub ~/cre/cre.sh -W depend=afterany:"${first_report_job}" -v family=${family},type=wes.synonymous)"
	echo "Synonymous WES Report Job ID: ${report_job}"
# generate wes.regular report for genomes
elif [ "$report_type" = "wes" ]; then
	report_job="$(qsub ~/cre/cre.sh -W depend=afterany:"${vcf2cre_job}" -v family=${family},type=wes.regular)"
	echo "Regular WES Report Job ID: ${report_job}"
# generate wgs report for genomes (i.e. panel, panel-flank)
elif [ "$report_type" = "wgs" ]; then
	report_job="$(qsub ~/cre/cre.sh -W depend=afterany:"${vcf2cre_job}" -v family=${family},type=wgs)"
  echo "WGS Report Job ID: ${report_job}"
# generate denovo report
elif [ "$report_type" = "denovo" ]; then
	report_job="$(qsub ~/cre/cre.sh -W depend=afterany:"${vcf2cre_job}" -v family=${family},type=denovo)"
  echo "WGS Report Job ID: ${report_job}"
fi

echo "The re-run subfolder will be cleaned up after the reports are created"
cleanup_job="$(qsub ~/cre/cleanup_run.sh -W depend=afterok:"${report_job}" -v family=${family})"

cd ../..
