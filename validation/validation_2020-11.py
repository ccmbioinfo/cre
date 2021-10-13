import argparse
import csv
from datetime import date
import pandas as pd
import sys 

def db_output_to_dict(db_output):
    db1 = {}
    db2 = {}
    columns = ["impact_severity", "clinvar_sig", "clinvar_pathogenic", "callers", "clinvar_status", "gnomad_af_popmax"]
    with open(db_output) as csvfile:
        reader = csv.DictReader(csvfile, delimiter='\t')
        for row in reader:
            alt_depths = [ int(row[col].strip(' ')) if not row[col] is None else row[col] for col in row if 'gt_alt_depths' in col ]
            depths = [ int(row[col].strip(' ')) if not row[col] is None else row[col] for col in row if 'gt_depths' in col ]
            if row['db'] == args.prefix1:
                variant = row['chrom'] + ':' +  row['start'] + ':' + row['end'] + ':' + row['ref'] + ':' + row['alt']
                db1[variant] = {column: row[column] for column in columns if column in row }
                db1[variant].update({'alt_depths': alt_depths, 'depths': depths})
            elif row['db'] == args.prefix2:
                variant = row['chrom'] + ':' +  row['start'] + ':' + row['end'] + ':' + row['ref'] + ':' + row['alt']
                db2[variant] = {column: row[column] for column in columns if column in row }
                db2[variant].update({'alt_depths': alt_depths, 'depths': depths})
    return db1, db2

def get_explanations(report1_var, report2_var):
    explanation = {}
    for variant in report1_var:
        report1_var[variant]['alt_depths'] = max(report1_var[variant]['alt_depths'])
        report1_var[variant]['depths'] = min(report1_var[variant]['depths'])

        if variant not in report2_var:
            if report1_var[variant]['callers'] == 'gatk-haplotype':
                explanation[variant] = 'Variant only called by GATK'
            else:   
                explanation[variant] = 'Variant not present in comparison database'
        # if only one caller and not GATK:
        elif report2_var[variant]['callers'] == 'freebayes' or report2_var[variant]['callers'] == 'samtools' or report2_var[variant]['callers'] == 'platypus':
            explanation[variant] = 'Variant only called by one of freebayes, samtools, or platypus in other report'
        # variants with non-null clinvar status and gnomad_af_popmax >1% are included in new reports
        elif not report1_var[variant]['clinvar_status'] in ['None',"0"] and report2_var[variant]['clinvar_status'] == "0" and float(report1_var[variant]['gnomad_af_popmax']) > 0.01:
            explanation[variant] = 'clinvar_status: %s'%(report1_var[variant]['clinvar_status'])
        # change in clinvar annotation
        elif report1_var[variant]['clinvar_sig'] != 'None' and report2_var[variant]['clinvar_sig'] == 'None':
            explanation[variant] = 'Change in clinvar_sig from %s to None'%report1_var[variant]['clinvar_sig']
        elif report1_var[variant]['clinvar_pathogenic'] != 'None' and report2_var[variant]['clinvar_pathogenic'] == 'None':
            explanation[variant] = 'Change in clinvar_pathogenic from %s to None'%report1_var[variant]['clinvar_pathogenic']
        # change in impact severity (vep)
        elif report1_var[variant]['impact_severity'] != 'LOW' and report2_var[variant]['impact_severity'] == 'LOW':
            explanation[variant] = 'Change in impact_severity from %s to LOW'%report1_var[variant]['impact_severity']

        #max(alt_depth)>= 3 (why variants are included)
        elif report1_var[variant]['alt_depths'] >= 3 and max(report2_var[variant]['alt_depths']) < 3 :
            explanation[variant] = 'Alt depth less than 3 in other report'
        elif report1_var[variant]['alt_depths'] == -1 and max(report2_var[variant]['alt_depths']) < 3 :
            explanation[variant] = 'Alt depth less than 3 in other report'
        #were very old reports not filtered by depth?
        elif 0 < report1_var[variant]['alt_depths'] < 3 and max(report2_var[variant]['alt_depths']) < 3 :
            explanation[variant] = 'Alt depth less than 3 in both'
        #https://github.com/ccmbioinfo/cre/blob/4c7ebe13f8a36251a0ca7f4100a424b9598c1f84/cre.gemini2txt.vcf2db.sh#L110
        #impact or callers affect inclusion/exclusion
        elif (report1_var[variant]['alt_depths'] == -1 or max(report2_var[variant]['alt_depths']) == -1 ) and ( not "gatk" in report1_var[variant]['callers']):
            explanation[variant] = 'Alt depths -1 and called by non-GATK callers'
        #very old report may have depth > 10 filter
        elif report1_var[variant]['alt_depths'] >= 3 and max(report2_var[variant]['depths']) < 10 :
             explanation[variant] = 'Max depth less than 10 in both reports but alt depth >3 in this report'
        #gnomad af popmax may have changed between reports
        elif float(report1_var[variant]['gnomad_af_popmax']) < 0.01 and float(report2_var[variant]['gnomad_af_popmax']) >= 0.01 :
            explanation[variant] = 'GnomAD AF popmax greater than 0.01 in the other report'
        else:
            explanation[variant] = 'Cannot explain'
    return explanation

def summarize_explanations(explanation_df, prefix):
    explanation_summary = explanation_df.groupby(['Explanation']).count().reset_index()
    explanation_summary.to_csv('validation_summary_unique_in_%s_%s_counts.csv'%(prefix,today), header=['Explanation', 'Count'], index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Determine cause for inclusion/exclusion of variants between two reports')
    parser.add_argument('-db_output1', help='output in format <prefix>.uniq.db.txt from gemini_compare.sh for first report', required=True)
    parser.add_argument('-db_output2', help='output in format <prefix>.uniq.db.txt from gemini_compare.sh for second report', required=True)
    parser.add_argument('-prefix1', help='prefix of first report, e.g. 428.wes.regular.2020-10-14', required=True)
    parser.add_argument('-prefix2', help='prefix of first report, e.g. 428.wes.regular.2020-10-19', required=True)
    args = parser.parse_args()

    # For variants unique to report 1, determine reason they were not included in report 2
    db1_unique = db_output_to_dict(args.db_output1)
    report1_var = db1_unique[0]
    report2_var = db1_unique[1]
    explanation_1 = get_explanations(report1_var,report2_var)
    
    today = date.today()
    today = today.strftime("%Y-%m-%d")
    explanation_1_df = pd.DataFrame.from_dict(explanation_1, orient='index').reset_index()
    if len(explanation_1_df) != 0:
        explanation_1_df.columns = ['Variant', 'Explanation']
        explanation_1_df.to_csv('validation_summary_unique_in_%s_%s.csv'%(args.prefix1,today), index=False)
        summarize_explanations(explanation_1_df, args.prefix1)
    # For variants unique to report 2, determine reason they were not included in report 1
    db2_unique = db_output_to_dict(args.db_output2)
    report1_var = db2_unique[0]
    report2_var = db2_unique[1]
    explanation_2 = get_explanations(report2_var, report1_var)
    
    explanation_2_df = pd.DataFrame.from_dict(explanation_2, orient='index').reset_index()
    if len(explanation_2_df) != 0:
        explanation_2_df.columns = ['Variant', 'Explanation']
        explanation_2_df.to_csv('validation_summary_unique_in_%s_%s.csv'%(args.prefix2,today), index=False)
        summarize_explanations(explanation_2_df, args.prefix2)
