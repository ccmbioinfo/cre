#!/usr/bin/env python3

import pandas as pd
import numpy as np
from glob import glob
from os.path import basename, join
from datetime import datetime

DAY = datetime.now().strftime("%Y-%m-%d")

VARIANT_WISE_REPORT = "%s.c4r.variant_wise.csv" % (DAY)
SAMPLE_WISE_REPORT = "%s.c4r.sample_wise.csv" % (DAY)
SEEN_IN_C4R_COUNTS = "seen_in_c4r_counts.txt"
SEEN_IN_C4R_SAMPLES = "seen_in_c4r_samples.txt"

ALL_SAMPLES_DFS = []
RESULTS_PATH = "/srv/shared/hpf/exomes/results/"
ANALYSIS_PATH = "/srv/shared/analyses/exomes/"
IGNORE_FOLDERS = set(join(RESULTS_PATH, folder) for folder in ["calx/", "misc/", "run_statistics/", "database/"])
IGNORE_FOLDERS_ANALYSIS = set(join(ANALYSIS_PATH, folder) for folder in ["bam_to_cram", "HPO", "GSO", "scripts"])
WES_REPORT_FIELDS = ['Position', 'Ref', 'Alt', 'Variation', 'Refseq_change', 'Depth', 'Gene', 'Conserved_in_20_mammals', \
'Sift_score', 'Polyphen_score', 'Cadd_score', 'Gnomad_af']

def get_latest_wes_report_paths():
	report_paths = []

	for family_prefix_dir in glob(join(RESULTS_PATH, "*x/")):
		if family_prefix_dir in IGNORE_FOLDERS:
			continue

		# sorted - ascending order
		for family in glob(join(family_prefix_dir, "*/")):
			# for reports generated with bcbio: 
			reports_bcbio = sorted([report for report in glob(join(family, "*wes*")) if "clinical" not in report \
			and "synonymous" not in report])
			# for reports generated with crg:
			report_dir = family +  f"/*/report/coding/*/*wes*"
			reports_crg2 = sorted([report for report in glob(report_dir) if "clinical" not in report \
			and "synonymous" not in report])
			reports = reports_bcbio + reports_crg2


			if len(reports) >= 1:
				report_paths.append(reports[-1]) #latest report since prefix is: yyyy-mm-dd
			else:
				print("No report files found for %s" % family)

	return report_paths

def get_report_paths_analysis_path():
	report_paths = []

	for family in glob(join(ANALYSIS_PATH, "*/")):
		if family in IGNORE_FOLDERS_ANALYSIS:
			continue

		# sorted - ascending order
		report_dir = family +  f"/*/report/coding/*/*wes*"
		reports = sorted([report for report in glob(report_dir) if "clinical" not in report \
		and "synonymous" not in report])

		if len(reports) >= 1:
			report_paths.append(reports[-1]) #latest report since prefix is: yyyy-mm-dd
		else:
			print("No report files found for %s" % family)
	
	return report_paths

paths = get_latest_wes_report_paths() + get_report_paths_analysis_path()

for report in paths:
	print("Parsing %s" % report)

	sep = "," if report.endswith(".csv") else "\t"

	try:
		df = pd.read_csv(report, sep=sep)
	except UnicodeDecodeError:
		print("UnicodeDecodeError on %s. Trying latin-1 decoding." % report)
		df = pd.read_csv(report, encoding='latin-1', sep=sep)

	zygosity_cols = [col for col in df.columns if col.startswith("Zygosity.")]
	burden_cols = [col for col in df.columns if col.startswith("Burden.")]
	alt_depths_cols = [col for col in df.columns if col.startswith("Alt_depths.")]

	samples = [col.replace("Zygosity.", "").strip() for col in zygosity_cols]

	df = df[[ col for col in df.columns if col in WES_REPORT_FIELDS \
	or col in zygosity_cols \
	or col in burden_cols \
	or col in alt_depths_cols ]]

	for sample in samples:

		sample_zygosity_col = "%s%s" % ('Zygosity.', sample)
		sample_burden_col = "%s%s" % ('Burden.', sample)
		sample_alt_depths_col = "%s%s" % ('Alt_depths.', sample)

		missing_report_fields = set(WES_REPORT_FIELDS) - set(df.columns)
		
		sample_df = df[ [field for field in WES_REPORT_FIELDS if field not in missing_report_fields] ].copy()

		for field in missing_report_fields:
			# first pandas function ignores nan
			sample_df[field] = np.nan

		for field in ["Zygosity", "Burden", "Alt_depths"]:
			sample_df[field] = np.nan

		if sample_zygosity_col in df.columns:
			sample_df["Zygosity"] = df[sample_zygosity_col]

		if sample_burden_col in df.columns:
			sample_df["Burden"] = df[sample_burden_col]

		if sample_alt_depths_col in df.columns:
			sample_df["Alt_depths"] = df[sample_alt_depths_col]

		# filter out homozygous reference variants or insufficent coverage variants from frequency
		sample_df = sample_df[ (sample_df["Zygosity"] != '-') & \
		(sample_df["Zygosity"] != 'Insufficient_coverage') & \
		(sample_df["Zygosity"] != 'Insufficient coverage')]

		sample_df["Sample"] = sample
		ALL_SAMPLES_DFS.append(sample_df)

df = pd.concat(ALL_SAMPLES_DFS).astype(str)

df.set_index(['Position', 'Ref', 'Alt']).to_csv(SAMPLE_WISE_REPORT)

df = df.groupby(['Position', 'Ref', 'Alt']).agg({
	'Variation' : 'first',
	'Refseq_change' : 'first',
	'Depth' : list,
	'Gene' : 'first',
	'Conserved_in_20_mammals' : 'first',
	'Sift_score' : 'first',
	'Polyphen_score' : 'first',
	'Cadd_score' : 'first',
	'Gnomad_af' : 'first',
	'Zygosity' : list,
	'Burden' : list,
	'Alt_depths' : list,
	'Sample' : list })

df['Frequency'] = df['Sample'].apply(lambda x: len(set(x))) # there could be multiple reports assoc with one family if there is a new analysis in the analysis dir
df['Zygosity'] = df['Zygosity'].apply(lambda g: '; '.join(g))
df['Burden'] = df['Burden'].apply(lambda g: '; '.join(g))
df['Alt_depths'] = df['Alt_depths'].apply(lambda g: '; '.join(g))
df['Sample'] = df['Sample'].apply(lambda g: '; '.join(g))
df['Depth'] = df['Depth'].apply(lambda g: '; '.join(g))

df.rename(columns={'Sample' : 'Samples'}, inplace=True)

df.to_csv(VARIANT_WISE_REPORT)

df = df.reset_index()
df['Position-Ref-Alt'] = df['Position'].str.cat(df[['Ref', 'Alt']], sep='-')
df[['Position-Ref-Alt', 'Frequency']].to_csv(SEEN_IN_C4R_COUNTS, sep='\t', index=False)
df[['Position-Ref-Alt', 'Samples']].to_csv(SEEN_IN_C4R_SAMPLES, sep='\t', index=False)

print("%s Database generation finished successfully." % DAY)
