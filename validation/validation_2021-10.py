import argparse
import csv
from datetime import date
import pandas as pd
import sys
from pathlib import Path


def db_output_to_dict(db_output):
    db1 = {}
    db2 = {}
    columns = [
        "impact_severity",
        "clinvar_sig",
        "clinvar_pathogenic",
        "callers",
        "clinvar_status",
        "gnomad_af_popmax",
    ]
    with open(db_output) as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")
        for row in reader:
            alt_depths = [
                int(row[col].strip(" ")) if not row[col] is None else row[col]
                for col in row
                if "gt_alt_depths" in col
            ]
            depths = [
                int(row[col].strip(" ")) if not row[col] is None else row[col]
                for col in row
                if "gt_depths" in col
            ]
            if row["db"] == args.prefix1:
                variant = (
                    row["chrom"]
                    + ":"
                    + row["start"]
                    + ":"
                    + row["end"]
                    + ":"
                    + row["ref"]
                    + ":"
                    + row["alt"]
                )
                db1[variant] = {
                    column: row[column] for column in columns if column in row
                }
                db1[variant].update({"alt_depths": alt_depths, "depths": depths})
            elif row["db"] == args.prefix2:
                variant = (
                    row["chrom"]
                    + ":"
                    + row["start"]
                    + ":"
                    + row["end"]
                    + ":"
                    + row["ref"]
                    + ":"
                    + row["alt"]
                )
                db2[variant] = {
                    column: row[column] for column in columns if column in row
                }
                db2[variant].update({"alt_depths": alt_depths, "depths": depths})
    return db1, db2


def get_explanations(report1_var, report2_var, report2_dir):
    # load tables of alt allele counts for each caller; parse AD if necessary; add max AD of samples
    tables = load_tables(report2_dir)
    gatk, freebayes, platypus = tables["gatk"], tables["freebayes"], tables["platypus"]
    # gatk alt depth cols
    if gatk is not None:
        gatk_ad_cols = get_gatk_ad_cols(gatk.columns)
        for col in gatk_ad_cols:
            gatk[col] = gatk.apply(lambda row: parse_ad(row[col]), axis=1)
        gatk["AD_max"] = gatk.apply(
            lambda row: get_max_ad([row[col] for col in gatk_ad_cols]), axis=1
        )
    # freebayes alt depth cols
    if freebayes is not None:
        freebayes_ad_cols = get_freebayes_ao_cols(freebayes.columns)
        for col in freebayes_ad_cols:
            freebayes[col] = freebayes.apply(lambda row: parse_ad(row[col]), axis=1)
        freebayes["AD_max"] = freebayes.apply(
            lambda row: get_max_ad([row[col] for col in freebayes_ad_cols]), axis=1
        )
    # platypus  alt depth cols
    if platypus is not None:
        platypus_ad_cols = get_platypus_nv_cols(platypus.columns)
        for col in platypus_ad_cols:
            platypus[col] = platypus.apply(lambda row: parse_ad(row[col]), axis=1)
        platypus["AD_max"] = platypus.apply(
            lambda row: get_max_ad([row[col] for col in platypus_ad_cols]), axis=1
        )

    tables = [table for table in [gatk, freebayes, platypus] if table is not None]

    explanation = {}
    for variant in report1_var:
        report1_var[variant]["alt_depths"] = max(report1_var[variant]["alt_depths"])
        report1_var[variant]["depths"] = min(report1_var[variant]["depths"])
        impact_severity = report1_var[variant]["impact_severity"]
        gnomad = float(report1_var[variant]["gnomad_af_popmax"])

        if variant not in report2_var:
            explanation[variant] = "Variant not present in comparison database"
        elif report2_var[variant]["callers"] in ["samtools", "freebayes", "platypus"]:
            explanation[
                variant
            ] = "Variant only called by one non-GATK caller in comparison database"
        # change in gnomAD AF:
        elif (float(report1_var[variant]["gnomad_af_popmax"]) < 0.01) & (
            float(report2_var[variant]["gnomad_af_popmax"]) >= 0.01
        ):
            explanation[
                variant
            ] = "GnomAD AF popmax greater than 0.01 in the other report"
        # change in impact severity (vep)
        elif (
            report1_var[variant]["impact_severity"] != "LOW"
            and report2_var[variant]["impact_severity"] == "LOW"
            and gnomad < 0.01
        ):
            impact_severity = report1_var[variant]["impact_severity"]
            explanation[
                variant
            ] = f"Change in impact_severity from {impact_severity} to LOW and gnomad_af_popmax < 0.01"
        # change in alt depth
        elif (
            report1_var[variant]["alt_depths"] >= 3
            and max(report2_var[variant]["alt_depths"]) < 3
        ):
            explanation[variant] = "Alt depth less than 3 in other report"
        elif (
            report1_var[variant]["alt_depths"] == -1
            and max(report2_var[variant]["alt_depths"]) < 3
        ):
            explanation[variant] = "Alt depth less than 3 in other report"
        # change in clinvar annotation, or alt depth for variant < 3 in tables
        elif gnomad < 0.01:
            if (
                report1_var[variant]["clinvar_pathogenic"] != "None"
                and report2_var[variant]["clinvar_pathogenic"] == "None"
            ):
                clin_path = report1_var[variant]["clinvar_pathogenic"]
                explanation[
                    variant
                ] = f"Change in clinvar_pathogenic from {clin_path} to None for variant with gnomad_af_popmax < 0.01 and impact_severity {impact_severity}"
            else:
                max_table_ad = get_variant_ad(variant, tables)
                if not max_table_ad or max_table_ad < 3:
                    explanation[
                        variant
                    ] = "Alt depth less than 3 in other report tables"
                else:
                    explanation[variant] = "Cannot explain"
        elif gnomad >= 0.01:
            if (
                report1_var[variant]["clinvar_pathogenic"] != "None"
                and report2_var[variant]["clinvar_pathogenic"] == "None"
            ):
                clin_path = report1_var[variant]["clinvar_pathogenic"]
                explanation[
                    variant
                ] = f"Change in clinvar_pathogenic from {clin_path} to None for variant with gnomad_af_popmax >= 0.01 and impact_severity {impact_severity}"
            elif report1_var[variant]["clinvar_pathogenic"] in [
                "Pathogenic",
                "Likely_pathogenic",
                "Conflicting_interpretations_of_pathogenicity",
            ] and report2_var[variant]["clinvar_pathogenic"] not in [
                "Pathogenic",
                "Likely_pathogenic",
                "Conflicting_interpretations_of_pathogenicity",
            ]:
                clin_path_1 = report1_var[variant]["clinvar_pathogenic"]
                clin_path_2 = report2_var[variant]["clinvar_pathogenic"]
                explanation[
                    variant
                ] = f"Change in clinvar_pathogenic from {clin_path_1} to {clin_path_2} for variant with gnomad_af_popmax >= 0.01 and impact_severity {impact_severity}"
            else:
                max_table_ad = get_variant_ad(variant, tables)
                if not max_table_ad or max_table_ad < 3:
                    explanation[
                        variant
                    ] = "Alt depth less than 3 in other report tables"
                else:
                    explanation[variant] = "Cannot explain"
        else:
            explanation[variant] = "Cannot explain"

    return explanation


def summarize_explanations(explanation_df, prefix):
    explanation_summary = explanation_df.groupby(["Explanation"]).count().reset_index()
    explanation_summary.to_csv(
        f"validation_summary_unique_in_{prefix}_{today}_counts.csv",
        header=["Explanation", "Count"],
        index=False,
    )


def load_tables(report_dir):
    # these are tables generated by cre that contain the alt depth count for the respective report variant
    # this is used by cre to populate alt depth counts for variants called by samtools and another caller (e.g. platypus)
    # because often the alt depth count for these variants is '-1' in the gemini db as samtools doesn't always have an alt depth field
    tables = list(Path(report_dir).glob("*table"))
    callers = ("gatk", "freebayes", "platypus")
    table_dict = dict.fromkeys(callers)
    for table in tables:
        filename = table.name
        if "gatk" in filename:
            gatk = pd.read_csv(table, sep="\t")
            table_dict["gatk"] = gatk
        elif "freebayes" in filename:
            freebayes = pd.read_csv(table, sep="\t")
            table_dict["freebayes"] = freebayes
        elif "platypus" in filename:
            platypus = pd.read_csv(table, sep="\t")
            table_dict["platypus"] = platypus
    return table_dict


def get_max_ad(ad_list):
    # get max alt depth of all samples in above tables
    ad_max = max(ad_list)
    return ad_max


def parse_ad(alt_obs):
    if pd.isna(alt_obs):
        alt_obs = 0
    # comma separated AO, AD values seem to be an artefact of decomposition; take max
    elif "," in alt_obs:
        alt_obs = max([int(ao) for ao in alt_obs.split(",")])
    else:
        alt_obs = int(alt_obs)

    return alt_obs


def get_freebayes_ao_cols(colnames):
    # AO: number of reads supporting alt allele
    ao_cols = []
    for col in colnames:
        if "AO" in col:
            ao_cols.append(col)
    return ao_cols


def get_platypus_nv_cols(colnames):
    # NV: number of reads supporting alt allele
    nv_cols = []
    for col in colnames:
        if "NV" in col:
            nv_cols.append(col)
    return nv_cols


def get_gatk_ad_cols(colnames):
    ad_cols = []
    for col in colnames:
        if "AD" in col:
            ad_cols.append(col)
    return ad_cols


def get_variant_ad(variant, tables):
    # get max alt depth from gatk, freebayes, platypus tables
    chr, start, end, ref, alt = variant.split(":")
    ad_list = []
    for table in tables:
        ad = table[
            (table["CHROM"] == chr)
            & (table["POS"] == int(end))
            & (table["REF"] == ref)
            & (table["ALT"] == alt)
        ]["AD_max"].values
        if len(ad) == 0:
            pass
        else:
            ad_list.append(int(ad[0]))
    if len(ad_list) == 0:
        return None
    else:
        return max(ad_list)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Determine cause for inclusion/exclusion of variants between two reports"
    )
    parser.add_argument(
        "-db_output1",
        help="output in format <prefix>.uniq.db.txt from gemini_compare.sh for first report",
        required=True,
    )
    parser.add_argument(
        "-db_output2",
        help="output in format <prefix>.uniq.db.txt from gemini_compare.sh for second report",
        required=True,
    )
    parser.add_argument(
        "-prefix1",
        help="prefix of first report, e.g. 428.wes.regular.2020-10-14",
        required=True,
    )
    parser.add_argument(
        "-prefix2",
        help="prefix of second report, e.g. 428.wes.regular.2020-10-19",
        required=True,
    )
    parser.add_argument(
        "-dir1",
        help="path to report directory for first report",
        required=True,
    )
    parser.add_argument(
        "-dir2",
        help="path to report directory for second report",
        required=True,
    )
    args = parser.parse_args()

    # For variants unique to report 1, determine reason they were not included in report 2
    db1_unique = db_output_to_dict(args.db_output1)
    report1_var = db1_unique[0]
    report2_var = db1_unique[1]
    explanation_1 = get_explanations(report1_var, report2_var, args.dir2)

    today = date.today()
    today = today.strftime("%Y-%m-%d")
    explanation_1_df = pd.DataFrame.from_dict(
        explanation_1, orient="index"
    ).reset_index()
    if len(explanation_1_df) != 0:
        explanation_1_df.columns = ["Variant", "Explanation"]
        explanation_1_df.to_csv(
            f"validation_summary_unique_in_{args.prefix1}_{today}.csv",
            index=False,
        )
        summarize_explanations(explanation_1_df, args.prefix1)
    # For variants unique to report 2, determine reason they were not included in report 1
    db2_unique = db_output_to_dict(args.db_output2)
    report1_var = db2_unique[0]
    report2_var = db2_unique[1]
    explanation_2 = get_explanations(report2_var, report1_var, args.dir1)

    explanation_2_df = pd.DataFrame.from_dict(
        explanation_2, orient="index"
    ).reset_index()
    if len(explanation_2_df) != 0:
        explanation_2_df.columns = ["Variant", "Explanation"]
        explanation_2_df.to_csv(
            f"validation_summary_unique_in_{args.prefix2}_{today}.csv",
            index=False,
        )
        summarize_explanations(explanation_2_df, args.prefix2)
