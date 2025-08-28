# Updates OMIM annotations for the cre pipeline
# Use: navigate to omim directory and run ./update_omim.sh
# Input:
#  * mimTitles.txt genemap2.txt mim2gene.txt
# OMIM last updated: 07-22-2020. place the above files into cre/data and run this script within cre/ 
# Output: 
#  * hgnc_$DATE.tsv
#    * a table used to map gene names with ensembl gene id's
#    * required file for cre.sh
#  * `prefix_hgnc_join_omim_phenos_yyyy-mm-dd.tsv`
#    * a table mapping each OMIM gene to all of its known HGNC gene symbol/previous symbols/synonyms
#    * allows for simple left-join to 'variant' column of CRE reports
#  * `prefix_remaining_unmapped_omim_with_info_yyyy-mm-dd.tsv`
#    * a list of the the unmapped OMIM genes

DATE=`date +%Y-%m-%d`
HGNC="hgnc_${DATE}.txt"
OMIM_PREFIX="OMIM"
CRE_DIR=$1

# get latest mappings from hgnc
# http params include default checkboxes, gene aliases, with ensembl gene id from hgnc, and ncbi gene id
wget \
"https://www.genenames.org/cgi-bin/download/custom?col=gd_hgnc_id&col=gd_app_sym&col=gd_app_name&col=gd_status&col=gd_prev_sym&col=gd_aliases&col=gd_pub_chrom_map&col=gd_pub_acc_ids&col=gd_pub_refseq_ids&col=gd_pub_eg_id&col=gd_pub_ensembl_id&status=Approved&status=Entry%20Withdrawn&hgnc_dbtag=on&order_by=gd_app_sym_sort&format=text&submit=submit"  \
-O ${CRE_DIR}/data/${HGNC}

Rscript ${CRE_DIR}/map_omim_hgnc.R ${CRE_DIR}/data/${HGNC} ${CRE_DIR}/data/mimTitles.txt ${CRE_DIR}/data/genemap2.txt ${CRE_DIR}/data/mim2gene.txt ${OMIM_PREFIX}

mv ${OMIM_PREFIX}_hgnc_join_omim_phenos_${DATE}.tsv ${CRE_DIR}/data/
mv ${OMIM_PREFIX}_remaining_unmapped_omim_with_info_${DATE}.tsv ${CRE_DIR}/data/

echo "Done!"
