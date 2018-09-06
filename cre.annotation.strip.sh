#!/bin/bash

bname=`basename $1 .vcf.gz`
vt rminfo -t CSQ,af_adj_exac_afr,af_adj_exac_amr,af_adj_exac_eas,af_adj_exac_fin,af_adj_exac_nfe,af_adj_exac_oth,af_adj_exac_sas,af_exac_all,common_pathogenic,max_aaf_all,num_exac_Het,num_exac_Hom,rs_ids,af_1kg_amr,af_1kg_eas,af_1kg_sas,af_1kg_afr,af_1kg_eur,af_1kg_all,fitcons,encode_consensus_gm12878,encode_consensus_h1hesc,encode_consensus_helas3,encode_consensus_hepg2,encode_consensus_huvec,encode_consensus_k562,dgv,hapmap1,hapmap2,gnomAD_AF,gnomAD_AFR_AF,gnomAD_AMR_AF,gnomAD_ASJ_AF,gnomAD_EAS_AF,gnomAD_FIN_AF,gnomAD_NFE_AF,gnomAD_OTH_AF \
	     $1 -o $bname.no_anno.vcf.gz
tabix $bname.no_anno.vcf.gz