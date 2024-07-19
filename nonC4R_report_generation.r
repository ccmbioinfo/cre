# variant report generator

#example usage: Rscript ~/cre/variant_list_annotation.r $family wgs /hpf/largeprojects/ccmbio/nhanafi/c4r/downloads/databases/

# store date to be used when writing files
#datetime <- format(Sys.time(),"%Y-%m-%d_%H-%M") # can use to get timestamp for testing
datetime <- format(Sys.time(),"%Y-%m-%d")
print(datetime)

# Rscript ~/cre/cre.vcf2.db.R <family> noncoding|default=NULL,coding <database path>
add_placeholder <- function(variants, column_name, placeholder){
    variants[,column_name] <- with(variants, placeholder)
    return(variants)
}

get_variants_from_file <- function (filename){
    variants <- read.delim(filename, stringsAsFactors = F)
    return(variants)
}

# deleted genotype2zygocity

# output : family.ensemble.txt
create_report <- function(family, samples, type){
    file <- paste0(family, ".variants.txt")
    variants <- get_variants_from_file(file)
    
    impact_file <- paste0(family, ".variant_impacts.txt")
    impacts <- get_variants_from_file(impact_file)
    

    
    #Column1 - Position
    variants$Position <- with(variants, paste(Chrom, Pos, sep = ':'))
    
    #Column2 - UCSC link
    sUCSC1 <- "=HYPERLINK(\"http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&hgt.out3=10x&position="
    sUCSC2 <- "\",\"UCSC_link\")"
    variants$UCSC_Link <- with(variants, paste(sUCSC1, Position, sUCSC2, sep = ''))

    # Column3 = GNOMAD_Link
    variants$GNOMAD_POS <- with(variants, paste(Chrom,Pos,Ref,Alt, sep='-'))
    sGNOMAD1 <- "=HYPERLINK(\"http://gnomad.broadinstitute.org/variant/"
    sGNOMAD2 <- "?dataset=gnomad_r3"
    sGNOMAD3 <- "\",\"GNOMAD_link\")"
    variants$GNOMAD_Link <- with(variants, paste(sGNOMAD1, GNOMAD_POS, sGNOMAD2, sGNOMAD3, sep = ''))

    # Columns 4,5: Ref,Alt

    # Column6 - Gene
    variants$Gene[variants$Gene == ""] <- NA
    
    # Column9 = gts
    
    # Column10 = Variation
    
    # Column11 =  Info
    variants <- add_placeholder(variants, "Info", "Info")
    
    for (i in 1:nrow(variants)){
        #debug: i=1  
        v_id <- variants[i,"Variant_id"]
        gene <- variants[i,"Gene"]
        #for WES reports we need only coding impacts in the info field, for WGS we need all
        if (coding){
            gene_impacts <- subset(impacts, variant_id == v_id & is_coding == 1,
                                   select = c("exon", "hgvsc", "hgvsp"))
        }else{
            gene_impacts <- subset(impacts, variant_id == v_id, 
                                   select = c("exon", "hgvsc", "hgvsp"))
        }
        
    gene_impacts$gene <- rep(gene, nrow(gene_impacts))
    
        gene_impacts$exon[gene_impacts$exon==''] <- 'NA'
        
        gene_impacts <- gene_impacts[c("gene", "exon", "hgvsc", "hgvsp")]
        
        if (nrow(gene_impacts) > 0){
            v_impacts <- paste0(gene_impacts$gene, ":exon", gene_impacts$exon,
                                ":", gene_impacts$hgvsc, ":", gene_impacts$hgvsp)
						# replace %3D with its url encoded character, "="
						v_impacts <- str_replace_all(v_impacts,"%3D","=")
            s_impacts <- paste(v_impacts, collapse = ",")
        }
        else s_impacts <- 'NA'
      
        variants[i,"Info"] <- s_impacts
    }
    
    # Column12 - Refseq_change
    variants <- add_placeholder(variants, "Refseq_change", "NA")
    
    # Columns 13,14 - Depth, Quality

    # Column 15 - Alt_depth - from v.gt_alt_depths - DELETED

    # Column 16 - Trio_coverage - fixed in merge_reports function - DELETED
    
    # Column17 = Ensembl_gene_id

    # Column18 = Gene_description
    gene_descriptions <- read.delim2(paste0(default_tables_path,"/ensembl_w_description.txt"), 
                                     stringsAsFactors = F)
    variants <- merge(variants, gene_descriptions, by.x = "Ensembl_gene_id",
                      by.y = "ensembl_gene_id", all.x = T)
    
    # Column19 - Omim_phenotype
    # Column20 - Omim_inheritance 
    omim_map_file <- paste0(default_tables_path,"/OMIM_hgnc_join_omim_phenos_2021-10-19.tsv")
    if(file.exists(omim_map_file)){
    # read in tsv
    hgnc_join_omim_phenos <- read.delim(omim_map_file, stringsAsFactors=FALSE)
    print(head(hgnc_join_omim_phenos))
    print("Just read the hgnc join omim phenos")
    # select only relevant columns from the key file (gene name (to be joined on), the mim inheritance, and the phenotypes)
    hgnc_omim <-  hgnc_join_omim_phenos %>%
        dplyr::select(gene_name, omim_phenotype, omim_inheritance) %>%
        mutate(gene_name = replace(gene_name, is.na(gene_name), ""))
    print("Successfully altered hgnc_omim")
    # assuming the column w/ gene name is 'gene_name'

    variants <-  left_join(variants, hgnc_omim, by = c("Gene" = "gene_name")) 
    }
    else{
    print(paste0("File not found: ",omim_map_file))
    }

    print("Successfully joined on OMIM file")

    # Column 21 = Orphanet
    # previous name - orphanet.deduplicated.txt
    orphanet_file_name <- paste0(default_tables_path,"/orphanet.txt")
       
    if (file.exists(orphanet_file_name)){
	    orphanet <- read.delim(orphanet_file_name, stringsAsFactors = F)  
	    variants <- merge(variants, orphanet, all.x = T)
    
	    variants$Orphanet[is.na(variants$Orphanet)] = 0
    }
    
    # Column 22 - Clinvar
    
    # Column 23 - Ensembl_transcript_id
    
    # Column 24 - AA_position
    # changing separator from / to _ because otherwise excel converts it into date
    variants[,"AA_position"] <- with(variants, gsub("/", "_", AA_position), fixed = T)
    
    # Column 25 - Exon
    variants[,"Exon"] <- with(variants,gsub("/", "_", Exon), fixed = T)
    
    # Column 26 - Protein_domains
    
    # Column 27, 28 = C4R_WES_counts, C4R_WES_samples - DELETED
    
    # Columns 29,30,31,32: HGMD
    for(hgmd_field in c("HGMD_id", "HGMD_gene", "HGMD_tag", "HGMD_ref")){
        variants <- add_placeholder(variants, hgmd_field, "NA")
    }

    # Column 33 - rsIds

    # population frequencies
    # Column34 = Gnomad_af
    # Column35 = Gnomad_af_popmax
    
    # Gnomad gene constraint scores
    # Column36 = Gnomad_oe_lof_score
    # Column37 = Gnomad_oe_mis_score
    gnomad_scores_file <- paste0(default_tables_path, "/gnomad_scores.csv")
    gnomad_scores <- read.csv(gnomad_scores_file, stringsAsFactors = F)
    variants <- merge(variants, gnomad_scores, all.x = T, all.y = F)

    # Column38 = Gnomad_ac
    # Column39 = Gnomad_hom
    for (field in c("Gnomad_ac","Gnomad_hom")){
        variants[,field] <- with(variants,gsub("-1", "0", get(field), fixed = T))
        variants[,field] <- with(variants,gsub("None", "0", get(field), fixed = T))
    }
    
    # Column41 - Conserved_in_30_mammals
    # Column 42? - SpliceAI (actually 47, these column indexes are no longer accurate)
    variants <- add_placeholder(variants, "SpliceAI_impact", "")
    for (i in 1:nrow(variants)){
        print(i)
        if (variants[i,"SpliceAI_score"] == ""){
            variants[i, "SpliceAI_impact"] <- "NA|NA|NA"
            variants[i, "SpliceAI_score"] <- 0
        } else {
            spliceai <- strsplit(variants[i,"SpliceAI_score"], ",", fixed = T)[[1]]
            score_list <- c("NA", "NA", 0, "NA")
            names(score_list) <- c("gene", "impact", "score", "pos")
            for (anno in spliceai){
                anno <- strsplit(anno, "|", fixed = T)[[1]]
                gene <- anno[2]
                DS_AG <- anno[3]
                DS_AL <- anno[4]	
                DS_DG <- anno[5]	
                DS_DL <- anno[6]	
                DP_AG <- anno[7]	
                DP_AL <- anno[8]	
                DP_DG <- anno[9]	
                DP_DL <- anno[10]	
                scores <- c(as.numeric(DS_AG), as.numeric(DS_AL), as.numeric(DS_DG), as.numeric(DS_DL))
                names(scores) <- c("acceptor_gain", "acceptor_loss", "donor_gain", "donor_loss")
                max_score <- max(scores)        
                for (name in names(scores)){
                    if (scores[name] == max_score){name_max_score <- name}
                }
                if (name_max_score == 0){
                    impact <- "NA"
                } else {
                    impact <- name_max_score
                }
                if (score_list["score"] < max_score){
                    score_list["score"] <- max_score
                    score_list["gene"] <- gene
                    score_list["impact"] <- impact
                    if (impact == "acceptor_gain"){
                        score_list["pos"] <- DP_AG
                        } else if (impact == "acceptor_loss"){
                        score_list["pos"] <- DP_AL
                        } else if (impact == "donor_gain"){
                        score_list["pos"] <- DP_DG
                        } else {
                        score_list["pos"] <- DP_DL
                        }
                    }
                }
            variants[i, "SpliceAI_impact"] <- paste(score_list["gene"], score_list["impact"], score_list["pos"], sep="|")
            variants[i, "SpliceAI_score"] <- score_list["score"]
            }
        }

    
    # pathogenicity scores
    # Column42 = sift
    # Column43 = polyphen
    # Column44 = cadd
    # Column45 = vest4
    for (i in 1:nrow(variants)){
        v_vest <- strsplit(variants[i,"Vest4_score"], ",", fixed = T)[[1]]
        variants[i, "Vest4_score"] <- max(v_vest)
    }
    
    # Column45 = revel
    
    # Column46 = Gerp
    
    # Column47 = Imprinting_status
    # Column48 = Imprinting_expressed_allele
    imprinting_file_name <- paste0(default_tables_path, "/imprinting.txt")
    imprinting <- read.delim(imprinting_file_name, stringsAsFactors = F)
    variants <- merge(variants, imprinting, all.x = T)
    
    # Column49 - pseudoautosomal
    pseudoautosomal_file_name <- paste0(default_tables_path, "/pseudoautosomal.txt")
    pseudoautosomal <- read.delim(pseudoautosomal_file_name, stringsAsFactors = F)
    variants <- merge(variants, pseudoautosomal, all.x = T)
    
    # Column50 - splicing
    variants <- add_placeholder(variants, "Splicing", "NA")
    if ("spliceregion" %in% colnames(impacts))
    {
        for (i in 1:nrow(variants)){
	          v_id <- variants[i,"Variant_id"]
	          splicing_impacts <- subset(impacts, variant_id == v_id,
	                                    select = c("maxentscan_diff","spliceregion"))
	          splicing_impacts <- subset(splicing_impacts, !is.na(maxentscan_diff))
	          splicing_impacts <- unique(splicing_impacts[order(splicing_impacts$maxentscan_diff),])
	          # capture the absolute difference - very weak site, or very strong site
	          # negative - strong alt, + weak alt.
	
	          s_splicing_field <- 0
	          
	          if (nrow(splicing_impacts) > 0){
	              strongest_alt_site <- head(splicing_impacts, n=1)
	              s_splicing_field <- strongest_alt_site$maxentscan_diff
	          }
	          
	          if (nrow(splicing_impacts) > 1){
	              weakest_alt_site <- tail(splicing_impacts, n=1)
	              s_splicing_field <- paste0(s_splicing_field, ";", weakest_alt_site$maxentscan_diff)
	          }
	          
	          variants[i,"Splicing"] <- s_splicing_field
        }
    }else print("VEP MaxEntScan annotation is missing")

    # Column 51: SpliceAI

    
    # Column 51: number of callers
    variants <- add_placeholder(variants, "Number_of_callers", "Number_of_callers")
    
    # Column 52: Old multiallelic
    variants$Old_multiallelic[variants$Old_multiallelic == "None"] <- "NA"

    # Column 53: UCE_100bp 
    # Column 54: UCE_200bp
    # Column 55: DNaseI_hypersensitive_site  
    # Column 56: CTCF_binding_site
    # Column 57: ENH_cellline_tissue
        
    # replace -1 with 0
    for (field in c("Gnomad_af", "Gnomad_af_popmax")){
        variants[,field] <- with(variants, gsub("-1", "0", get(field), fixed = T))
        variants[,field] <- with(variants, gsub("None", "0", get(field), fixed = T))
    }

    print(sort(colnames(variants)))
    select_and_write2(variants, samples, paste0(family, ".create_report"), type)
    return(variants)
}

# load tables
load_tables <- function(){
  hgmd.csv <- paste0(c4r_database_path,"/hgmd_hg38.csv")
  
  if (file.exists(hgmd.csv)){
    hgmd <- read.csv(hgmd.csv,stringsAsFactors = F,header = F)
    colnames(hgmd) <- c("chrom","pos","HGMD_id","ref","alt","HGMD_gene","HGMD_tag","author",
                        "allname","vol","page","year","pmid")
    hgmd$superindex <- with(hgmd,paste0(chrom,':',pos,'-',ref,'-',alt))
    hgmd$HGMD_ref <- with(hgmd,paste(author,allname,vol,page,year,"PMID:",pmid,sep = ' '))
    hgmd <<- hgmd[,c("superindex","HGMD_id","HGMD_gene","HGMD_tag","HGMD_ref")]
  }else{
    print("No HGMD database")
  }
}

# writes in CSV format
select_and_write2 <- function(variants, samples, prefix, type){
    print(colnames(variants))
    if (type == 'wgs' || type == 'denovo'){
        noncoding_cols <- c("DNaseI_hypersensitive_site", "CTCF_binding_site", "ENH_cellline_tissue", "TF_binding_sites",
                           "GreenDB_variant_type", "GreenDB_closest_gene", "GreenDB_controlled_gene")
        noncoding_scores <- c("ncER_score", "ReMM_score", "LINSIGHT_score")
        }
    else {
        noncoding_cols <- c()
        wgs_counts <- c()
        noncoding_scores <- c()
        }
    variants <- variants[c(c("Position", "UCSC_Link", "GNOMAD_Link", "Ref", "Alt"),
                          c("Gene"),
                          c("gts", "Variation", "Info", "Refseq_change", "Depth", "Quality"),
                          c("Ensembl_gene_id", "Gene_description", "omim_phenotype", "omim_inheritance",
                            "Orphanet", "Clinvar"), c("HGMD_id", "HGMD_gene", "HGMD_tag", "HGMD_ref",
                            "Gnomad_af_popmax", "Gnomad_af", "Gnomad_ac", "Gnomad_hom",
                            "Ensembl_transcript_id", "AA_position", "Exon", "Protein_domains", "rsIDs",
                            "Gnomad_oe_lof_score", "Gnomad_oe_mis_score", "Exac_pli_score", "Exac_prec_score", "Exac_pnull_score",
                            "Conserved_in_30_mammals", "SpliceAI_impact", "SpliceAI_score", "Sift_score", "Polyphen_score", "Cadd_score", "Vest4_score", "Revel_score", "Gerp_score", "AlphaMissense"),
                           noncoding_scores,
                            c("Imprinting_status", "Imprinting_expressed_allele", "Pseudoautosomal", "Gnomad_male_ac",
                            "Old_multiallelic", "UCE_100bp", "UCE_200bp"), noncoding_cols)]
  
    variants <- variants[order(variants$Position),]

    if (type == 'denovo'){
        variants <- variants[variants$C4R_WGS_counts < 10,]
    }
    
    write.csv(variants, paste0(prefix,".csv"), row.names = F)
}

fix_column_name <- function(column_name){
    if(grepl("^[0-9]", column_name)){
        column_name <- paste0("X", column_name)
    }
    return(column_name)
}

# merges ensembl, gatk-haplotype reports
merge_reports <- function(family, samples, type){
    ensemble_file <- paste0(family, ".create_report.csv")
    ensemble <- read.csv(ensemble_file, stringsAsFactors = F)
    ensemble$superindex <- with(ensemble, paste(Position, Ref, Alt, sep = '-'))
    
    for (i in 1:nrow(ensemble)){
        v_impacts <- strsplit(ensemble[i,"Info"], "," , fixed = T)[[1]]
	    for (impact in v_impacts){
            if (grepl(":NM_", impact, fixed = T)){
                v_subimpacts <- strsplit(impact, ":", fixed=T)[[1]]
                ensemble[i,"Refseq_change"] <- paste0(v_subimpacts[3], ":", v_subimpacts[4], ":", v_subimpacts[6])
                break
            }
        }
    }
    
    ensemble_table_file <- paste0(family, ".table")
    if (file.exists(ensemble_table_file)){
        ensemble_table <- read.delim(ensemble_table_file, stringsAsFactors = F)
        ensemble_table$superindex <- with(ensemble_table, paste(paste0(CHROM,":",POS), REF, ALT, sep = '-'))
        ensemble_table[c("CHROM", "POS", "REF", "ALT")] <- NULL
        for (i in 1:nrow(ensemble_table)){
            if(!is.na(ensemble_table[i, "CALLERS"])){
        	      v_callers <- strsplit(ensemble_table[i, "CALLERS"],",")[[1]]
        	      ensemble_table[i, "Number_of_callers"] <- length(v_callers)
            }else ensemble_table[i,"Number_of_callers"] <- NA
        }
        ensemble_table["CALLERS"] <- NULL
        ensemble$Number_of_callers <- NULL
        #two variant callers called one genotype, two another - two genotypes, creates two records at the same site
        ensemble <- merge(ensemble, ensemble_table, by.x = "superindex", 
                          by.y = "superindex",all.x = T, all.y = F)
    }

    gatk_file <- paste0(family,"-gatk-haplotype-annotated-decomposed.table")
    if (file.exists(gatk_file)){
        gatk <- read.delim(gatk_file, stringsAsFactors = F)
        gatk$superindex <- with(gatk, paste(paste0(CHROM, ":", POS), REF, ALT, sep = '-'))
        gatk[c("CHROM","POS","REF","ALT")] <- NULL
    
        ensemble <- merge(ensemble, gatk, by.x = "superindex", by.y = "superindex", all.x = T, all.y = F)
    
        ensemble$Depth <- ensemble$DP
        n_sample <- 1
        prefix <- ""

    
        for (sample in samples){
            ensemble[c("DP", paste0(fix_column_name(sample), ".DP"), 
                       paste0(fix_column_name(sample),".AD"))] <- NULL
        }

      for (sample in samples){
          ensemble[c("DP", paste0(fix_column_name(sample),".DP"), 
                     paste0(fix_column_name(sample),".AO"))] <- NULL
      }
    
      for (sample in samples){
          ensemble[c("TC", paste0(fix_column_name(sample), ".NV"), paste0(fix_column_name(sample),".NR"))] <-  NULL
      }
    }
    # don't apply filter here
    filtered_ensemble = ensemble
    select_and_write2(filtered_ensemble, samples, paste0(family, ".merge_reports"),type)
}


annotate_hgmd <- function(family,samples,type){
  variants <- read.csv(paste0(family, ".merge_reports.csv"), stringsAsFactors = F)
  variants$superindex <- with(variants, paste(Position, Ref, Alt, sep='-'))
  
  variants$HGMD_gene <- NULL
  variants$HGMD_id <- NULL
  variants$HGMD_ref <- NULL
  variants$HGMD_tag <- NULL
  variants <- merge(variants, hgmd, by.x = "superindex", 
                    by.y = "superindex", all.x = T, all.y = F)
  variants$HGMD_gene <- NULL
  
  hgmd.genes <- as.data.frame(unique(sort(hgmd$HGMD_gene)))
  hgmd.genes <- cbind(hgmd.genes, hgmd.genes)
  colnames(hgmd.genes) <- c("index", "HGMD_gene")
  variants <- merge(variants, hgmd.genes, by.x = "Gene", by.y = "index",
                    all.x = T, all.y = F)
  return(variants)
}

# clean report columns for mock vcf
clean_and_output_report <- function(out){
  out$Depth = NA
  out$Quality = NA
  
  column_order = c('Gene','Position','UCSC_Link','GNOMAD_Link','Ref','Alt','gts',
                   'Variation','Info','Refseq_change','Depth','Quality',
                   'Ensembl_gene_id','Gene_description','omim_phenotype','omim_inheritance',
                   'Orphanet','Clinvar','HGMD_id','HGMD_tag','HGMD_ref','HGMD_gene',
                   'Gnomad_af_popmax','Gnomad_af','Gnomad_ac','Gnomad_hom',
                   'Ensembl_transcript_id','AA_position','Exon','Protein_domains','rsIDs',
                   'Gnomad_oe_lof_score','Gnomad_oe_mis_score','Exac_pli_score','Exac_prec_score','Exac_pnull_score',
                   'Conserved_in_30_mammals','SpliceAI_impact','SpliceAI_score','Sift_score','Polyphen_score','Cadd_score',
                   'Vest4_score','Revel_score','Gerp_score', 'AlphaMissense', 'ncER_score', 'ReMM_score',               
                   'LINSIGHT_score', 'Imprinting_status','Imprinting_expressed_allele','Pseudoautosomal',
                   'Gnomad_male_ac','Old_multiallelic','UCE_100bp','UCE_200bp','DNaseI_hypersensitive_site','CTCF_binding_site',
                   'ENH_cellline_tissue','TF_binding_sites')
  out = out[,column_order]
  
  write.csv(out, paste0(family,".", type, ".", datetime,".csv"), row.names = F)
}


library(stringr)
library(data.table)
library(plyr)
library(dplyr)
default_tables_path <- "~/cre/data"

# R substitutes "-" with "." in sample names in columns so fix this in samples.txt
# sample names starting with letters should be prefixed by X in *.table
# for correct processing. most of them start with numbers, and R adds X automatically

args <- commandArgs(trailingOnly = T)
print(args)
family <- args[1]

coding <- if(is.null(args[2])) T else F
coding <- F

type <- if(is.na(args[2])) '' else args[2]

c4r_database_path <- args[3] 

print(paste0("Running cre.vcf2db.R with inputs: ", family, coding, type))
setwd(family)

samples <- unlist(read.table("samples.txt", stringsAsFactors = F))
samples <- gsub("-", ".", samples)
print(samples)

print("Loading tables")
load_tables()
print("Creating report, generates *create_report.csv")
create_report(family,samples,type)
print("Merging reports, merge *create_report.csv with *-gatk-haplotype-annotated-decomposed.table if available. Mock vcf input does not have this table")
merge_reports(family,samples,type)
print("Merging in HGMD data")
out = annotate_hgmd(family,samples,type)
colnames(out)

clean_and_output_report(out)

cat("\n", "Done\n\n")

setwd("..")
