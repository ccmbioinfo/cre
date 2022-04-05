# variant report generator

# store date to be used when writing files
#datetime <- format(Sys.time(),"%Y-%m-%d_%H-%M") # can use to get timestamp for testing
datetime <- format(Sys.time(),"%Y-%m-%d")

# Rscript ~/cre/cre.vcf2.db.R <family> noncoding|default=NULL,coding <database path>
add_placeholder <- function(variants, column_name, placeholder){
    variants[,column_name] <- with(variants, placeholder)
    return(variants)
}

get_variants_from_file <- function (filename){
    variants <- read.delim(filename, stringsAsFactors = F)
    return(variants)
}

# returns Hom / Het / - (for HOM reference)
genotype2zygocity <- function (genotype_str, ref, alt_depth){
    # test
    # genotype_str = "A|A|B"
    # genotype_str = "./." - call not possible
    # genotype_str = "TCA/."
    # genotype_str = "G"
    # genotype_str = "A/A"
    # greedy
    genotype_str <- gsub("|", "/", genotype_str, fixed = T)
    genotype_str <- gsub("./.", "Insufficient_coverage", genotype_str, fixed = T)
    #genotype_str <- gsub("/.","NO_CALL",genotype_str,fixed=T)
      
    if(grepl("Insufficient_coverage", genotype_str)){
      result <- genotype_str
    }else if(alt_depth == 0){
      result <- '-'
    }else{
        ar <- strsplit(genotype_str, "/", fixed = T)
        len <- length(ar[[1]])
        if (len == 2){
            if (ar[[1]][1] == ar[[1]][2]){
                if (ar[[1]][1] == ref)
                    result <- "-"
                else
                    result <- "Hom"
            }else result <- "Het"
        }else result <- genotype_str
    }
    return(result)
}

# output : family.ensemble.txt
create_report <- function(family, samples, type){
    file <- paste0(family, ".variants.txt")
    variants <- get_variants_from_file(file)
    
    impact_file <- paste0(family, ".variant_impacts.txt")
    impacts <- get_variants_from_file(impact_file)
    
    # temporarily due to https://github.com/quinlan-lab/vcf2db/issues/48
    # fixed in vcf2db
    #transcripts_genes = read.csv("~/cre/data/genes.transcripts.csv")
    #variants$Ensembl_gene_id=NULL
    #variants$Ensembl_transcript_id1=variants$Ensembl_transcript_id
    
    #for (i in 1:nrow(variants)){
    #    variant_id = variants[i,"Variant_id"]
    #    #variant_id="4489"
    #    variant_impacts = subset(impacts, variant_id == variant_id & gene == variants[i,"Gene"])
    #    variant_impacts = variant_impacts[order(variant_impacts$transcript),]
        
    #    if (nrow(variant_impacts)>0){
    #        variants[i,"Ensembl_transcript_id1"] = variant_impacts[1,"transcript"]
    #    }
    #    ar = strsplit(variants[i,"Ensembl_transcript_id1"],".",fixed=T)
    #    variants[i,"Ensembl_transcript_id1"] = ar[[1]][1]
    #}
    
    #variants = merge(variants,transcripts_genes,by.x='Ensembl_transcript_id1',by.y="Ensembl_transcript_id",all.x=T)
    
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
    
    # Column 6 - Zygosity, column 8 - Burden
    # use new loader vcf2db.py - with flag  to load plain text
    # for genotype and depth - Noah
    # otherwise have to decode BLOB 
    # snappy decompression
    # https://github.com/arq5x/gemini/issues/700
    # https://github.com/lulyon/R-snappy
    for(sample in samples){
        #DEBUG: gene = IL20RA
        #sample=samples[1]
        
        zygocity_column_name <- paste0("Zygosity.", sample)
        #t = lapply(variants[,paste0("gts.",sample),"Ref"],genotype2zygocity)
        #t = lapply(variants[,paste0("gts.",sample),"Ref"],genotype2zygocity)
        t <- unlist(mapply(genotype2zygocity, variants[,paste0("gts.",sample)], 
                           variants[,"Ref"], variants[,paste0("gt_alt_depths.",sample)]))
        variants[,zygocity_column_name] <- unlist(t)
    
        burden_column_name <- paste0("Burden.", sample)
        # calculating Burden using gene rather then Ensembl_gene_id - request from Matt
        t <- subset(variants, 
                    get(zygocity_column_name) == 'Hom' | get(zygocity_column_name) == 'Het',
                    select = c("Gene", zygocity_column_name))
        # counts from plyr
        df_burden <- plyr::count(t, "Gene")
	df_burden <- subset(df_burden, Gene!='None')    
        colnames(df_burden)[2] <- burden_column_name
        variants <- merge(variants, df_burden, all.x = T)
        variants[,burden_column_name][is.na(variants[,burden_column_name])] <- 0
        variants[,burden_column_name][is.na(variants$Gene)] <-0
    }
    
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

    # Column 15 - Alt_depth - from v.gt_alt_depths
    # when multiple callers used, AD is not set and fixed in merge_reports function
    for(sample in samples){
        new_name <- paste0("Alt_depths.", sample)
        setnames(variants, paste0("gt_alt_depths.", sample), new_name)
    }

    # Column 16 - Trio_coverage - fixed in merge_reports function
    variants <- add_placeholder(variants, "Trio_coverage", "")
    n_sample <- 1
    prefix <- ""
    
    #order gts column in the same way as in samples
    variants$gts <- ""
    for(sample in samples){
        column <- paste0("gt_depths.", sample)
        
        if (n_sample>1) prefix <- "/"
        
        variants$Trio_coverage <- with(variants, paste0(Trio_coverage, prefix, get(column)))
      
        column <- paste0("gts.", sample)
        
        if (n_sample>1) prefix <- ","
        
        variants$gts <- with(variants,paste0(gts, prefix, get(column)))
      
        n_sample <- n_sample+1
    }
    
    # Column17 = Ensembl_gene_id

    # Column18 = Gene_description
    gene_descriptions <- read.delim2(paste0(default_tables_path,"/ensembl_w_description.txt"), 
                                     stringsAsFactors = F)
    variants <- merge(variants, gene_descriptions, by.x = "Ensembl_gene_id",
                      by.y = "ensembl_gene_id", all.x = T)
    
    # Column19 - Omim_phenotype
    # Column20 - Omim_inheritance 
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
    
    # Column 27, 28 = C4R_WES_counts, C4R_WES_samples
    variants <- add_placeholder(variants, "C4R_WES_counts", "C4R_WES_counts")
    variants <- add_placeholder(variants, "C4R_WES_samples", "C4R_WES_samples")
    
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
    for (field in c("Trio_coverage", "Gnomad_af", "Gnomad_af_popmax")){
        variants[,field] <- with(variants, gsub("-1", "0", get(field), fixed = T))
        variants[,field] <- with(variants, gsub("None", "0", get(field), fixed = T))
    }

    for (field in c(paste0("Alt_depths.",samples))){
        variants[,field] <- with(variants, gsub("-1", NA, get(field), fixed = T))  
    }

    print(sort(colnames(variants)))
    select_and_write2(variants, samples, paste0(family, ".create_report"), type)
}

# writes in CSV format
select_and_write2 <- function(variants, samples, prefix, type)
{
    print(colnames(variants))
    if (type == 'wgs' || type == 'denovo'){
        noncoding_cols <- c("DNaseI_hypersensitive_site", "CTCF_binding_site", "ENH_cellline_tissue", "TF_binding_sites")
        wgs_counts <- c("C4R_WGS_counts", "C4R_WGS_samples")
        variants$C4R_WGS_counts[variants$C4R_WGS_counts == "None"] <- 0 
        variants$C4R_WGS_counts <- as.integer(variants$C4R_WGS_counts)
        variants$C4R_WGS_samples[variants$C4R_WGS_samples == "None"] <- 0
        }
    else {
        noncoding_cols <- c()
        wgs_counts <- c()
        }
    variants <- variants[c(c("Position", "UCSC_Link", "GNOMAD_Link", "Ref", "Alt"),
                          paste0("Zygosity.", samples),
                          c("Gene"),
                          paste0("Burden.", samples),
                          c("gts", "Variation", "Info", "Refseq_change", "Depth", "Quality"),
                          paste0("Alt_depths.", samples),
                          c("Trio_coverage", "Ensembl_gene_id", "Gene_description", "omim_phenotype", "omim_inheritance",
                            "Orphanet", "Clinvar",
                            "C4R_WES_counts", "C4R_WES_samples"), wgs_counts, c("HGMD_id", "HGMD_gene", "HGMD_tag", "HGMD_ref",
                            "Gnomad_af_popmax", "Gnomad_af", "Gnomad_ac", "Gnomad_hom",
                            "Ensembl_transcript_id", "AA_position", "Exon", "Protein_domains", "rsIDs",
                            "Gnomad_oe_lof_score", "Gnomad_oe_mis_score", "Exac_pli_score", "Exac_prec_score", "Exac_pnull_score",
                            "Conserved_in_30_mammals", "SpliceAI_impact", "SpliceAI_score", "Sift_score", "Polyphen_score", "Cadd_score", "Vest4_score", "Revel_score", "Gerp_score",
                            "Imprinting_status", "Imprinting_expressed_allele", "Pseudoautosomal", "Gnomad_male_ac",
                            "Number_of_callers", "Old_multiallelic", "UCE_100bp", "UCE_200bp"), noncoding_cols)]
  
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
        ensemble$Trio_coverage <- ""
    
        for(sample in samples){
            #R fixes numerical column names with X?
            #what if sample is not numerical
            column <- fix_column_name(sample)  
            column <- paste0(column,".DP")
        
	        #prefix changed to _ from / because otherwise excel converts the field into date
            if (n_sample > 1) prefix <- "_"
        
            ensemble$Trio_coverage <- with(ensemble, paste0(Trio_coverage, prefix, get(column)))
      
            column <- paste0("Alt_depths.", sample)
            column_gatk <- fix_column_name(sample)
            column_gatk <- paste0(column_gatk, ".AD")
        
            ensemble[,column] <- ensemble[,column_gatk]
      
            n_sample <- n_sample + 1
        }
    
        for (i in 1:nrow(ensemble)){
            for (sample in samples){
                field <- paste0("Alt_depths.", sample)
                #when combining reports from vcfs called elsewere there may be no AD field, just -1
                if (grepl(",", ensemble[i,field])){
            	    ensemble[i, field] <- strsplit(ensemble[i,field], ",", fixed = T)[[1]][2]
            	}
            }
        }
    
        for (sample in samples){
            ensemble[c("DP", paste0(fix_column_name(sample), ".DP"), 
                       paste0(fix_column_name(sample),".AD"))] <- NULL
        }
    }

    freebayes_file <- paste0(family,"-freebayes-annotated-decomposed.table")
    if(file.exists(freebayes_file)){
        freebayes <- read.delim(freebayes_file, stringsAsFactors = F)
        freebayes$superindex <- with(freebayes, paste(paste0(CHROM,":",POS), REF, ALT, sep = '-'))
        freebayes[c("CHROM","POS","REF","ALT")] <- NULL

        ensemble <- merge(ensemble, freebayes, by.x = "superindex", 
                          by.y = "superindex", all.x = T, all.y = F)
        for (i in 1:nrow(ensemble)){
            #if(grepl("NA",ensemble[i,"Trio_coverage"]))
            #wrong: a variant may be called by gatk with 10/10/NA,
            #and freebayes will destroy coverage info
            if (str_count(ensemble[i,"Trio_coverage"], "NA") == length(samples)){
                ensemble[i, "Depth"] <- ensemble[i, "DP"]
                for (sample in samples){
                    field_depth <- paste0("Alt_depths.", sample)
                    field_bayes <- paste0(fix_column_name(sample), ".AO")
                    #field_bayes = paste0(sample,".AO")
                    ensemble[i, field_depth] <- ensemble[i, field_bayes]
                }
                n_sample <- 1
                prefix <- ""
                ensemble[i, "Trio_coverage"] <- ""
            
                for(sample in samples){
                    column <- paste0(fix_column_name(sample),".DP")
                    if (n_sample > 1) prefix <- "_"
                        ensemble[i, "Trio_coverage"] <- paste(ensemble[i,"Trio_coverage"],
                                                              ensemble[i,column], sep = prefix)
              
                    n_sample <- n_sample+1
                }
          }
      }
      for (sample in samples){
          ensemble[c("DP", paste0(fix_column_name(sample),".DP"), 
                     paste0(fix_column_name(sample),".AO"))] <- NULL
      }
    }

    platypus_file <- paste0(family, "-platypus-annotated-decomposed.table")
    if(file.exists(platypus_file)){
        platypus <- read.delim(platypus_file, stringsAsFactors = F)
        platypus$superindex <- with(platypus, paste(paste0(CHROM,":",POS), REF, ALT, sep = '-'))
        platypus[c("CHROM", "POS", "REF", "ALT")] <- NULL
        ensemble <- merge(ensemble, platypus, by.x = "superindex", by.y = "superindex", 
                          all.x = T, all.y = F)
    
        for (i in 1:nrow(ensemble)){
          #if(grepl("NA",ensemble[i,"Trio_coverage"])) - wrong, may be 10/10/NA in gatk
          #if (ensemble[i,"Trio_coverage"]=="NA/NA/NA")
            if (str_count(ensemble[i,"Trio_coverage"],"NA") == length(samples)){
                ensemble[i,"Depth"] <- ensemble[i,"TC"]
                for (sample in samples){
                    field_depth <- paste0("Alt_depths.", sample)
                    field_bayes <- paste0(fix_column_name(sample), ".NV")
          
                    #sometimes freebayes has 10,10,10 for decomposed alleles
                    if (grepl(",", ensemble[i,field_bayes])){
                        ensemble[i,field_depth] <- strsplit(ensemble[i,field_bayes], ",", fixed = T)[[1]][1]
                    }else{ #teja
                    ensemble[i, field_depth] <- ensemble[i, field_bayes]
                    } #teja
                }
                n_sample <- 1
                prefix <- ""
                ensemble[i, "Trio_coverage"] <- ""
        
                for(sample in samples){
                    column <- paste0(fix_column_name(sample), ".NR")
                    if (n_sample > 1) prefix <- "_"
                    #sometimes freebayes has 10,10,10 for decomposed alleles
                    if (grepl(",",ensemble[i,column])){
                        cov_value <- strsplit(ensemble[i,column], ",", fixed = T)[[1]][1]
                    }else cov_value <- ensemble[i,column]
                    
                    ensemble[i, "Trio_coverage"] <- paste(ensemble[i, "Trio_coverage"], cov_value, sep = prefix)
                    n_sample <- n_sample + 1
                }
          }
      }
    
      for (sample in samples){
          ensemble[c("TC", paste0(fix_column_name(sample), ".NV"), paste0(fix_column_name(sample),".NR"))] <-  NULL
      }
    }
    
    #don't use samtools file by default!
    samtools_file <- paste0(family,"-samtools-annotated-decomposed.table")
    if(file.exists(samtools_file)){
        samtools <- read.delim(samtools_file, stringsAsFactors = F)
        samtools$superindex <- with(samtools, paste(paste0(CHROM, ":", POS), REF, ALT, sep = '-'))
        samtools[c("CHROM", "POS", "REF", "ALT")] = NULL
        ensemble <- merge(ensemble, samtools, by.x = "superindex", 
                         by.y="superindex", all.x = T, all.y = F)
      
        for (i in 1:nrow(ensemble)){
            ensemble[i, "Depth"] = ensemble[i,"DP"]
            for (sample in samples){
                field_depth <- paste0("Alt_depths.", sample)
                field_samtools <- paste0(fix_column_name(sample), ".DP")
                ensemble[i, field_depth] <- ensemble[i, field_samtools]
            }
            ensemble[i, "Trio_coverage"] <- ""
        }
        for (sample in samples){
          ensemble[c("DP", paste0(fix_column_name(sample),".DP"))] <- NULL
          #samtools does not discriminate between insufficient coverage (cannot call) and no_call =reference
          field <- paste0("Zygosity.", sample)
          ensemble[,field] <- with(ensemble, gsub("Insufficient_coverage", 
                                   "-", get(field), fixed=T))
        }
    }
    
    ensemble[,"Trio_coverage"] <- with(ensemble,gsub("NA", "0", get("Trio_coverage"), fixed = T))  
   
    for (i in 1:nrow(ensemble)){
        if (is.na(ensemble[i, "Depth"])){
            l <- strsplit(ensemble[i, "Trio_coverage"],"_")[[1]]
            ensemble[i, "Depth"] <- sum(as.integer(l))
        }
        sample_index <- 1
        for (sample in samples){
            field_depth <- paste0("Alt_depths.", sample)
            parsed_alt_depth <- parse_ad(ensemble[i,field_depth])
            ensemble[i,field_depth] <- parsed_alt_depth
            # fix the zygosity after the alternate depths are set
            # if ad is 0 make zygosity -
            zygocity_column_name <- paste0("Zygosity.", sample)
            # split by comma to grab the sample's gt
            #gts <- data.frame(do.call('rbind', strsplit(as.character(ensemble$gts), ",",fixed=TRUE)))
            #sample_gt <- gts[]
            #print(i)
            #print(ensemble[i,"Position"])
            #print(sample_index)
            #print("before")
            #print(ensemble[i,zygocity_column_name])
            #print("gt")
            gts <- unlist(strsplit(ensemble[i,"gts"],","))
            #print(gts[sample_index])
            fixed_zygosity <- genotype2zygocity(gts[sample_index],ensemble[i,"Ref"],ensemble[i,field_depth])
            #print("after")
            #print(fixed_zygosity)
            ensemble[i,zygocity_column_name] <- fixed_zygosity
            sample_index <- sample_index + 1
        }
    }

    # after the alt depths columns are fixed, remove all variants that don't pass the alt depth >= 3 filter
    filtered_ensemble <- dplyr::filter_at(ensemble, paste0("Alt_depths.",samples), any_vars(as.integer(.) >= 3))
    select_and_write2(filtered_ensemble, samples, paste0(family, ".merge_reports"),type)	
}

parse_ad <- function(ad_cell) {
  if (is.na(ad_cell)){
    alt_depth <- as.integer(0)
  }
  else if (grepl(",",ad_cell)){
    # there can be multiple ad values reported here. use the largest
    alt_depths <- unlist(strsplit(ad_cell, ","))
    alt_depth <- 0
    for (a in alt_depths) {
      if (!(is.na(a)) && (as.integer(a) > alt_depth)){alt_depth <- as.integer(a)}
    }
  }
  else{alt_depth <- as.integer(ad_cell)}
  return(alt_depth)
}

annotate_w_care4rare <- function(family,samples,type){
    variants <- read.csv(paste0(family, ".merge_reports.csv"), stringsAsFactors = F)
  
    variants$superindex <- with(variants, paste(Position, Ref, Alt, sep='-'))
    
    if(exists("seen_in_c4r_counts")){
        variants <- merge(variants, seen_in_c4r_counts, by.x = "superindex", 
                          by.y = "Position.Ref.Alt", all.x = T)
        variants$C4R_WES_counts <- variants$Frequency
        variants$Frequency <- NULL
    }
    
    variants$C4R_WES_counts[is.na(variants$C4R_WES_counts)] <- 0
    
    if(exists("seen_in_c4r_samples")){
        variants <- merge(variants,seen_in_c4r_samples,by.x = "superindex", 
                          by.y = "Position.Ref.Alt", all.x = T)
        variants$C4R_WES_samples <- variants$Samples
    }
    
    variants$C4R_WES_samples[is.na(variants$C4R_WES_samples)] <- 0        
		# truncate column if it has more than 30000 variants
		variants$C4R_WES_samples <- strtrim(variants$C4R_WES_samples, 30000)
		    
    if (exists("hgmd")){
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
    }
    
    select_and_write2(variants, samples, paste0(family, ".", type, ".", datetime), type)
}

load_tables <- function(debug = F){
    print(paste0("Debug:", debug))
    #debug
    if (debug == T){
        seen_in_c4r_counts.txt <- "seen_in_c4r_counts.txt"
        seen_in_c4r_samples.txt <- "seen_in_c4r_samples.txt"
        hgmd.csv <- "hgmd_hg38.csv"
    }else{
        seen_in_c4r_counts.txt <- paste0(c4r_database_path,"/seen_in_c4r_counts.txt")    
        seen_in_c4r_samples.txt <- paste0(c4r_database_path,"/seen_in_c4r_samples.txt")
        hgmd.csv <- paste0(c4r_database_path,"/hgmd_hg38.csv")
    }
    
    if (file.exists(seen_in_c4r_counts.txt)){
        seen_in_c4r_counts <<- read.delim(seen_in_c4r_counts.txt, stringsAsFactors=F)
    }else{
        print("No C4R counts found")
    }
    
    if (file.exists(seen_in_c4r_samples.txt)){
        seen_in_c4r_samples <<- read.delim(seen_in_c4r_samples.txt, stringsAsFactors=F)
    }else{
        print("No C4R samples found")
    }
    
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

# creates clinical report - more conservative filtering and less columns
clinical_report <- function(project,samples,type){
    report_file_name <- paste0(project, ".", type, ".", datetime,".csv")
    full_report <- read.csv(report_file_name, header = T, stringsAsFactors = F)
    
    full_report$max_alt <- with(full_report, pmax(get(paste0("Alt_depths.", samples))))

    for (i in 1:nrow(full_report)){
        for (sample in samples){
            field_depth <- paste0("Alt_depths.", sample)
            parsed_alt_depth <- parse_ad(full_report[i,field_depth])
            full_report[i,field_depth] <- parsed_alt_depth
        }
    }

    # for clinical, only keep variants where one of the alt depths was >= 20
    full_report <- dplyr::filter_at(full_report, paste0("Alt_depths.",samples), any_vars(as.integer(.) >= 20))    
    filtered_report <- subset(full_report, 
               Quality > 1000 & Gnomad_af_popmax < 0.005 & C4R_WES_counts < 6,
               select = c("Position", "GNOMAD_Link", "Ref", "Alt", "Gene", paste0("Zygosity.", samples), 
                          paste0("Burden.",samples),
                          paste0("Alt_depths.",samples),
                        "Variation", "Info", "Refseq_change", "omim_phenotype", "omim_inheritance",
                        "Orphanet", "Clinvar", "C4R_WES_counts",
                        "Gnomad_af_popmax", "Gnomad_af", "Gnomad_ac", "Gnomad_hom",
                        "Sift_score", "Polyphen_score", "Cadd_score", "Vest4_score", "Revel_score",
                        "Imprinting_status", "Pseudoautosomal", "Gnomad_male_ac", "UCE_100bp","UCE_200bp")
               )
    
    # recalculate burden using the filtered report
    for(sample in samples){
        zygosity_column_name <- paste0("Zygosity.", sample)
        burden_column_name <- paste0("Burden.", sample)
        t <- subset(filtered_report, 
                    get(zygosity_column_name) == 'Hom' | get(zygosity_column_name) == 'Het',
                    select = c("Gene", zygosity_column_name))
        # count is from plyr
        df_burden <- plyr::count(t, "Gene")    
        colnames(df_burden)[2] <- burden_column_name
        filtered_report[,burden_column_name] <- NULL
        filtered_report <- merge(filtered_report, df_burden, all.x = T)
        filtered_report[,burden_column_name][is.na(filtered_report[, burden_column_name])] <- 0
        filtered_report[,burden_column_name][is.na(filtered_report$Gene)] <-0
    }
    
    #order columns
    filtered_report <- filtered_report[c("Position", "GNOMAD_Link", "Ref", "Alt", "Gene", paste0("Zygosity.", samples), 
      paste0("Burden.", samples),
      "Variation", "Info", "Refseq_change", "omim_phenotype", "omim_inheritance",
      "Orphanet", "Clinvar", "C4R_WES_counts",
      "Gnomad_af_popmax", "Gnomad_af", "Gnomad_ac", "Gnomad_hom",
      "Sift_score", "Polyphen_score", "Cadd_score", "Vest4_score", "Revel_score",
      "Imprinting_status", "Pseudoautosomal", "Gnomad_male_ac", "UCE_100bp", "UCE_200bp")]

    write.csv(filtered_report, paste0(project, ".clinical.", type, ".", datetime, ".csv"), row.names = F)
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

debug <- F

print(paste0("Running cre.vcf2db.R with inputs: ", family, coding, type))
setwd(family)

samples <- unlist(read.table("samples.txt", stringsAsFactors = F))
samples <- gsub("-", ".", samples)

print("Loading tables")
load_tables(debug)
print("Creating report")
create_report(family,samples,type)
print("Merging reports")
merge_reports(family,samples,type)
print("Annotating Reports")
annotate_w_care4rare(family,samples,type)
print("Writing Clinical Report")
clinical_report(family,samples,type)

setwd("..")
