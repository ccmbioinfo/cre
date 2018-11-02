setwd("~/Desktop/work")

# variant_report.csv is the output of cre has Ensembl_gene_id, 
# panel.csv = gene panel with ensembl_gene_id column
filter_variants = function(variant_report.csv,panel.csv,output.csv)
{
    variants = read.csv(variant_report.csv, stringsAsFactors = F)    
    panel = read.csv(panel.csv, stringsAsFactors=F)
    
    variants.panel = variants[variants$Ensembl_gene_id %in% panel$ensembl_gene_id,]
    
    write.csv(variants.panel,output.csv, row.names = F)
}

# when the panel is from genomics England
# https://panelapp.genomicsengland.co.uk/panels/60/
filter_variants_genomics_england_panel = function(variant_report.csv,panel.tsv)
{
    variants = read.csv(variant_report.csv, stringsAsFactors = F)    
    panel = read.delim("Periodic fever syndromes.tsv", stringsAsFactors=F)

    panel = panel[,c("Gene.Symbol","Model_Of_Inheritance","Phenotypes","UserRatings_Green_amber_red","EnsemblId.GRch37.")]

    colnames(panel)=c("PanelAPP.Gene.Symbol","PanelAPP.Model_Of_Inheritance","PanelAPP.Phenotypes",
                  "PanelAPP.UserRatings_Green_amber_red", "PanelAPP.EnsemblId.GRch37")

    variants.ens = variants[variants$Ensembl_gene_id %in% panel$PanelAPP.EnsemblId.GRch37,]
    variants.gene = variants[variants$Gene %in% panel$PanelAPP.Gene.Symbol,]
}