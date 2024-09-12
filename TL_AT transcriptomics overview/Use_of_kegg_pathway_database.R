library(limma)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(KEGGREST)

kegg_pathway <- getKEGGPathwayNames("hsa")

gene_list <- getGeneKEGGLinks("hsa")
gene_list$Symbol <- mapIds(
  org.Hs.eg.db,
  gene_list$GeneID,
  column = "SYMBOL",
  keytype = "ENTREZID"
)

sum(duplicated(gene_list$GeneID)) #28,215

length(unique(gene_list$PathwayID)) #359