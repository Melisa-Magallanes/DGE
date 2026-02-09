library(readr)
library(DGEobj.utils)
library(dplyr)
library(stringr)
library(textshape)
library(edgeR)

## ====== Paths ======
setwd("")

## ====== Leer conteos ======
mercenariafc <- read.delim("counts_merc_30vs10.txt", comment.char = "#", check.names = FALSE)
studydesign <- read.delim("studydesign.txt")

countdata <- mercenariafc %>%
  column_to_rownames("Geneid") %>% # turn the geneid column into rownames
  rename_all(str_remove, ".bam") %>% # remove the ".bam" from the column names
  rename_all(str_remove, "_sorted") %>%
  select(-Length)%>%
  select(-Chr)%>%
  select(-Start)%>%
  select(-End)%>%
  select(-Strand)%>%
  as.matrix()

librarySizes <- colSums(countdata)
barplot(librarySizes, 
        names=names(librarySizes), 
        las=2, 
        main="Barplot of library sizes")

sampleLabels <- studydesign$Sample
myDGEList <- DGEList(countdata)
cpm <- cpm(myDGEList)

colSums(cpm)
library(tidyr)
log2.cpm <- cpm(myDGEList, log=TRUE)

log2.cpm.df <- as_tibble(log2.cpm, rownames = "geneID")
log2.cpm.df

table(rowSums(myDGEList$counts==0)==6)
keepers <- rowSums(cpm>1)>=3
myDGEList.filtered <- myDGEList[keepers,]
dim(myDGEList.filtered)
log2.cpm.filtered <- cpm(myDGEList.filtered, log=TRUE)
log2.cpm.filtered.df <- as_tibble(log2.cpm.filtered, rownames = "geneID")
myDGEList.filtered.norm <- calcNormFactors(myDGEList.filtered, method = "TMM")
log2.cpm.filtered.norm <- cpm(myDGEList.filtered.norm, log=TRUE)
cpm.filtered.norm <- cpm(myDGEList.filtered.norm)
log2.cpm.filtered.norm.df <- as_tibble(log2.cpm.filtered.norm, rownames = "geneID")

group <- factor(studydesign$Condition)
pca.res <- prcomp(t(log2.cpm.filtered.norm), scale.=F, retx=T)
#look at the PCA result (pca.res) that you just created
ls(pca.res)
summary(pca.res) # Prints variance summary for all principal components.
pca.res$rotation #$rotation shows you how much each gene influenced each PC (called 'scores')
pca.res$x # 'x' shows you how much each sample influenced each PC (called 'loadings')
#note that these have a magnitude and a direction (this is the basis for making a PCA plot)
screeplot(pca.res) # A screeplot is a standard way to view eigenvalues for each PCA
pc.var<-pca.res$sdev^2 # sdev^2 captures these eigenvalues from the PCA result
pc.per<-round(pc.var/sum(pc.var)*100, 1) # we can then use these eigenvalues to calculate the percentage variance explained by each PC
pc.per

pca.res.df <- as_tibble(pca.res$x)

library(ggplot2)
ggplot(pca.res.df) +
  aes(x=PC1, y=PC2) +
  aes(x=PC1, y=PC2, label=sampleLabels, fill = group) +
  #geom_point(aes(shape = specie, color = group), size = 8)+ 
  geom_label() +
  geom_text()+
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot") +
  theme_bw()

mydata.df <- log2.cpm.filtered.norm.df %>% 
  mutate(Control.AVG = (SRR24733783	+SRR24733767	+SRR24733782)/3,
         Low.AVG = (SRR24733771+ SRR24733770+ SRR24733764)/3,
         #now make columns comparing each of the averages above that you're interested in
         LogFC = (Control.AVG - Low.AVG)) %>% 
  mutate_if(is.numeric, round, 2)

design <- model.matrix(~0+group)

v.DEGList.filtered.norm <- voom(myDGEList.filtered.norm, design, plot = TRUE)
# fit a linear model to your data
fit <- lmFit(v.DEGList.filtered.norm, design)

# Contrast matrix ----
contrast.matrix <- makeContrasts(adaptation1 = groupCONTROL - groupLOW,
                                 levels=design)

fits <- contrasts.fit(fit, contrast.matrix)
#get bayesian stats for your linear model fit
ebFit <- eBayes(fits)

myTopHits <- topTable(ebFit, adjust ="BH", coef=1, number=40000, lfc= 1,sort.by="logFC", p.value= 0.05)

myTopHits.df <- myTopHits %>%
  as_tibble(rownames = "geneID")


write.table(myTopHits.df$geneID, "DEGs_mercenaria_geneIDs.txt",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

diamond <- read.delim("diamond.tsv", header = FALSE)
colnames(diamond) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                       "qstart", "qend", "sstart", "send", "evalue", "bitscore")

# Quedarte con el mejor hit por bitscore
diamond_best <- diamond %>%
  group_by(qseqid) %>%
  arrange(desc(bitscore), evalue) %>%
  slice(1) %>%
  ungroup()

##Extraer nombre de la proteína de UniProt
diamond_best <- diamond_best %>%
  mutate(
    # Extrae el ID UniProt (entre el primer y segundo "|")
    uniprot_id = sub("^[a-z]{2}\\|([^|]+)\\|.*", "\\1", sseqid),
    
    # Extrae el nombre de la proteína (parte entre el segundo "|" y el "_")
    protein_name = sub("^[a-z]{2}\\|[^|]+\\|([^_]+)_.*", "\\1", sseqid),
    
    # Extrae el código de especie (lo que está después del "_")
    species = sub(".*_([A-Z0-9]+)$", "\\1", sseqid)
  )

loc2best <- read.delim("loc2bestprot.tsv", header = FALSE,
                       col.names = c("geneID", "protein_id"))

deg_annot <- myTopHits.df %>%
  left_join(loc2best, by = "geneID") %>%
  left_join(diamond_best, by = c("protein_id" = "qseqid")) %>%
  arrange(adj.P.Val, desc(abs(logFC)))

uniprot_info <- read_tsv("D:/PASSER/BIVALVES/TRANSCRIPTOMAS/CGIGAS/uniprotkb_bivalvia.tsv")
uniprot_info <- uniprot_info %>%
  rename(
    uniprot_id = Entry,                      # o ntry según tu archivo
    entry_name = `Entry Name`,
    protein_names = `Protein names`,
    gene_names = `Gene Names`,
    organism = Organism,
    reviewed = Reviewed
  )

annotated_DEGs <- deg_annot %>%
  left_join(uniprot_info, by = "uniprot_id")

# Guardar resultado final
write.table(annotated_DEGs, "mercenaria_DEGs_with_Uniprot.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
# --- Paquetes ---
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr); library(tidyr)
  library(clusterProfiler)  # para GSEA
  library(fgsea)         # opcional
})


# Helpers para detectar columnas con nombres distintos
pick_col <- function(candidates, df=annotated_DEGs, required=TRUE){
  found <- candidates[candidates %in% names(df)]
  if (length(found)==0) {
    if (required) stop("No encuentro columnas: ", paste(candidates, collapse=", "))
    return(NA_character_)
  }
  found[1]
}

col_uniprot <- pick_col("uniprot_id")
col_logfc   <- pick_col("logFC")
col_q       <- pick_col("adj.P.Val")
col_go_ids  <- pick_col(c("Gene Ontology IDs"), required=FALSE)


# --- 2) Limpiar/estandarizar GO (IDs y, si hay, nombres) ---
if (is.na(col_go_ids) && is.na(col_go_nm)) {
  stop("No encuentro columnas de GO en tu archivo. Exportá desde UniProt con 'Gene Ontology (GO) IDs' y/o 'Gene Ontology (GO)'.")
}

df_go <- annotated_DEGs %>%
  transmute(
    uniprot_id = .data[[col_uniprot]],
    go_ids_raw = if (!is.na(col_go_ids)) .data[[col_go_ids]] else NA_character_)

# separadores posibles ; , |
split_any <- function(x) str_split(x, pattern = "\\s*[;,|]\\s*")

df_go_long <- annotated_DEGs %>%
  select(uniprot_id = all_of(col_uniprot), go_ids_raw = all_of(col_go_ids)) %>%
  filter(!is.na(go_ids_raw), go_ids_raw != "") %>%
  mutate(go_ids_list = str_split(go_ids_raw, "\\s*[;,|]\\s*")) %>%
  unnest_longer(go_ids_list, values_to = "go_id") %>%
  mutate(go_id = trimws(go_id)) %>%
  filter(str_detect(go_id, "^GO:")) %>%
  distinct(uniprot_id, go_id)


# --- 3) Construir TERM2GENE y TERM2NAME (para clusterProfiler::GSEA) ---
TERM2GENE <- df_go_long %>% select(go_id, uniprot_id) %>% distinct()
TERM2NAME <- TERM2GENE %>% select(go_id) %>% distinct() %>%
  mutate(go_name = go_id)


library(GO.db); library(AnnotationDbi)
go_names2 <- AnnotationDbi::select(GO.db,
                                   keys = unique(TERM2NAME$go_id),
                                   keytype = "GOID", columns = "TERM") %>%
  transmute(go_id = GOID, go_term = TERM)

TERM2NAME <- TERM2NAME %>%
  left_join(go_names2, by = "go_id") %>%
  mutate(go_name = dplyr::coalesce(go_term, go_name, go_id)) %>%
  select(go_id, go_name) %>%
  distinct()


# --- 4) Vector rankeado para GSEA (signed score) ---
rank_tbl <- annotated_DEGs %>%
  filter(!is.na(.data[[col_uniprot]]),
         is.finite(.data[[col_logfc]]),
         is.finite(.data[[col_q]])) %>%
  transmute(
    gene  = .data[[col_uniprot]],
    score = sign(.data[[col_logfc]]) * -log10(pmax(.data[[col_q]], 1e-300))
  ) %>%
  group_by(gene) %>% summarise(score = mean(score), .groups = "drop")

geneVec <- rank_tbl$score
names(geneVec) <- rank_tbl$gene
geneVec <- sort(geneVec, decreasing = TRUE)

# --- 5) Correr GSEA (clusterProfiler) ---
gsea <- GSEA(geneList   = geneVec,
             TERM2GENE  = TERM2GENE,
             TERM2NAME  = TERM2NAME,
             pvalueCutoff = 0.1,
             minGSSize  = 5,
             maxGSSize  = 500,
             verbose    = FALSE)

gsea_df <- as.data.frame(gsea)
write.table(gsea_df, "GSEA_GO_results.tsv", sep="\t", quote=FALSE, row.names=FALSE)

