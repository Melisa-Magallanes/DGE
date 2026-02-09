## ====== Paquetes ======
library(dplyr)
library(readr)
library(stringr)
library(DEGseq)
library(ggplot2)
library(tidyr)

## ====== Paths ======
setwd("")

## ====== Leer conteos ======
cgigasfc <- read.delim("counts_cgigas2.txt", comment.char = "#", check.names = FALSE)

# Dejar sólo columnas necesarias: Geneid + dos muestras
counts <- cgigasfc %>%
  select(Geneid, `SRR334273_sorted.bam`, `SRR334269_sorted.bam`) %>%
  rename(LOW = `SRR334269_sorted.bam`,
         HIGH = `SRR334273_sorted.bam`)

# Asegurar enteros no negativos
counts$LOW  <- pmax(0L, as.integer(round(counts$LOW)))
counts$HIGH <- pmax(0L, as.integer(round(counts$HIGH)))

# QC rápido de tamaños de librería
libSizes <- colSums(counts[, c("LOW","HIGH")])
barplot(libSizes, las=2, main="Library sizes (raw counts)")

## ====== (Opcional) Exploratorio con CPM log2 para gráficos, SIN tests ======
cpm_mat <- apply(counts[,c("LOW","HIGH")], 2, function(x) 1e6 * x / sum(x))
log2cpm <- log2(cpm_mat + 1)
log2cpm_df <- cbind(geneID = counts$Geneid, as.data.frame(log2cpm))
log2cpm_long <- pivot_longer(as_tibble(log2cpm_df), cols = c(LOW, HIGH),
                             names_to = "sample", values_to = "expression")
ggplot(log2cpm_long, aes(sample, expression, fill = sample)) +
  geom_violin(show.legend=FALSE, trim=FALSE) +
  stat_summary(fun="median", geom="point", shape=95, size=10) +
  theme_bw() + labs(title="log2(CPM+1) (exploratorio, sin test)")

## ====== DEGseq (MARS) con 1 muestra por condición ======
# Formato de entrada de DEGexp: data.frame con columna de genes + columna de expresión
# Usaremos dos data.frames "gem1/gem2" con las mismas filas (mismos genes)
gem1 <- counts %>% select(Geneid, LOW)  %>% rename(gene = Geneid, exp = LOW)
gem2 <- counts %>% select(Geneid, HIGH) %>% rename(gene = Geneid, exp = HIGH)

# Ejecutar MARS (MA-plot-based method with Random Sampling)
# qValue=0.1 es típico en DEGseq; cambiá a 0.05 si querés más estricto
dir.create("DEGseq_cgigas", showWarnings = FALSE)

deg_out <- DEGexp(
  geneExpMatrix1 = gem1, geneCol1 = 1, expCol1 = 2, groupLabel1 = "LOW",
  geneExpMatrix2 = gem2, geneCol2 = 1, expCol2 = 2, groupLabel2 = "HIGH",
  method = "MARS",
  outputDir = "DEGseq_cgigas",
  qValue = 1,           # umbral de q-value
  thresholdKind = 5,      # usa qValue y fold-change como criterios
  foldChange = 0.1
)

# Ruta al archivo producido por DEGexp
res_path <- file.path("DEGseq_cgigas", "output_score.txt")

res <- read.delim(res_path, check.names = FALSE)

# Harmonizar nombres y añadir log2FC
res <- res %>%
  dplyr::rename(
    geneID = GeneNames,
    count_LOW = value1,
    count_HIGH = value2,
    log2FC = `log2(Fold_change)`,
    zscore = `z-score`,
    pvalue = `p-value`,
    q_BH = `q-value(Benjamini et al. 1995)`,
    q_Storey = `q-value(Storey et al. 2003)`
  )

res <- res %>% arrange(q_BH, desc(abs(log2FC)))


# Vista rápida de significativos por tu umbral (q<=0.1 & |log2FC|>=1)
sig <- res %>% filter(q_BH <= 0.001, abs(log2FC) >= 1)
head(sig)
# Guardar tabla ordenada
write.table(sig, file = "DEGseq_cgigas_results_stringent.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
## ====== Gráficos rápidos (MA y Volcano) ======
# MA-plot
ma_df <- counts %>%
  mutate(A = 0.5 * log2((LOW + 1) * (HIGH + 1)),
         M = log2((HIGH + 1) / (LOW + 1))) %>%
  left_join(res %>% select(geneID,q_BH), by = c("Geneid" = "geneID")) %>%
  mutate(sig = !is.na(q_BH) & q_BH <= 0.1)

ggplot(ma_df, aes(A, M, color = sig)) +
  geom_point(alpha=0.4, size=1) +
  scale_color_manual(values=c("FALSE"="grey70","TRUE"="red")) +
  geom_hline(yintercept = c(-1,1), linetype="dashed") +
  theme_bw() + labs(title="MA-plot (DEGseq MARS)", y="M=log2(HIGH/LOW)")

# Volcano
vol_df <- sig %>% mutate(sig = q_BH <= 0.1 & abs(log2FC) >= 1)
ggplot(vol_df, aes(log2FC, -log10(q_BH), color=sig)) +
  geom_point(alpha=0.5, size=1) +
  scale_color_manual(values=c("FALSE"="grey70","TRUE"="red")) +
  theme_bw() + labs(title="Volcano (DEGseq MARS)")


# Cargar tu tabla DIAMOND
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

# Combinar DEGs (sig) con el mapping LOC → NP/XP
deg_annot <- sig %>%
  left_join(loc2best, by = "geneID") %>%
  left_join(diamond_best, by = c("protein_id" = "qseqid")) %>%
  arrange(q_BH, desc(abs(log2FC)))

uniprot_info <- read_tsv("uniprotkb_bivalvia.tsv")

# Renombrar para estandarizar
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
write.table(annotated_DEGs, "cgigas_DEGs_with_Uniprot.tsv",
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

col_uniprot <- pick_col(c("uniprot_id"))
col_logfc   <- pick_col(c("log2FC"))
col_q       <- pick_col(c("q_BH"))
# Columnas posibles con GO IDs y/o nombres
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
             pvalueCutoff = 0.07,
             minGSSize  = 5,
             maxGSSize  = 500,
             verbose    = FALSE)

gsea_df <- as.data.frame(gsea)
write.table(gsea_df, "GSEA_GO_results.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# --- 6) (Opcional) Visualizaciones rápidas ---
if (nrow(gsea_df) > 0) {
  suppressPackageStartupMessages({library(enrichplot); library(ggplot2)})
  pdf("GSEA_GO_plots.pdf", width=8, height=6)
  print(dotplot(gsea, showCategory=20))
  print(ridgeplot(gsea, showCategory=20))
  dev.off()
}

# --- 7) Resumen de cobertura ---
message("Genes (UniProt) en ranking: ", length(geneVec))
message("Genes con al menos un GO:   ", length(intersect(names(geneVec), unique(TERM2GENE$uniprot_id))))
message("Términos GO en sets:       ", length(unique(TERM2GENE$go_id)))

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

