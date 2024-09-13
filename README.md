# mG4H-OFF-Bakshi-Wende-2024

**Title:** Transient Elevation in Cellular Glucose Uptake Exacerbates Pressure Overload-Induced Cardiac Hypertrophy and Dysfunction\
**Lead author:** Sayan Bakshi, M.Sc.\
**Corresponding author:** Dr. Adam R. Wende, Ph.D.\
**Place:** Department of Pathology, University of Alabama at Birmingham, Birmingham, Alabama, USA

# RNA-seq Analysis

## Libraries
```{R}
basic.packages <- c("dplyr", "tidyr", "openxlsx", "stringr", "ggplot2", "tibble")
viz.packages <- c("pheatmap", "RColorBrewer", "corrplot", "plotly")
gg.packages <- c("ggrepel", "ggpubr", "ggfortify", "gridExtra", "ggbeeswarm", "ggnewscale")
RNAseq.packages <- c(
  "biomaRt", "DESeq2", "vsn", "clusterProfiler", "enrichR",
  "fgsea", "pathview", "ReactomePA", "circlize", "enrichplot"
)
packages <- c(basic.packages, viz.packages, gg.packages, RNAseq.packages)
pacman::p_load(packages, character.only = TRUE)
```

## Input
Sample Groups: Con (Con-ON), Tg (mG4H-ON), DCon (Con-OFF), Toff (mG4H-OFF)
```{R}
## RNASeq Data Input
#### (Code modified from https://www.biostars.org/p/241602/)
raw.input <- list.files("../2_Input/Raw_Counts/", full.names = T)
counts.files <- lapply(raw.input, read.table, skip = 4) # Removing first 4 rows of summary data
counts <- as.data.frame(sapply(counts.files, function(x) x[, 4]))
raw.input %<>% gsub("../2_Input/Raw_Counts/", "", .) %>%
  gsub("_[A-Z]+_L[0-9]+[.]ReadsPerGene[.]out[.]tab", "", .)
colnames(counts) <- raw.input
row.names(counts) <- counts.files[[1]]$V1
counts <- counts[, order(colnames(counts))] %>%
  dplyr::rename(., "DCon2" = "Dcon2") # fixing Dcon2 naming error to DCon2 for consistency

## Sample Data Input
Sample.data <- read.xlsx("../2_Input/mG4H.OFF_Sample.Info.xlsx",
  sheet = "sample.info.v02"
) %>% dplyr::arrange(.$Name)

## Setting Factor Levels
sample.info <- Sample.data %>%
  dplyr::mutate(mG4H = factor(.$mG4, levels = c("not", "Tg"))) %>%
  dplyr::mutate(GT = factor(.$GT, levels = c("Con", "mG4H"))) %>%
  dplyr::mutate(Group_num = factor(.$Group_num, levels = c(
    "1", "2", "3", "4"
  ))) %>%
  dplyr::mutate(Group = factor(.$Group, levels = c(
    "Con.Veh", "mG4H.Veh",
    "Con.Veh.OFF", "mG4H.Veh.OFF"
  )))

## Result Annotation
mm39 <- useMart(
  biomart = "ENSEMBL_MART_ENSEMBL",
  dataset = "mmusculus_gene_ensembl"
)
bm <- getBM(attributes = c(
  "ensembl_gene_id", "external_gene_name", "ensembl_gene_id_version",
  "gene_biotype", "ensembl_transcript_id",
  "transcription_start_site", "transcript_length", "transcript_count",
  "refseq_mrna", "entrezgene_id", "strand", "chromosome_name",
  "start_position", "end_position", "description"
), mart = mm39)
```

## DESeq2 analysis
```{R}
## DESeq2 dataset preparation
cts <- as.matrix(counts[, as.factor(sample.info$Name)])
coldata <- sample.info %>%
  dplyr::select(c("Group", "mG4H", "GT", "Treat", "Group_num"))
rownames(coldata) <- sample.info$Name

## Differential Expression Analysis
dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = coldata,
  design = ~Group
)

## Prefiltering:
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

dds <- estimateSizeFactors(dds)
keep <- rowSums(counts(dds) >= 5) >= 3
dds <- dds[keep, ]

## DESeq2 analysis:
dds <- estimateDispersions(dds)
dds <- nbinomWaldTest(dds)

## Results Table:
res_TgvCon <- results(dds, contrast = c("Group", "mG4H.Veh", "Con.Veh"))
res_ToffvDCon <- results(dds, contrast = c("Group", "mG4H.Veh.OFF", "Con.Veh.OFF"))

## Extracting transformed values
rld <- rlog(dds, blind = F)
```

## Grouped PCA
```{R}
sample.PCA <- sample.info
mat <- assay(rld)

pca <- prcomp(t(mat))
Var <- pca$sdev^2 / sum(pca$sdev^2)
percentVar <- round(100 * Var, digits = 2)

pca <- as.data.frame(pca$x)
pca$Name <- rownames(pca)

sample.PCA.2 <- sample.PCA %>%
  dplyr::mutate(Group = paste(sample.PCA$GT, sample.PCA$Transgene, sep = "-"))

pca.comp <- left_join(pca, sample.PCA.2, by = "Name")
row.names(pca.comp) <- pca.comp$Name

ax_text <- list(family = "Arial", size = 25, color = "black")
tick_text <- list(family = "Arial", size = 15, color = "black")

p <- plot_ly(pca.comp,
  x = ~PC1, y = ~PC2, z = ~PC3,
  color = ~ pca.comp$Group,
  colors = c("#56b4e9", "black", "#1b9e77", "#d95f02"),
  marker = list(size = 20),
  text = paste(pca.comp$GT, pca.comp$Transgene, sep = "-")
) %>%
  add_markers() %>%
  layout(
    scene = list(
      xaxis = list(
        title = paste0("PC1: ", percentVar[1], "% variance"),
        zerolinewidth = 4, zerolinecolor = "darkgrey",
        linecolor = "darkgrey", linewidth = 4,
        titlefont = ax_text, tickfont = tick_text
      ),
      yaxis = list(
        title = paste0("PC2: ", percentVar[2], "% variance"),
        zerolinewidth = 4, zerolinecolor = "darkgrey",
        linecolor = "darkgrey", linewidth = 4,
        titlefont = ax_text, tickfont = tick_text
      ),
      zaxis = list(
        title = paste0("PC3: ", percentVar[3], "% variance"),
        zerolinewidth = 4, zerolinecolor = "darkgrey",
        linecolor = "darkgrey", linewidth = 4,
        titlefont = ax_text, tickfont = tick_text
      )
    ),
    annotations = list(
      x = 1.13,
      y = 1.03,
      text = "Group",
      xref = "1",
      yref = "0",
      showarrow = T,
      plot_bgcolor = "black"
    )
  )

p
```

## Downstream analysis
```{R}
## Annotation of normalized counts
Norm_counts <- counts(dds, normalized = T)
res$ensembl_gene_id <- rownames(res)

res_2 <- merge(res, Norm_counts, by = 0)
rownames(res_2) <- res_2$Row.names
res_2$ensembl_gene_id <- gsub("\\..*", "", res_2$ensembl_gene_id)
res_tmp.0 <- res_2[, c(2:8)]
res_tmp <- res_2[, c(9:ncol(res_2))]
res_select <- res_tmp[which(colnames(res_tmp) %in% phenodata$Name)]
res_select$ensembl_gene_id <- rownames(res_select)
res_select$ensembl_gene_id <- gsub("\\..*", "", res_select$ensembl_gene_id)
res_3 <- merge(res_tmp.0, res_select, by = "ensembl_gene_id")

res_annot <- merge(bm, res_3, by = c("ensembl_gene_id")) %>%
  mutate(across(everything(), ~ ifelse(. == "", NA, .)))
res_annot.2 <- res_annot %>% dplyr::select(-c(
  ensembl_transcript_id,
  transcription_start_site,
  transcript_length,
  transcript_count
))
res_clean <- res_annot.2[!duplicated(res_annot.2$ensembl_gene_id), ]
res_clean$log2FoldChange.ABS <- abs(res_clean$log2FoldChange)
res_clean$FoldChange.ABS <- 2^(res_clean$log2FoldChange.ABS)
res_clean$FoldChange <- ifelse(res_clean$log2FoldChange < 0, -(res_clean$FoldChange.ABS),
  res_clean$FoldChange.ABS
)

## Data Filtering
res_p05FC1.5 <- dplyr::filter(res_clean, abs(log2FoldChange) > log2(1.5) & pvalue < 0.05)
res_q0.1FC1.5 <- dplyr::filter(res_clean, abs(log2FoldChange) > log2(1.5) & padj < 0.1)
```

## Heatmap visualization
```{R}
hm_p05FC1.5 <- hm.df %>%
  dplyr::filter(pvalue < 0.05, abs(log2FoldChange) > log2(1.5)) %>%
  `rownames<-`(.$ensembl_gene_id) %>%
  dplyr::select(match) %>%
  dplyr::arrange(.) %>%
  as.matrix(.)

hm_q0.1FC1.5 <- hm.df %>%
  dplyr::filter(padj < 0.1, abs(log2FoldChange) > log2(1.5)) %>%
  `rownames<-`(.$ensembl_gene_id) %>%
  dplyr::select(match) %>%
  dplyr::arrange(.) %>%
  as.matrix(.)

ann_colors <- list(mycolor)
names(ann_colors) <- factor

callback <- function(hc, mat) {
  sv <- svd(t(mat))$v[, 1]
  dend <- reorder(as.dendrogram(hc), wts = sv)
  as.hclust(dend)
}

## p-value and absFC significant (p<0.05; |FC|>1.5)
hm_p05FC1.5.Plot <- pheatmap(hm_p05FC1.5,
  cluster_rows = T, show_rownames = FALSE,
  cluster_cols = T, annotation_colors = ann_colors,
  clustering_callback = callback,
  cutree_rows = 2, cutree_cols = 2,
  annotation_col = ann_col, scale = "row", fontsize = 7,
  main = paste0(comparison, "_", "p<0.05_|FC|>1.5", sep = ""),
  colorRampPalette(c("dodgerblue3", "white", "firebrick3"))(50)
)

## q-value and absFC significant (q<0.1; |FC|>1.5)
hm_q0.1FC1.5.Plot <- pheatmap(hm_q0.1FC1.5,
  cluster_rows = T, show_rownames = FALSE,
  cluster_cols = T, annotation_colors = ann_colors,
  clustering_callback = callback,
  cutree_rows = 2, cutree_cols = 2,
  annotation_col = ann_col, scale = "row", fontsize = 7,
  main = paste0(comparison, "_", "q<0.1_|FC|>1.5", sep = ""),
  colorRampPalette(c("dodgerblue3", "white", "firebrick3"))(50)
)
```

## Volcano plot visualization
```{R}
# Adapted from https://www.biostars.org/p/282295/

FCCutoff <- log2(1.5)
PCutoff <- 0.05
TopGenes <- 30

volcanotable <- res_clean %>%
  `rownames<-`(.$ensembl_gene_id) %>%
  dplyr::filter(!is.na(log2FoldChange) & !is.na(pvalue) & !is.na(padj))

volcanotable$Significance <- "NS"
volcanotable$Significance[(abs(volcanotable$log2FoldChange) > FCCutoff)] <- "FCSig"
volcanotable$Significance[(volcanotable$pvalue < PCutoff)] <- "PSig"
volcanotable$Significance[(volcanotable$pvalue < PCutoff) &
  (abs(volcanotable$log2FoldChange) > FCCutoff)] <- "FC_PSig"
table(volcanotable$Significance)
volcanotable$Significance <- factor(volcanotable$Significance,
  levels = c("NS", "FCSig", "PSig", "FC_PSig")
)
volcano.color <- c(NS = "black", FCSig = "#C0C0C0", PSig = "#0000FF", FC_PSig = "red2")
xmax <- max(abs(volcanotable$log2FoldChange), na.rm = T)
ymax <- max(-log10(volcanotable$pvalue), na.rm = T)

# Top PSig Genes
Ptop <- volcanotable %>%
  dplyr::select(ensembl_gene_id, external_gene_name, pvalue, log2FoldChange) %>%
  subset(pvalue < PCutoff & abs(log2FoldChange) > FCCutoff) %>%
  data.frame(.)

Ptop <- Ptop %>%
  arrange(pvalue) %>%
  `rownames<-`(.$ensembl_gene_id) %>%
  dplyr::slice(1:TopGenes)

Ptop$external_gene_name <- ifelse(is.na(Ptop$external_gene_name) == TRUE,
  Ptop$ensembl_gene_id, Ptop$external_gene_name
)
rownames(Ptop) <- Ptop$external_gene_name

# Top FCSig Genes
FCtop <- volcanotable %>%
  dplyr::select(ensembl_gene_id, external_gene_name, pvalue, log2FoldChange) %>%
  subset(pvalue < PCutoff & abs(log2FoldChange) > FCCutoff) %>%
  data.frame(.)

FCtop <- FCtop %>%
  dplyr::mutate(AbsFC = abs(FCtop$log2FoldChange)) %>%
  arrange(desc(AbsFC)) %>%
  `rownames<-`(.$ensembl_gene_id) %>%
  dplyr::slice(1:TopGenes) %>%
  dplyr::select(-AbsFC)

FCtop$external_gene_name <- ifelse(is.na(FCtop$external_gene_name) == TRUE,
  FCtop$ensembl_gene_id, FCtop$external_gene_name
)
rownames(FCtop) <- FCtop$external_gene_name

# Top PSig & FCSig Genes Combined
PFCtop <- rbind(Ptop, FCtop) %>%
  dplyr::filter(duplicated(ensembl_gene_id) == FALSE)

Volcano.PSig <- ggplot(volcanotable, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(fill = factor(Significance), color = factor(Significance)),
    colour = "black", shape = 21, size = 1
  ) +
  theme_bw(base_size = 24) +
  scale_fill_manual(values = volcano.color) +
  scale_color_manual(
    values = volcano.color,
    labels = c(
      NS = "NS", FCSig = paste("LogFC>|", FCCutoff, "|", sep = ""),
      PSig = paste("P-value P<", PCutoff, sep = ""),
      FC_PSig = paste("P-value P<", PCutoff,
        " & LogFC>|", FCCutoff, "|",
        sep = ""
      )
    )
  ) +
  theme(
    legend.background = element_rect(),
    plot.title = element_text(angle = 0, size = 12, face = "bold", vjust = 1),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 0, size = 12, vjust = 1),
    axis.text.y = element_text(angle = 0, size = 12, vjust = 1),
    axis.title = element_text(size = 12),
    legend.position = "top", legend.key = element_blank(),
    legend.key.size = unit(0.5, "cm"), legend.text = element_text(size = 8),
    title = element_text(size = 8), legend.title = element_blank()
  ) +
  guides(colour = guide_legend(override.aes = list(size = 2.5))) +
  scale_x_continuous(limits = c(-xmax, xmax)) +
  scale_y_continuous(limits = c(0, ymax)) +
  xlab(bquote(~ Log[2] ~ "fold change")) +
  ylab(bquote(~ -Log[10] ~ italic(P))) +
  ggtitle(paste0(comparison, "_", "Volcano Plot", sep = "")) +
  geom_vline(xintercept = c(-FCCutoff, FCCutoff), linetype = "longdash", colour = "black", size = 0.4) +
  geom_hline(yintercept = -log10(PCutoff), linetype = "longdash", colour = "black", size = 0.4)

Volcano.PSig_2 <- Volcano.PSig + geom_text_repel(
  data = PFCtop[c(1:nrow(PFCtop)), ],
  aes(label = rownames(PFCtop[c(1:nrow(PFCtop)), ])),
  size = 3, segment.size = 0.5,
  min.segment.length = 0, seed = 42, box.padding = 0.5
)

Volcano.PSig_2
```

## Pathway Enrichment Analysis (EnrichR)
```{R}
enrichrdbs <- listEnrichrDbs()
my_enrichrdbs <- c(
  "KEGG_2019_Mouse", "BioPlanet_2019", "GO_Biological_Process_2021",
  "Reactome_2016", "WikiPathways_2019_Mouse"
)

enrichr_sheet <- c(
  "KEGG_2019_Mouse", "BioPlanet_2019", "GO_Bio_Process_2021",
  "Reactome_2016", "WikiPath_2019_Mouse"
)

DEG_enriched <- enrichr(res_p05FC1.5$external_gene_name, my_enrichrdbs)
wb <- createWorkbook()
for (i in 1:length(enrichr_sheet)) {
  addWorksheet(wb = wb, sheetName = enrichr_sheet[i], gridLines = T)
  writeData(wb = wb, sheet = i, x = DEG_enriched[[my_enrichrdbs[i]]])
}
```

# RRBS Analysis
The R code and list of R packages repurposed for this analysis can be found in:
Article PMCID: PMC10292965
File name: Supplemental File S1
File link: https://doi.org/10.6084/m9.figshare.20419122

# Merging datasets
```{R}
## Merging both RNAseq datasets (selecting Q < 0.1 & |FC| > 1.5 DEG list)
DEG_TgvCon.ToffDCon_Merge <- merge(Tg.v.Con, Toff.v.DCon, by = "ensembl_gene_id")
DEG_Tg.v.Con.only <- anti_join(Tg.v.Con, Toff.v.DCon, by = "ensembl_gene_id")
DEG_Toff.v.DCon.only <- anti_join(Toff.v.DCon, Tg.v.Con, by = "ensembl_gene_id")

## Merging both Toff.v.DCon datasets - RNAseq & Methylseq
## (selecting Q < 0.1 & |FC| > 1.5 / |MD| > 5 DEG / DMC list)
ToffDCon_Merge <- inner_join(Toff.v.DCon_RNA, Toff.v.DCon_DNAm,
  by = c("external_gene_name" = "annot.symbol")
)
names(ToffDCon_Merge) <- gsub(names(ToffDCon_Merge),
  pattern = "\\.x", replacement = "\\.RNA"
) %>% gsub(names(ToffDCon_Merge), pattern = "\\.y", replacement = "\\.DNAm")

ToffDCon_Merge.df <- ToffDCon_Merge[!duplicated(ToffDCon_Merge$tag), ]
ToffDCon_RNA.only <- anti_join(Toff.v.DCon_RNA, Toff.v.DCon_DNAm,
  by = c("external_gene_name" = "annot.symbol")
)
ToffDCon_DNAm.only <- anti_join(Toff.v.DCon_DNAm, Toff.v.DCon_RNA,
  by = c("annot.symbol" = "external_gene_name")
)
```

# Session info
```{R}
sessionInfo()
```

R version 4.1.2 (2021-11-01)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 22621)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252   
[3] LC_MONETARY=English_United States.1252 LC_NUMERIC=C                          
[5] LC_TIME=English_United States.1252    

attached base packages:
[1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] org.Mm.eg.db_3.14.0         org.Hs.eg.db_3.14.0         AnnotationDbi_1.56.2       
 [4] enrichplot_1.14.2           circlize_0.4.15             ReactomePA_1.38.0          
 [7] pathview_1.34.0             fgsea_1.20.0                enrichR_3.0                
[10] clusterProfiler_4.2.2       vsn_3.62.0                  DESeq2_1.34.0              
[13] SummarizedExperiment_1.24.0 Biobase_2.54.0              MatrixGenerics_1.6.0       
[16] matrixStats_0.62.0          GenomicRanges_1.46.1        GenomeInfoDb_1.30.1        
[19] IRanges_2.28.0              S4Vectors_0.32.3            BiocGenerics_0.40.0        
[22] biomaRt_2.50.3              ggnewscale_0.4.7            ggbeeswarm_0.6.0           
[25] gridExtra_2.3               ggfortify_0.4.14            ggpubr_0.4.0               
[28] ggrepel_0.9.1               corrplot_0.92               RColorBrewer_1.1-3         
[31] pheatmap_1.0.12             tibble_3.1.7                ggplot2_3.3.6              
[34] stringr_1.4.0               openxlsx_4.2.5              tidyr_1.2.0                
[37] dplyr_1.0.9                

loaded via a namespace (and not attached):
  [1] pacman_0.5.1           utf8_1.2.2             R.utils_2.11.0         tidyselect_1.1.2      
  [5] RSQLite_2.2.14         grid_4.1.2             BiocParallel_1.28.3    scatterpie_0.1.7      
  [9] munsell_0.5.0          preprocessCore_1.56.0  withr_2.5.0            colorspace_2.0-3      
 [13] GOSemSim_2.20.0        filelock_1.0.2         knitr_1.39             rstudioapi_0.13       
 [17] ggsignif_0.6.3         DOSE_3.20.1            labeling_0.4.2         KEGGgraph_1.54.0      
 [21] GenomeInfoDbData_1.2.7 polyclip_1.10-0        bit64_4.0.5            farver_2.1.0          
 [25] downloader_0.4         vctrs_0.4.1            treeio_1.18.1          generics_0.1.2        
 [29] xfun_0.31              BiocFileCache_2.2.1    R6_2.5.1               graphlayouts_0.8.0    
 [33] locfit_1.5-9.5         bitops_1.0-7           cachem_1.0.6           gridGraphics_0.5-1    
 [37] DelayedArray_0.20.0    assertthat_0.2.1       scales_1.2.0           ggraph_2.0.5          
 [41] beeswarm_0.4.0         gtable_0.3.0           affy_1.72.0            tidygraph_1.2.1       
 [45] rlang_1.0.2            genefilter_1.76.0      GlobalOptions_0.1.2    splines_4.1.2         
 [49] rstatix_0.7.0          lazyeval_0.2.2         checkmate_2.1.0        broom_0.8.0           
 [53] BiocManager_1.30.18    yaml_2.3.5             reshape2_1.4.4         abind_1.4-5           
 [57] backports_1.4.1        qvalue_2.26.0          tools_4.1.2            ggplotify_0.1.0       
 [61] affyio_1.64.0          ellipsis_0.3.2         Rcpp_1.0.8.3           plyr_1.8.7            
 [65] progress_1.2.2         zlibbioc_1.40.0        purrr_0.3.4            RCurl_1.98-1.7        
 [69] prettyunits_1.1.1      viridis_0.6.2          magrittr_2.0.3         data.table_1.14.2     
 [73] DO.db_2.9              reactome.db_1.77.0     R.cache_0.15.0         hms_1.1.1             
 [77] patchwork_1.1.1        evaluate_0.15          xtable_1.8-4           XML_3.99-0.10         
 [81] shape_1.4.6            compiler_4.1.2         crayon_1.5.1           shadowtext_0.1.2      
 [85] R.oo_1.25.0            htmltools_0.5.2        ggfun_0.0.6            geneplotter_1.72.0    
 [89] aplot_0.1.6            DBI_1.1.3              tweenr_1.0.2           dbplyr_2.2.0          
 [93] MASS_7.3-54            rappdirs_0.3.3         Matrix_1.3-4           car_3.1-0             
 [97] cli_3.3.0              R.methodsS3_1.8.2      parallel_4.1.2         igraph_1.3.1          
[101] pkgconfig_2.0.3        xml2_1.3.3             ggtree_3.2.1           annotate_1.72.0       
[105] vipor_0.4.5            XVector_0.34.0         yulab.utils_0.0.4      digest_0.6.29         
[109] graph_1.72.0           Biostrings_2.62.0      rmarkdown_2.14         fastmatch_1.1-3       
[113] tidytree_0.3.9         curl_4.3.2             graphite_1.40.0        rjson_0.2.21          
[117] lifecycle_1.0.1        nlme_3.1-153           jsonlite_1.8.0         carData_3.0-5         
[121] viridisLite_0.4.0      limma_3.50.3           fansi_1.0.3            pillar_1.7.0          
[125] lattice_0.20-45        KEGGREST_1.34.0        fastmap_1.1.0          httr_1.4.3            
[129] survival_3.2-13        GO.db_3.14.0           glue_1.6.2             zip_2.2.0             
[133] png_0.1-7              bit_4.0.4              Rgraphviz_2.38.0       ggforce_0.3.3         
[137] stringi_1.7.6          blob_1.2.3             memoise_2.0.1          styler_1.8.1          
[141] ape_5.6-2             
