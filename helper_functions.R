root.dir <- "/Users/pengja/PennNeurosurgery Dropbox/Jackie Peng/SilverbushLab/Projects/Jackie_post-mortem/"
`%ni%` <- Negate(`%in%`)
library(Orthology.eg.db)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(tidyverse)
library(readxl)
library(AUCell)
library(ggplot2)
library(pals)
library(Seurat)
library(ggpubr)
library(rstatix)
library(cowplot)
library(patchwork)
library(ggrepel)
library(stringr)
library(scales)
library(FNN)

##### Plotting #####
correlation.method.use <- "spearman"

# Set a theme and adjust text sizes globally
theme_set(
  theme_bw() + 
    theme(
      axis.text = element_text(size = 16),      
      axis.title = element_text(size = 20, facet = "bold"),     
      plot.title = element_text(size = 24),
      plot.subtitle = element_text(size = 20),
      legend.title = element_text(size = 20), 
      legend.text = element_text(size = 16),
      strip.text = element_text(size = 20)
    )
)

# UMAP axes
umap.axes <- list(
  labs(x = "UMAP 1", y = "UMAP 2"),
  theme(
    axis.text  = element_blank(),
    axis.ticks = element_blank()
  )
)


patients <- c("C1027419", "C1037505", "C1060383", "C1061121", "C107625", "C15498",  "C176874",  "C2399853", "C2542041", "C2542533", "C2751264", "C34809",  "C547104",  "C70848",  "C714384",  "C799746")
patients.color.scale <- scale_color_manual(breaks=patients, values=as.character(polychrome(length(patients))))
patients.fill.scale <- scale_fill_manual(breaks=patients, values=as.character(polychrome(length(patients))))

timepoints <- c("Initial CNS Tumor", "Recurrence/Progressive (Non-Autopsy)", "Progressive (Autopsy)") 
timepoints.color.scale <- scale_color_manual(breaks=timepoints, values=c("#FFB000", "#FE6100", "#785EF0"))
timepoints.fill.scale <- scale_fill_manual(breaks=timepoints, values=c("#FFB000", "#FE6100", "#785EF0"))

timepoints.broad <- c("Surgical", "Postmortem")
timepoints.broad.colors <- c("#FE8900", "#785EF0")
timepoints.broad.color.scale <- scale_color_manual(breaks=timepoints.broad, values=timepoints.broad.colors)
timepoints.broad.fill.scale <- scale_fill_manual(breaks=timepoints.broad, values=timepoints.broad.colors)


# published sussman2024 cell types
cell.types <- c("Endothelial Cell", "Macrophage/Microglia", "Mature Neuron", "Mature Oligodendrocyte", "Neuroglial", "Pericyte", "T Cell")
cell.type.color.scale <- scale_color_manual(breaks=c(cell.types, NA), values=c(unname(trubetskoy(length(cell.types))), "#888888"))

# cell types re-annotated
cell.types.anotation_l1 <- c("Glial", "Immune", "Malignant cluster", "Neurons", "Vascular")
cell.types.anotation_l1.colors <- c(palette.colors(palette = "Okabe-Ito")[2:(length(cell.types.anotation_l1)+1)])
cell.types.anotation_l1.w.quality.colors <- c(cell.types.anotation_l1.colors, "#888888", "#888888")


cell.types.anotation_l1.color.scale <- scale_color_manual(breaks=cell.types.anotation_l1, values=cell.types.anotation_l1.w.quality.colors)

cell.types.anotation_l1.w.quality <- c(cell.types.anotation_l1, "Low quality", "Filtered out")
cell.types.anotation_l1.w.quality.color.scale <- scale_color_manual(breaks=cell.types.anotation_l1.w.quality, values=cell.types.anotation_l1.w.quality.colors)
cell.types.anotation_l1.w.quality.fill.scale <- scale_fill_manual(breaks=cell.types.anotation_l1.w.quality, values=cell.types.anotation_l1.w.quality.colors)
names(cell.types.anotation_l1.w.quality.colors) <- cell.types.anotation_l1.w.quality



cell.types.anotation_l2.w.quality.colors <- c(
  # Neurons (Teal/Cyan family)
  "Excitatory neurons"   = "#00868B",  # Dark Cyan
  "Inhibitory neurons"   = "#8EE5EE",  # Light Cyan
  
  # Glial (Green/Yellow family)
  "Astrocytes"           = "#228B22",  # Forest Green
  "Oligodendrocytes"     = "#A6D854",  # Light Apple Green
  "OPC"                  = "#F0E442",  # Bright Yellow
  
  # Immune (Blue/Indigo family)
  "Myeloid"              = "#084594",  # Deep Royal Blue
  "Lymphoid"             = "#377EB8",  # Medium Blue
  "Other immune"         = "#A6CEE3",  # Very Pale Blue
  
  # Vascular (Pink/Magenta family)
  "Endothelial"          = "#7A0177",  # Dark Magenta
  "Mural"                = "#CC79A7",  # Muted Rose
  "Vascular"             = "#CC79A7",  # Same as Mural
  
  # Malignant (The Anchor - Vermillion)
  "Malignant cluster"    = "#D55E00",  # High-visibility "Caution" color
  
  # QC
  "Low quality"          = "#EEEEEE",   # Very light "Whitesmoke" gray,
  "Filtered out"         = "#888888"
)


cell.types.anotation_l2.w.quality.colors.2 <- c(
  # Neurons (Teal/Cyan family)
  "Excitatory neurons"   = "#00868B",  # Dark Cyan
  "Inhibitory neurons"   = "#8EE5EE",  # Light Cyan
  
  # Glial (Green/Yellow family)
  "Astrocytes"           = "#228B22",  # Forest Green
  "Oligodendrocytes"     = "#A6D854",  # Light Apple Green
  "OPC"                  = "#F0E442",  # Bright Yellow
  
  # Immune (Blue/Indigo family)
  "Myeloid"              = "#084594",  # Deep Royal Blue
  "Lymphoid"             = "#377EB8",  # Medium Blue
  "Other immune"         = "#A6CEE3",  # Very Pale Blue
  
  # Vascular (Pink/Magenta family)
  "Endothelial"          = "#7A0177",  # Dark Magenta
  "Mural"                = "#CC79A7",  # Muted Rose
  "Vascular"             = "#CC79A7",  # Same as Mural
  
  # Malignant (The Anchor - Vermillion)
  "Malignant cluster"    = "grey", 
  
  # QC
  "Low quality"          = "#888888",  
  "Filtered out"         = "black"
)


cell.types.anotation_l2.w.quality <- names(cell.types.anotation_l2.w.quality.colors)
cell.types.anotation_l2.w.quality.color.scale <- scale_color_manual(breaks=cell.types.anotation_l2.w.quality, values=cell.types.anotation_l2.w.quality.colors)
cell.types.anotation_l2.w.quality.fill.scale <- scale_fill_manual(breaks=cell.types.anotation_l2.w.quality, values=cell.types.anotation_l2.w.quality.colors)
cell.types.anotation_l2.w.quality.fill.scale.2 <- scale_fill_manual(breaks=cell.types.anotation_l2.w.quality, values=cell.types.anotation_l2.w.quality.colors.2)


cell.types.anotation_l2 <- setdiff(cell.types.anotation_l2.w.quality, c("Low quality", "Filtered out"))
cell.types.anotation_l2.colors <- cell.types.anotation_l2.w.quality.colors[cell.types.anotation_l2]
cell.types.anotation_l2.color.scale <- scale_color_manual(breaks=cell.types.anotation_l2, values=cell.types.anotation_l2.colors)



# for our dataset
sample.ids <- c("P1 A3", "P1 A5", "P1 A8",  "P2 S1", "P2 S5", "P2 S5 10X", "P2 S5 Invent", "P2 S11", "P2 S8", "P2 S16", "P2 S16 10X", "P2 S16 Invent", "P4 A15", "P4 A35")
sample.ids.color.scale <- scale_color_manual(breaks = sample.ids, values = as.character(polychrome(length(sample.ids))))

patient.sections <- c("P1 A3",  "P1 A5", "P1 A8", "P2 S1", "P2 S5", "P2 S11", "P2 S8", "P2 S16", "P4 A15", "P4 A35")
patient.sections.color.scale <- scale_color_manual(breaks = patient.sections, values = as.character(polychrome(length(patient.sections))))


sample.descs.colors <- c(
  "GBM 1 Core"       = "#0B3C5D",
  "GBM 1 Edge"       = "#3288BD",
  "GBM 1 Distal"     = "#A6CEE3",
  
  "GBM 2 Core"       = "#8B0000",
  "GBM 2 Brainstem"  = "#49006A",
  "GBM 2 Edge"       = "#D7301F",
  "GBM 2 Distal Inf."  = "#7A0177",
  "GBM 2 Distal"     = "#F4A582",
  
  "IDH-mut 1 Edge"   = "#1B9E77",
  "IDH-mut 1 Distal" = "#A6DBA0"
)

sample.descs <- names(sample.descs.colors)
sample.descs.color.scale <- scale_color_manual(breaks = sample.descs, values = as.character(sample.descs.colors))
sample.descs.fill.scale <- scale_fill_manual(breaks = sample.descs, values = as.character(sample.descs.colors))


sample.descs.by.method.colors <- c(
  "GBM 1 Core FACS"       = "#0B3C5D",
  "GBM 1 Edge FACS"       = "#3288BD",
  "GBM 1 Distal FACS"     = "#A6CEE3",
  
  "GBM 2 Core FACS"       = "#8B0000",
  "GBM 2 Brainstem FACS"  = "#49006A",
  "GBM 2 Edge FACS"       = "#D7301F",
  "GBM 2 Edge Invent"       = "#D7301F",
  "GBM 2 Edge 10X"       = "#D7301F",
  "GBM 2 Distal Inf. FACS"  = "#7A0177",
  "GBM 2 Distal FACS"     = "#F4A582",
  "GBM 2 Distal Invent"     = "#F4A582",
  "GBM 2 Distal 10X"     = "#F4A582",
  
  "IDH-mut 1 Edge FACS"   = "#1B9E77",
  "IDH-mut 1 Distal FACS" = "#A6DBA0"
)

sample.descs.by.method <- names(sample.descs.by.method.colors)
sample.descs.by.method.color.scale <- scale_color_manual(breaks = sample.descs.by.method, values = as.character(sample.descs.by.method.colors))
sample.descs.by.method.fill.scale <- scale_fill_manual(breaks = sample.descs.by.method, values = as.character(sample.descs.by.method.colors))


sample.prep.method.colors <- c(
  "FACS"      = "#F8766D",
  "Invent"    = "#00BA38",
  "10X"       = "#619CFF"
)
sample.prep.methods <- names(sample.prep.method.colors)
sample.prep.method.color.scale <- scale_color_manual(breaks = sample.prep.methods, values = as.character(sample.prep.method.colors))


healthy.sections <- c("P1 A8", "P2 S16", "P2 S16 10X", "P2 S16 Invent", "P4 A35") # P2 S8 may have a few malignant cells


malignant.states.colors = c(
  "AC-like (Neftel)"          = "#2E7D32",
  "MES-like (Neftel)"         = "#C0392B",  
  "NPC-like (Neftel)"         = "#7B4FA3", 
  "OPC-like (Neftel)"         = "#4F81BD", 
  "undifferentiated (Tirosh)" = "#8C6D62", 
  "stem-like (Tirosh)"        = "#A23B72", 
  "AC-like (Tirosh)"          = "#7CB342", 
  "OC-like (Tirosh)"          = "#6BAED6",
  "Low quality"               = "#EEEEEE"
)

#neftel.colors.hexcode <- malignant.states.colors[1:4]
#tirosh.colors.hexcode <- malignant.states.colors[5:8]
malignant.states <- names(malignant.states.colors)
malignant.state.color.scale <- scale_color_manual(breaks = names(malignant.states.colors), values = as.character(malignant.states.colors))
malignant.state.fill.scale <- scale_fill_manual(breaks = names(malignant.states.colors), values = as.character(malignant.states.colors))


cell.types.l2.w.neftel <- c(cell.types.anotation_l2, malignant.states)
cell.types.l2.w.neftel.color.scale <- scale_color_manual(breaks = cell.types.l2.w.neftel, values = c(cell.types.anotation_l2.colors, malignant.states.colors))
cell.types.l2.w.neftel.fill.scale <- scale_fill_manual(breaks = cell.types.l2.w.neftel, values = c(cell.types.anotation_l2.colors, malignant.states.colors))


sample.location.colors = c(
  "Core"         = "#081D40",
  "Brainstem"    = "#2171B5", 
  "Edge"         = "#41B6C4", 
  "Distal Inf."  = "#C7E9B4",
  "Distal"      = "#E6F7A0"
)
sample.locations <- names(sample.location.colors)
sample.location.color.scale <- scale_color_manual(breaks = sample.locations, values = as.character(sample.location.colors))
sample.location.fill.scale <- scale_fill_manual(breaks = sample.locations, values = as.character(sample.location.colors))



###### Quality thresholds #####
percent.mt.cutoff <- 10
percent.mt.cutoff.filter.out <- 20
nCount_RNA.cutoff.filter.out <- 500
nCount_RNA.cutoff.low <- 1000
nCount_RNA.cutoff.high <- 100000
nFeature_RNA.cutoff.low <- 500
nFeature_RNA.cutoff.high <- 10000
post.mortem.noise.cutoff <- 0.8






###### Gene signatures #####

#' function for reading in GSEA gene sets
#' 
#' reads in a .grp file from GSEA and returns a vector of genes in the gene set
#'
#' @param filepath the path to the .grp file to read in
read.in.gene.set <- function(filepath) {
  gene.set <- read.csv(filepath, skip=2, header = F) %>% pull(V1)
  return(gene.set)
}


# GSEA signatures
gsea.dir <- paste0(root.dir, "ExternalData/", "GSEA_gene_sets/")
gene.signatures.gsea <- list()
gene.signatures.gsea[["apoptosis.signature"]] <- read.in.gene.set(paste0(gsea.dir, "HALLMARK_APOPTOSIS.v2025.1.Hs.grp"))
gene.signatures.gsea[["necrosis.signature"]] <- read.in.gene.set(paste0(gsea.dir, "REACTOME_REGULATED_NECROSIS.v2025.1.Hs.grp"))
gene.signatures.gsea[["autophagy.signature"]] <- read.in.gene.set(paste0(gsea.dir, "REACTOME_AUTOPHAGY.v2025.1.Hs.grp"))
gene.signatures.gsea[["hypoxia.signature"]] <- read.in.gene.set(paste0(gsea.dir, "HALLMARK_HYPOXIA.v2025.1.Hs.grp"))
gene.signatures.gsea[["inflammatory.signature"]] <- read.in.gene.set(paste0(gsea.dir, "HALLMARK_INFLAMMATORY_RESPONSE.v2025.1.Hs.grp"))
gene.signatures.gsea[["interferon.alpha.signature"]] <- read.in.gene.set(paste0(gsea.dir, "HALLMARK_INTERFERON_ALPHA_RESPONSE.v2025.1.Hs.grp"))
gene.signatures.gsea[["interferon.gamma.signature"]] <- read.in.gene.set(paste0(gsea.dir, "HALLMARK_INTERFERON_GAMMA_RESPONSE.v2025.1.Hs.grp"))
gene.signatures.gsea[["interferon.signature"]] <- read.in.gene.set(paste0(gsea.dir, "REACTOME_INTERFERON_SIGNALING.v2025.1.Hs.grp"))
gene.signatures.gsea[["g2m.signature"]] <- read.in.gene.set(paste0(gsea.dir, "HALLMARK_G2M_CHECKPOINT.v2026.1.Hs.grp"))



gene.signatures.gsea.celltype <- list()


gene.signatures.gsea.celltype[["CELLULAR_RESPONSE_TO_STRESS"]] <- read.in.gene.set(paste0(gsea.dir, "GOBP_CELLULAR_RESPONSE_TO_STRESS.v2026.1.Hs.grp"))
gene.signatures.gsea.celltype[["CELLULAR_RESPONSE_TO_HEAT"]] <- read.in.gene.set(paste0(gsea.dir, "cell_type/", "GOBP_CELLULAR_RESPONSE_TO_HEAT.v2026.1.Hs.grp"))
gene.signatures.gsea.celltype[["PHAGOCYTOSIS"]] <- read.in.gene.set(paste0(gsea.dir, "cell_type/", "GOBP_PHAGOCYTOSIS.v2026.1.Hs.grp"))
gene.signatures.gsea.celltype[["T_CELL_DIFFERENTIATION"]] <- read.in.gene.set(paste0(gsea.dir, "cell_type/", "GOBP_T_CELL_DIFFERENTIATION.v2026.1.Hs.grp"))
gene.signatures.gsea.celltype[["VASCULAR_PROCESS"]] <- read.in.gene.set(paste0(gsea.dir, "cell_type/", "GOBP_VASCULAR_PROCESS_IN_CIRCULATORY_SYSTEM.v2026.1.Hs.grp"))
gene.signatures.gsea.celltype[["EXTRACELLULAR_MATRIX_ORGANIZATION"]] <- read.in.gene.set(paste0(gsea.dir, "cell_type/", "GOBP_POSITIVE_REGULATION_OF_EXTRACELLULAR_MATRIX_ORGANIZATION.v2026.1.Hs.grp"))
gene.signatures.gsea.celltype[["HYPOXIA"]] <- read.in.gene.set(paste0(gsea.dir, "HALLMARK_HYPOXIA.v2025.1.Hs.grp"))
gene.signatures.gsea.celltype[["LEUKOCYTE_MIGRATION"]] <- read.in.gene.set(paste0(gsea.dir, "cell_type/", "GOBP_LEUKOCYTE_MIGRATION.v2026.1.Hs.grp"))
gene.signatures.gsea.celltype[["VASCULAR_TRANSPORT"]] <- read.in.gene.set(paste0(gsea.dir, "cell_type/", "GOBP_VASCULAR_TRANSPORT.v2026.1.Hs.grp"))
gene.signatures.gsea.celltype[["T_CELL_ACTIVATION"]] <- read.in.gene.set(paste0(gsea.dir, "cell_type/", "GOBP_T_CELL_ACTIVATION.v2026.1.Hs.grp"))
gene.signatures.gsea.celltype[["GOBP_WOUND_HEALING"]] <- read.in.gene.set(paste0(gsea.dir, "cell_type/", "GOBP_WOUND_HEALING.v2026.1.Hs.grp"))
gene.signatures.gsea.celltype[["ANGIOGENESIS"]] <- read.in.gene.set(paste0(gsea.dir, "cell_type/", "HALLMARK_ANGIOGENESIS.v2026.1.Hs.grp"))
gene.signatures.gsea.celltype[["EPITHELIAL_MESENCHYMAL_TRANSITION"]] <- read.in.gene.set(paste0(gsea.dir, "cell_type/", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.v2026.1.Hs.grp"))






# FGFR signaling
fgfr.gene.signatures <- list()

fgfr.gene.signatures[["HALLMARK_PI3K_AKT_MTOR_SIGNALING"]] <- read.in.gene.set(paste0(gsea.dir, "FGFR/", "HALLMARK_PI3K_AKT_MTOR_SIGNALING.v2025.1.Hs.grp"))
fgfr.gene.signatures[["KEGG_MAPK_SIGNALING_PATHWAY"]] <- read.in.gene.set(paste0(gsea.dir, "FGFR/", "KEGG_MAPK_SIGNALING_PATHWAY.v2026.1.Hs.grp"))
fgfr.gene.signatures[["HALLMARK_KRAS_SIGNALING_UP"]] <- read.in.gene.set(paste0(gsea.dir, "FGFR/", "HALLMARK_KRAS_SIGNALING_UP.v2026.1.Hs.grp"))
fgfr.gene.signatures[["HALLMARK_MYC_TARGETS"]] <- read.in.gene.set(paste0(gsea.dir, "FGFR/", "HALLMARK_MYC_TARGETS_V1.v2026.1.Hs.grp"))
fgfr.gene.signatures[["HALLMARK_E2F_TARGETS"]] <- read.in.gene.set(paste0(gsea.dir, "FGFR/", "HALLMARK_E2F_TARGETS.v2026.1.Hs.grp"))
fgfr.gene.signatures[["GOBP_PHOSPHOLIPASE_C"]] <- read.in.gene.set(paste0(gsea.dir, "FGFR/", "GOBP_PHOSPHOLIPASE_C_ACTIVATING_G_PROTEIN_COUPLED_RECEPTOR_SIGNALING_PATHWAY.v2026.1.Hs.grp"))
fgfr.gene.signatures[["BIOCARTA_PLC_PATHWAY"]] <- read.in.gene.set(paste0(gsea.dir, "FGFR/", "BIOCARTA_PLC_PATHWAY.v2026.1.Hs.grp"))
fgfr.gene.signatures[["HALLMARK_IL6_JAK_STAT3_SIGNALING"]] <- read.in.gene.set(paste0(gsea.dir, "FGFR/", "HALLMARK_IL6_JAK_STAT3_SIGNALING.v2026.1.Hs.grp"))
fgfr.gene.signatures[["HALLMARK_G2M_CHECKPOINT"]] <- read.in.gene.set(paste0(gsea.dir, "HALLMARK_G2M_CHECKPOINT.v2026.1.Hs.grp"))


 
# fgfr.gene.signatures[["REACTOME_FGFR2_LIGAND_BINDING_AND_ACTIVATION"]] <- read.in.gene.set(paste0(gsea.dir, "FGFR/", "REACTOME_FGFR2_LIGAND_BINDING_AND_ACTIVATION.v2025.1.Hs.grp"))
# fgfr.gene.signatures[["REACTOME_SIGNALING_BY_FGFR2"]] <- read.in.gene.set(paste0(gsea.dir, "FGFR/", "REACTOME_SIGNALING_BY_FGFR2.v2025.1.Hs.grp"))
# #fgfr.gene.signatures[["REACTOME_DOWNSTREAM_SIGNALING_OF_ACTIVATED_FGFR2"]] <- read.in.gene.set(paste0(gsea.dir, "FGFR/", "REACTOME_DOWNSTREAM_SIGNALING_OF_ACTIVATED_FGFR2.v2025.1.Hs.grp"))
# #fgfr.gene.signatures[["REACTOME_SIGNALING_BY_FGFR2_IN_DISEASE"]] <- read.in.gene.set(paste0(gsea.dir, "FGFR/", "REACTOME_SIGNALING_BY_FGFR2_IN_DISEASE.v2025.1.Hs.grp"))
# 
# fgfr.gene.signatures[["REACTOME_FGFR1_LIGAND_BINDING_AND_ACTIVATION"]] <- read.in.gene.set(paste0(gsea.dir, "FGFR/", "REACTOME_FGFR1_LIGAND_BINDING_AND_ACTIVATION.v2025.1.Hs.grp"))
# fgfr.gene.signatures[["REACTOME_SIGNALING_BY_FGFR1"]] <- read.in.gene.set(paste0(gsea.dir, "FGFR/", "REACTOME_SIGNALING_BY_FGFR1.v2025.1.Hs.grp"))
# 
# fgfr.gene.signatures[["REACTOME_SIGNALING_BY_FGFR"]] <- read.in.gene.set(paste0(gsea.dir, "FGFR/", "REACTOME_SIGNALING_BY_FGFR.v2025.1.Hs.grp"))
# 
# fgfr.gene.signatures[["REACTOME_MAPK_TARGETS_NUCLEAR_EVENTS_MEDIATED_BY_MAP_KINASES"]] <- read.in.gene.set(paste0(gsea.dir, "FGFR/", "REACTOME_MAPK_TARGETS_NUCLEAR_EVENTS_MEDIATED_BY_MAP_KINASES.v2025.1.Hs.grp"))
# 
# fgfr.gene.signatures[["KEGG_MEDICUS_REFERENCE_FGF_FGFR_RAS_ERK_SIGNALING_PATHWAY"]] <- read.in.gene.set(paste0(gsea.dir, "FGFR/", "KEGG_MEDICUS_REFERENCE_FGF_FGFR_RAS_ERK_SIGNALING_PATHWAY.v2025.1.Hs.grp"))
# fgfr.gene.signatures[["KEGG_MEDICUS_REFERENCE_FGF_FGFR_PI3K_SIGNALING_PATHWAY"]] <- read.in.gene.set(paste0(gsea.dir, "FGFR/", "KEGG_MEDICUS_REFERENCE_FGF_FGFR_PI3K_SIGNALING_PATHWAY.v2025.1.Hs.grp"))


# Brain cell type markers McKenzie 2018
celltype.markers.df <- read_excel(paste0(root.dir, "ExternalData/", "McKenzie_2018_celltype_markers/", "41598_2018_27293_MOESM2_ESM.xlsx"), sheet = 1, skip = 2)
top.n.celltype.markers <- 100

brain_signatures <- list()
brain_signatures_all <- list() # no cutoff
for (x in unique(celltype.markers.df$Celltype)) {
  celltype <- case_when(
    x == "ast" ~ "Astrocytes",
    x == "end" ~ "Vascular",
    x == "mic" ~ "Myeloid",
    x == "neu" ~ "Neurons",
    x == "oli" ~ "Oligodendrocytes",
    x == "opc" ~ "OPC",
  )
  top.genes <- celltype.markers.df %>% filter(Celltype == x) %>% arrange(desc(grand_mean)) %>% head(top.n.celltype.markers) %>% pull(gene)
  brain_signatures[[celltype]] <- top.genes
  brain_signatures_all[[celltype]] <- celltype.markers.df %>% filter(Celltype == x) %>% arrange(desc(grand_mean)) %>% pull(gene)
}


# GBM signatures
signatures.dir="/Users/pengja/PennNeurosurgery Dropbox/Jackie Peng/SilverbushLab/Shared_resources/Signatures/Neftel_GBM_cell_states/SingleCell_GBM/"

ac_sig <- readLines(paste0(signatures.dir, 'ac_genes.txt'))
mes1_sig <- readLines(paste0(signatures.dir, 'mes1_genes.txt'))
mes2_sig <- readLines(paste0(signatures.dir, 'mes2_genes.txt'))
npc1_sig <- readLines(paste0(signatures.dir, 'npc1_genes.txt'))
npc2_sig <- readLines(paste0(signatures.dir, 'npc2_genes.txt'))
opc_sig <- readLines(paste0(signatures.dir, 'opc_genes.txt'))

gbm.signatures <- list()
gbm.signatures[["AC.like"]] <- ac_sig
gbm.signatures[["MES.like"]] <- c(mes1_sig, mes2_sig)
gbm.signatures[["NPC.like"]] <- c(npc1_sig, npc2_sig)
gbm.signatures[["OPC.like"]] <- opc_sig

celltype.signatures <- c(brain_signatures, gbm.signatures)



# Demyelination and remyelination signature

# convert mouse to human genes
# x = list of gene symbols in mouse
mouse.to.human <- function(x) {
  mouse <- mapIds(org.Mm.eg.db, x, "ENTREZID", "ALIAS")
  mapped <- AnnotationDbi::select(Orthology.eg.db, mouse,"Homo.sapiens", "Mus.musculus")
  human <- mapIds(org.Hs.eg.db, as.character(mapped[,2]), "SYMBOL","ENTREZID")
  human <- do.call(c, lapply(human, function(x) if(is.null(x)) return(NA) else return(x)))
  cbind(x, mapped, human)
}

read.in.myelin.signature <- function(myelination.df, start.col, end.col, log2fc.cutoff = 0.5, pval.cutoff = 0.05, top.n.genes = 100) {
  
  myelination.df <- myelination.df[,start.col:end.col]
  colnames(myelination.df) <- c("gene", "base.mean", "log2fc", "pvalue", "padj")
  myelination.df$pvalue <- as.numeric(myelination.df$pvalue)
  myelination.df$padj <- as.numeric(myelination.df$padj)
  mouse.genes <- myelination.df %>% filter(log2fc > log2fc.cutoff, padj < pval.cutoff) %>% arrange(padj) %>% head(top.n.genes) %>% pull(gene)
  
  human.genes <- mouse.to.human(mouse.genes)$human
  human.genes <- human.genes[!is.na(human.genes)]
  
  return(human.genes)
}

# acute myelination
acute.myelination.df <- read_excel(paste0(root.dir, "ExternalData/", "Samtani_2023_de_remyelination/", "Table_1.XLSX"), sheet = 2, skip = 1)

# chronic myelination
chronic.myelination.df <- read_excel(paste0(root.dir, "ExternalData/", "Samtani_2023_de_remyelination/", "Table_2.XLSX"), sheet = 1, skip = 0)


human.genes.acute.demyelination <- read.in.myelin.signature(acute.myelination.df, 1, 5)
human.genes.acute.remyelination <- read.in.myelin.signature(acute.myelination.df, 13, 17)
human.genes.acute.de.re.myelination <- read.in.myelin.signature(acute.myelination.df, 7, 11)

human.genes.chronic.demyelination <- read.in.myelin.signature(chronic.myelination.df, 1, 5)
human.genes.chronic.remyelination <- read.in.myelin.signature(chronic.myelination.df, 25, 29)
human.genes.chronic.de.re.myelination <- read.in.myelin.signature(chronic.myelination.df, 7, 11)

myelination.signatures <- list()
myelination.signatures[["acute.demyelination"]] <- human.genes.acute.demyelination
myelination.signatures[["acute.remyelination"]] <- human.genes.acute.remyelination
myelination.signatures[["acute.de.re.myelination"]] <- human.genes.acute.de.re.myelination
myelination.signatures[["chronic.demyelination"]] <- human.genes.chronic.demyelination
myelination.signatures[["chronic.remyelination"]] <- human.genes.chronic.remyelination
myelination.signatures[["chronic.de.re.myelination"]] <- human.genes.chronic.de.re.myelination



##### Reactive oligodendrocyte signatures #####
#markers.dir <- paste0(root.dir, "Results/", "my_data/", "autopsy/", "all_postmortem_samples_cellbender.FALSE_denoise.heat.shock.signature_postmortem.sensitive_postmortem.signature_filter.FALSE_batch.correction.None/", "markers/", "Oligodendrocytes", "/")
markers.dir <- paste0(root.dir, "Figures/Figure 4/markers/Oligodendrocytes/")


# oligo vs tumor
#oligo.tumor.markers.summary.df <- read.csv(paste0(markers.dir, "my_data_", "Oligodendrocytes", ".vs.malignant_summary_use.latent.vars=", "heat.shock.signature_postmortem.sensitive_postmortem.signature", "_",  "MAST", ".tsv"), sep = "\t")
#oligo.vs.tumor.degs <- oligo.tumor.markers.summary.df %>% filter(n.patients >= 2, n.tumor.samples >= 5) %>% pull(gene)

oligo.tumor.markers.summary.df <- read.csv(paste0(markers.dir, "summary_DEGs_matched_Oligodendrocytes_vs_malignant_denoise=TRUE.tsv"), sep = "\t")
oligo.vs.tumor.degs <- oligo.tumor.markers.summary.df %>% filter(n.patients >= 2, n.samples >= 5) %>% pull(gene)

# reactive oligo markers
#summary.df <- read.csv(paste0(markers.dir, "my_data_all.markers.tumor.vs.normal_summary_use.latent.vars=", "heat.shock.signature_postmortem.sensitive_postmortem.signature", "_", "MAST", "_downsample=", "FALSE", ".tsv"), sep = "\t")
#reactive.oligo.filter.df <- summary.df %>% filter(n.patients >= 2, n.tumor.samples >= 4) %>% arrange(median.p_val_adj)

summary.df <- read.csv(paste0(markers.dir, "summary_DEGs_matched_sample_vs_distal_denoise=TRUE.tsv"), sep = "\t")
reactive.oligo.filter.df <- summary.df %>% filter(n.patients >= 2, n.samples >= 4) %>% arrange(median.p_val_adj)

#reactive.oligo.specific.filter.df <- reactive.oligo.filter.df %>% filter(gene %in% oligo.vs.rest.degs)
reactive.oligo.specific.filter.df <- reactive.oligo.filter.df %>% filter(gene %in% brain_signatures_all$Oligodendrocytes)
reactive.oligo.non.tumor.filter.df <- reactive.oligo.filter.df %>% filter(gene %in% oligo.vs.tumor.degs)


# Oligo signatures
oligo.gene.signatures <- list()
oligo.gene.signatures[["myelin.genes"]] <- c("PLP1", "MBP", "MOBP", "MOG", "MAG", "CNP")
oligo.gene.signatures[["reactive.oligo.general"]] <- reactive.oligo.filter.df %>% dplyr::filter(median.avg_log2FC > 0.5) %>% head(100) %>% pull(gene)
oligo.gene.signatures[["reactive.oligo.non.tumor"]] <- reactive.oligo.non.tumor.filter.df$gene
oligo.gene.signatures[["reactive.oligo.specific"]] <- reactive.oligo.specific.filter.df$gene
oligo.gene.signatures[["oligo.like"]]  <- oligo.gene.signatures[["reactive.oligo.general"]] 


# Jakel 2019 signatures
oligo.jakel2019.signatures <- list()
oligo.jakel2019.signatures[["OPCs"]] <- c("SOX6", "VCAN", "SNTG1", "CA10", "MEGF11", "HIF3A", "PCDH15", "LHFPL3", "MEGF11", "NXPH1", "PTPRZ1", "ITPR2", "PDGFRA")
oligo.jakel2019.signatures[["COPs"]] <- c("DYSF", "PLEKHG3", "RYR2", "CHRM5", "LINC01608", "PLP1", "CNP", "ST18", "MAG", "MOBP")
oligo.jakel2019.signatures[["imOLG"]] <- c("ARHGAP24", "MEF2C", "LRMDA", "APOE", "CD74", "DOCK8", "PLXDC2", "ELL2", "APBB1IP", "HLA-DRA", "C3", "PTPRC") # original: c("ARHGAP24", "MEF2C", "C10orf11", "APOE", "CD74", "DOCK8", "PLXDC2", "ELL2", "APBB1IP", "HLA.DRA", "C3", "PTPRC")
oligo.jakel2019.signatures[["Oligo1"]] <- c("MYO1D", "FAM13C", "PDE1A", "RBFOX1", "RASGRF1", "PLEKHG1", "FBXO32", "AFF3")
oligo.jakel2019.signatures[["Oligo2"]] <- c("HSPA1A", "CDH19")
oligo.jakel2019.signatures[["Oligo3"]] <- c("TPPP")
oligo.jakel2019.signatures[["Oligo4"]] <- c("LINC01608")
oligo.jakel2019.signatures[["Oligo5"]] <- c("TNFRSF21", "GJB1", "TUBB2B", "KLK6")
oligo.jakel2019.signatures[["Oligo6"]] <- c("OPALIN", "LINC00844", "COL4A5", "EFHD1", "LIPA", "LPAR1", "CNDP1")


# Mikolajewicz signatures
mikolajewicz.signatures.df <- read.csv(paste0(root.dir, "ExternalData/", "Mikolajewicz_2026_GBM_OL/", "supplemental/", "mmc4.csv"))
mikolajewicz.signatures <- list()
for (i in 1:length(colnames(mikolajewicz.signatures.df))) {
  genes <- mikolajewicz.signatures.df[,i]
  genes <- genes[genes != ""]
  mikolajewicz.signatures[[colnames(mikolajewicz.signatures.df)[i]]] <- genes
}


# From Brad Zuchero
fgf.gene.signatures <- list()
df <- read.csv(paste0(root.dir, "ExternalData/", "Weigel_Cooper_2026_FGF/", "__top100_all_markers_per_cluster_oligo.csv"))

for (cluster.i in unique(df$cluster)) {
  mouse.genes <- df %>% filter(cluster == cluster.i) %>% pull(gene)
  human.genes <- mouse.to.human(mouse.genes)$human
  human.genes <- human.genes[!is.na(human.genes)]
  fgf.gene.signatures[[cluster.i]] <- human.genes
}


##### Postmortem signatures #####
mt.genes <- c("MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND4L", "MT-ND5", "MT-ND6", "MT-CYB", "MT-CO1", "MT-CO2", "MT-CO3", "MT-ATP6", "MT-ATP8")

signatures.rank.df <- read.csv(paste0(root.dir, "Results/tensor_decomposition/Sussman_2024_pHGG/l1_top_10000_vf/non_negative_parafac_5000_iter_tolerance_e_32_n=100_1_to_15_seed_777/3 runs combined/lof 0.05 20 neighbors/rank_8", "/gene_signatures.tsv"), sep = "\t")

gene.signatures <- list()
gene.signatures[["postmortem.glial"]] <- signatures.rank.df[,6][1:100]
gene.signatures[["postmortem.signature"]] <- signatures.rank.df[,7][1:100]
gene.signatures[["postmortem.sensitive"]] <- signatures.rank.df[,8][1:100]
gene.signatures[["heat.shock.signature"]] <- c("HSPA1A", "HSPA1B", "HSPA6", "HSPH1", "HSPB1") 

# core hsp: c("HSPA6", "HSPA1A", "HSPA1B", "HSPA4L", "HSPH1", "HSPD1", "HSPB1", "HSPB8", "HSPE1")
#all heat shock 70 family genes c("HSPA1A", "HSPA1B", "HSPA1L", "HSPA2", "HSPA5", "HSPA6", "HSPA8", "HSPA9", "HSPA12A", "HSPA12B", "HSPA13", "HSPA14")

celltype.genes <- unlist(c(brain_signatures_all, gbm.signatures)) %>% as.character()
gene.signatures$postmortem.signature <- setdiff(gene.signatures$postmortem.signature, c(celltype.genes))
gene.signatures$postmortem.sensitive <- setdiff(gene.signatures$postmortem.sensitive, c(celltype.genes))
#gene.signatures$postmortem.glial <- setdiff(gene.signatures$postmortem.glial, c(brain_signatures_all$Neurons)) 

gene.signatures[["malignant.factor"]] <- signatures.rank.df[,1][1:100]
gene.signatures[["neurons.factor"]] <- signatures.rank.df[,2][1:100]
gene.signatures[["vascular.factor"]] <- signatures.rank.df[,3][1:100]
gene.signatures[["glial.factor"]] <- signatures.rank.df[,4][1:100]
gene.signatures[["immune.factor"]] <- signatures.rank.df[,5][1:100]

gene.signatures[["oligo.like"]] <- oligo.gene.signatures[["reactive.oligo.general"]] 

print(gene.signatures)


#old:
#postmortem.signatures.df <- read.csv(paste0(root.dir, "Results/", "Signatures/", "postmortem.signatures.tsv"), sep = "\t")

all.gene.signatures <- gene.signatures
all.gene.signatures[["postmortem.glial.factor.100"]] <- signatures.rank.df[,6][1:100]
all.gene.signatures[["postmortem.glial.factor.50"]] <- signatures.rank.df[,6][1:50]
all.gene.signatures[["postmortem.signature.factor"]] <- signatures.rank.df[,7][1:100]
all.gene.signatures[["postmortem.sensitive.factor"]] <- signatures.rank.df[,8][1:100]

signatures.rank.9.df <- read.csv(paste0(root.dir, "Results/tensor_decomposition/Sussman_2024_pHGG/l1_top_10000_vf/non_negative_parafac_5000_iter_tolerance_e_32_n=100_1_to_15_seed_777/3 runs combined/lof 0.05 20 neighbors/rank_9", "/gene_signatures.tsv"), sep = "\t")
signatures.rank.10.df <- read.csv(paste0(root.dir, "Results/tensor_decomposition/Sussman_2024_pHGG/l1_top_10000_vf/non_negative_parafac_5000_iter_tolerance_e_32_n=100_1_to_15_seed_777/3 runs combined/lof 0.05 20 neighbors/rank_10", "/gene_signatures.tsv"), sep = "\t")

all.gene.signatures[["heat.shock.signature.9"]] <- signatures.rank.9.df[,8][1:100]
all.gene.signatures[["heat.shock.signature.10"]] <- signatures.rank.10.df[,8][1:100]
all.gene.signatures[["heat.shock.signature.overlap.factor"]] <- intersect(all.gene.signatures[["heat.shock.signature.9"]], all.gene.signatures[["heat.shock.signature.10"]])

liharska.2025.postmortem.df <- read_excel(paste0(root.dir, "ExternalData/", "Liharska_2025_Postmortem/", "41380_2025_3163_MOESM8_ESM.xlsx"), sheet = 2)
gene.mapping <- mapIds(org.Hs.eg.db, keys = liharska.2025.postmortem.df$gene, keytype = "ENSEMBL", column="SYMBOL")
gene.mapping.df <- data.frame(symbol = gene.mapping, gene = names(gene.mapping))

liharska.2025.postmortem.df <- liharska.2025.postmortem.df %>% inner_join(gene.mapping.df, by = "gene")

liharska.2025.postmortem.df.filter <- liharska.2025.postmortem.df %>% filter(!is.na(symbol)) # only use genes with valid gene symbol
liharska.2025.liv.degs <- liharska.2025.postmortem.df.filter %>% filter(`DE status` == "LIV DEG") %>% arrange(logFC) %>% pull(symbol)
liharska.2025.pm.degs <- liharska.2025.postmortem.df.filter %>% filter(`DE status` == "PM DEG") %>% arrange(desc(logFC)) %>% pull(symbol)

all.gene.signatures[["liv.degs"]] <- liharska.2025.liv.degs[1:100]
all.gene.signatures[["pm.degs"]] <- liharska.2025.pm.degs[1:100]
all.gene.signatures[["liv.degs.all"]] <- liharska.2025.liv.degs
all.gene.signatures[["pm.degs.all"]] <- liharska.2025.pm.degs





##### Helper functions #####

# doublet finder pipeline
find.doublets <- function(seurat.obj, pcs.use.doublet.finder, resolution.use.doublet.finder = 0.8, qc.output.dir = paste0(results.dir, output.dir, "qc/")) {
  set.seed(123)
  
  # pre-processing
  seurat.obj <- FindVariableFeatures(seurat.obj)
  seurat.obj <- ScaleData(seurat.obj)
  
  seurat.obj <- RunPCA(seurat.obj)
  p <- ElbowPlot(seurat.obj, ndims = 50)
  pdf(paste0(qc.output.dir, "elbowplot_before_doublet_removal.pdf"), width = 8, height = 5)
  print(p)
  dev.off()
  seurat.obj <- FindNeighbors(seurat.obj, dims = pcs.use.doublet.finder)
  seurat.obj <- FindClusters(seurat.obj, resolution = resolution.use.doublet.finder)
  seurat.obj <- RunUMAP(seurat.obj, dims = pcs.use.doublet.finder)
  
  
  # find doublets
  sweep.res.list <- paramSweep(seurat.obj, PCs = pcs.use.doublet.finder, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)
  mpK <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  annotations <- seurat.obj@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)  ## ex: annotations <- sample14@meta.data$ClusteringResults
  
  nExp_poi <- round(0.075 * nrow(seurat.obj@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  seurat.obj <- doubletFinder(seurat.obj, PCs = pcs.use.doublet.finder, pN = 0.25, pK = mpK, nExp = nExp_poi.adj, sct = FALSE)
  seurat.obj$doubFind_res = seurat.obj@meta.data %>% dplyr::select(contains('DF.classifications'))
  seurat.obj$doubFind_score = seurat.obj@meta.data %>% dplyr::select(contains('pANN'))
  
  
  # umap with results
  p.clusters <- DimPlot(seurat.obj, group.by = "seurat_clusters", shuffle = T, label = T)
  p.doubfind.res <- DimPlot(seurat.obj, group.by = "doubFind_res", shuffle = T)
  p.doubFind.score <- FeaturePlot(seurat.obj, features = "doubFind_score", order = T)
  #p.celltype <- DimPlot(seurat.obj, group.by = "predicted.annotation_l2", shuffle = T) + cell.types.l2.color.scale
  p <- p.clusters + p.doubfind.res + p.doubFind.score
  
  pdf(paste0(qc.output.dir, "doublet_finder_umap.pdf"), width = 18, height = 6)
  print(p)
  dev.off()
  
  return(seurat.obj)
}



#' function for using AUCell to score gene sets in a Seurat object
#' 
#' a wrapper for AUCell
#' 
#' @param seurat.obj the seurat object you want to add the AUCell scores to
#' @param gene.signatures a list of gene signatures to pass into AUCell
add.aucell.score <- function(seurat.obj, gene.signatures) {
  expr_mat <- GetAssayData(seurat.obj, slot = "data")
  cells_rankings <- AUCell_buildRankings(expr_mat, plotStats = FALSE)
  
  cells_AUC <- AUCell_calcAUC(gene.signatures, cells_rankings)
  
  # Extract scores into matrix (cells × signatures)
  auc_scores <- as.data.frame(getAUC(cells_AUC)) %>% t()
  seurat.obj <- AddMetaData(seurat.obj, auc_scores)
  return(seurat.obj)
}



#' function for adding module score and then correcting the metadata column name (removes extra digits)
#' 
#' a wrapper for Seurat's AddModuleScore
#' 
#' @param seurat.obj the seurat object you want to add the module scores to
#' @param gene.signatures a list of gene signatures to pass into AddModuleScore
add.module.score <- function(seurat.obj, gene.signatures, assay.use = NULL, slot.use = "data", verbose = T) {
  
  if (is.null(assay.use)) {
    assay.use <- DefaultAssay(seurat.obj)
  }
  
  # overwrite existing metadata
  for (g in names(gene.signatures)) {
    if (g %in% colnames(seurat.obj@meta.data)) {
      seurat.obj@meta.data[g] <- NULL
    }
  }
  
  if (verbose) {
    seurat.obj <- AddModuleScore(
      object = seurat.obj,
      features = gene.signatures,
      name = names(gene.signatures),
      search = T,
      assay = assay.use,
      slot = slot.use
    )
  } else {
    suppressWarnings(
      seurat.obj <- AddModuleScore(
        object = seurat.obj,
        features = gene.signatures,
        name = names(gene.signatures),
        search = T,
        assay = assay.use,
        slot = slot.use
      )
    )
  }
  

  i <- 1
  for (g in names(gene.signatures)) {
    colnames(seurat.obj@meta.data) <- gsub(paste0(g, i), g, colnames(seurat.obj@meta.data))
    
    if (i == 11) {
      colnames(seurat.obj@meta.data) <- gsub(paste0(g, 1), g, colnames(seurat.obj@meta.data))
    }
    i <- i + 1
  }
  
  return(seurat.obj)
}







# method to assign a neftel type to each cell on the neftel plot
# used in the score.neftel function
get.neftel.type <- function(row) {
  y <- as.numeric(row[['stemness']])
  if (y > 0) {
    diff <- as.numeric(row[['opc_sig']]) - as.numeric(row[['npc_sig']])
    if (diff > 0) {
      return('OPC-like')
    } else {
      return('NPC-like')
    }
  } else {
    diff <- as.numeric(row[['ac_sig']]) - as.numeric(row[['mes_sig']])
    if (diff > 0) {
      return('AC-like')
    } else {
      return('MES-like')
    }
  }
}

# method to compute x axis placements of each cell on the neftel plot
# used in the score.neftel function
get.neftel.x <- function(row) {
  y <- as.numeric(row[['stemness']])
  if (y > 0) {
    diff <- as.numeric(row[['opc_sig']]) - as.numeric(row[['npc_sig']])
  } else {
    diff <- as.numeric(row[['ac_sig']]) - as.numeric(row[['mes_sig']])
  }
  
  if (diff > 0) {
    x <- log2(abs(diff) + 1)
  } else {
    x <- -log2(abs(diff) + 1)
  }
  return(x)
}


#' function for classifying neftel types for each cell and assigning an x and y coordinate for each cell on the neftel plot
#'
#' @param seurat.obj the seurat object to assign neftel types for
#' @param assay.use assay to use for AddModuleScore for neftel scoring. Default = "RNA"
#' 
score.neftel <- function(seurat.obj, signatures.dir="/Users/pengja/PennNeurosurgery Dropbox/Jackie Peng/SilverbushLab/Shared_resources/Signatures/Neftel_GBM_cell_states/SingleCell_GBM/", assay.use = NULL) {
  
  if (is.null(assay.use)) {
    assay.use <- DefaultAssay(seurat.obj)
  }
  
  # import signatures
  ac_sig <- readLines(paste0(signatures.dir, 'ac_genes.txt'))
  mes1_sig <- readLines(paste0(signatures.dir, 'mes1_genes.txt'))
  mes2_sig <- readLines(paste0(signatures.dir, 'mes2_genes.txt'))
  npc1_sig <- readLines(paste0(signatures.dir, 'npc1_genes.txt'))
  npc2_sig <- readLines(paste0(signatures.dir, 'npc2_genes.txt'))
  opc_sig <- readLines(paste0(signatures.dir, 'opc_genes.txt'))
  
  signatures <- list(ac_sig, mes1_sig, mes2_sig, npc1_sig, npc2_sig, opc_sig)
  names <- c("ac_sig.2", "mes1_sig.2", "mes2_sig.2", "npc1_sig.2", "npc2_sig.2", "opc_sig.2")
  
  seurat.obj <- AddModuleScore(seurat.obj, features = signatures, name = names, search=TRUE, assay = assay.use, slot='data')
  seurat.obj$ac_sig <- seurat.obj$ac_sig.21
  seurat.obj$opc_sig <- seurat.obj$opc_sig.26
  
  # compute stemness
  seurat.obj$mes_sig <- apply(seurat.obj@meta.data, 1, function(row) as.numeric(max(row[['mes1_sig.22']], row[['mes2_sig.23']])))
  seurat.obj$npc_sig <- apply(seurat.obj@meta.data, 1, function(row) as.numeric(max(row[['npc1_sig.24']], row[['npc2_sig.25']])))
  seurat.obj$stemness <- apply(seurat.obj@meta.data, 1, function(row) as.numeric(max(as.numeric(row[['opc_sig']]), as.numeric(row[['npc_sig']])) - max(as.numeric(row[['ac_sig']]), as.numeric(row[['mes_sig']]))))
  
  seurat.obj$neftel_x <- apply(seurat.obj@meta.data, 1, get.neftel.x)
  seurat.obj$neftel_y <- seurat.obj$stemness
  seurat.obj$neftel_type <- apply(seurat.obj@meta.data, 1, get.neftel.type)
  return(seurat.obj)
}


##### Code for calculating IDH mutant glioma malignant state based on Tirosh et al. signatures #####

# method from PC1-associated genes and lineage scores, Tirosh et al. 2016
calculate.lineage.score <- function(row, tumor.type="IDH-A") {
  
  if (tumor.type == "IDH-O") {
    ac.score <- as.numeric(row[["AC-like IDH-O4"]])
    oc.score <- as.numeric(row[["OC-like IDH-O5"]])
  } else if (tumor.type == "IDH-A") {
    ac.score <- as.numeric(row[["AC-like IDH-A2"]])
    oc.score <- as.numeric(row[["OC-like IDH-A3"]])
  } else {
    stop(paste("Unknown tumor.type:", tumor.type))
  }
  
  return(c(ac.score, oc.score))
}

# same as method above, but if AC score > OC score, converts the score to negative for plotting
calculate.lineage.score.plot <- function(row, tumor.type="IDH-A", ac.score.col="ac.score", oc.score.col="oc.score") {
  
  ac.score <- as.numeric(row[[ac.score.col]])
  oc.score <- as.numeric(row[[oc.score.col]])
  
  if (is.na(ac.score) || is.na(oc.score)) { # when calculating on cells with no existing metadata
    return(NA)
  }
  
  # a little bit artificial: according to methods, if lineage score < 0, 
  # assigne a random number between 0 and 0.15
  if (max(ac.score, oc.score) < 0) {
    lineage.score <- runif(1, min = 0, max = 0.15)
    return(lineage.score)
  }
  
  if (ac.score > oc.score) {
    lineage.score <- -ac.score
  } else {
    lineage.score <- oc.score
  }
  return(as.numeric(lineage.score))
}

# method from PC2 and PC3 associated genes and stemness scores, Tirosh et al. 2016
calculate.stemness.score <- function(row) {
  ac.score <- as.numeric(row[['ac.score']])
  oc.score <- as.numeric(row[['oc.score']])
  lineage.score <- max(ac.score, oc.score)
  stemness.score <- as.numeric(row[['stem-like1']]) - lineage.score
  return(stemness.score)
}

calculate.cell.assignment <- function(row, threshold=0.25) {
  stemness.score <- as.numeric(row[['stemness.score']])
  ac.score <- as.numeric(row[['ac.score']])
  oc.score <- as.numeric(row[['oc.score']])
  lineage.score <- max(ac.score, oc.score)
  
  if (max(stemness.score, lineage.score) > threshold) {
    if (stemness.score - lineage.score > threshold) {
      cell.assignment <- "stem-like"
    } else if ((ac.score - stemness.score > threshold) && (ac.score - oc.score > threshold)) {
      cell.assignment <- "AC-like"
    } else if ((oc.score - stemness.score > threshold) && (oc.score - ac.score > threshold)) {
      cell.assignment <- "OC-like"
    } else {
      cell.assignment <- "undifferentiated"
    }
  } else {
    cell.assignment <- "undifferentiated"
  }
  
  return(cell.assignment)
}

score.idh.mutant.state <- function(seurat.obj, signatures.dir = "/Users/pengja/PennNeurosurgery Dropbox/Jackie Peng/SilverbushLab/Shared_resources/Signatures/Tirosh_IDHmut_cell_states/") {
  idh.a.stemness.genes <- read_lines(paste0(signatures.dir, "ASTROCYTOMA_stemness_genes.txt"))
  idh.o.stemness.genes <- read_lines(paste0(signatures.dir, "OLIGODENDROGLIOMA_stemness.csv"))
  stopifnot(identical(idh.a.stemness.genes, idh.o.stemness.genes)) # stemness genes should be the same
  
  idh.a.stemness.genes <- replace(idh.a.stemness.genes, idh.a.stemness.genes=="HN1", "JPT1")  
  idh.a.stemness.genes <- replace(idh.a.stemness.genes, idh.a.stemness.genes=="ATP5J", "ATP5PF")  
  
  idh.a.ac.genes <- read_lines(paste0(signatures.dir, "ASTROCYTOMA_ac_genes.txt"))
  idh.a.ac.genes <- replace(idh.a.ac.genes, idh.a.ac.genes=="PPAP2B", "PLPP3")
  
  idh.a.oc.genes <- read_lines(paste0(signatures.dir, "ASTROCYTOMA_opc_genes.txt")) # no alternative symbol for SOX-OT
  idh.a.oc.genes <- replace(idh.a.oc.genes, idh.a.oc.genes=="LPPR1", "PLPPR1")
  
  idh.o.ac.genes <- read_lines(paste0(signatures.dir, "OLIGODENDROGLIOMA_astro.csv"))
  idh.o.oc.genes <- read_lines(paste0(signatures.dir, "OLIGODENDROGLIOMA_oligo.csv"))
  
  idh.glioma.gene.signatures <- list()
  idh.glioma.gene.signatures[['stem-like']] <- idh.a.stemness.genes
  idh.glioma.gene.signatures[['AC-like IDH-A']] <- idh.a.ac.genes
  idh.glioma.gene.signatures[['OC-like IDH-A']] <- idh.a.oc.genes
  idh.glioma.gene.signatures[['AC-like IDH-O']] <- idh.o.ac.genes
  idh.glioma.gene.signatures[['OC-like IDH-O']] <- idh.o.oc.genes
  
  seurat.obj <- AddModuleScore(seurat.obj, features=idh.glioma.gene.signatures, name=names(idh.glioma.gene.signatures), slot="data", nbin=25, ctrl = 100, seed = 1)
  
  output <- apply(seurat.obj@meta.data, 1, calculate.lineage.score)
  seurat.obj$ac.score <- as.numeric(output[1,])
  seurat.obj$oc.score <- as.numeric(output[2,]) 
  seurat.obj$stemness.score <- apply(seurat.obj@meta.data, 1, calculate.stemness.score)
  seurat.obj$lineage.score <- apply(seurat.obj@meta.data, 1, function(row) max(as.numeric(row[['ac.score']]), as.numeric(row[['oc.score']])))
  seurat.obj$lineage.score.plot <- apply(seurat.obj@meta.data, 1, calculate.lineage.score.plot)
  seurat.obj$idh.mutant.state <- apply(seurat.obj@meta.data, 1, calculate.cell.assignment)
  
  return(seurat.obj)
}






# clustering entropy (unweighted and weighted)
calculate.clustering.entropy <- function(clusters, annotations) {
  df <- data.frame(clusters, annotations)
  
  # entropy per cluster
  entropy.per.cluster <- df %>%
    dplyr::group_by(clusters) %>%
    summarise(
      cluster.size = n(),
      entropy = {
        counts <- table(annotations)
        entropy.empirical(counts, unit = "log2")
      },
      .groups = "drop"
    )
  
  # mean entropy unweighted and weighted
  unweighted.entropy <- mean(entropy.per.cluster$entropy)
  weighted.entropy <- weighted.mean(entropy.per.cluster$entropy, w = entropy.per.cluster$cluster.size)
  
  return (list(
    unweighted.entropy = unweighted.entropy,
    weighted.entropy = weighted.entropy,
    entropy.per.cluster = entropy.per.cluster
  ))
}





##### Plotting helper functions #####

format.p.value <- function(p) {
  p.label <- paste0(
    "p ",
    ifelse(
      p < 1e-6,
      "< 1e-06",
      ifelse(
        p < 0.001,
        paste0("= ", formatC(p, format = "e", digits = 2)),
        paste0("= ", formatC(p, format = "f", digits = 3))
      )
    )
  )
  return(p.label)
}

# Function to add the number of cells per cell type/group in the legend
# Returns a color scale with counts that can be added to a DimPlot
get.count.labels <- function(metadata.col, mapping) {
  cell.counts <- table(metadata.col)
  cell.counts <- cell.counts[intersect(names(mapping), names(cell.counts))]
  
  labels.w.counts <- paste0(
    names(cell.counts),
    " (n=", as.integer(cell.counts), ")"
  )
  names(labels.w.counts) <- names(cell.counts)
  color.scale.w.counts <- scale_color_manual(
    labels = labels.w.counts,
    breaks = names(labels.w.counts),
    values = mapping[names(cell.counts)]
  )
  return(color.scale.w.counts)
}



plot_composition <- function(
    metadata,
    x.var,                     
    fill.var,                  
    fill.levels = NULL,
    x.levels = NULL,
    fill.scale = NULL,
    label.threshold = 5,
    label.digits = 1,
    xlab = NULL,
    ylab = "Percent of sample (%)",
    legend.title = NULL,
    rotate.x = TRUE,
    ylim.max = 105,
    add.text = TRUE
) {
  
  df <- metadata %>%
    dplyr::mutate(!!sym(fill.var) := as.character(.data[[fill.var]])) %>%
    count(.data[[x.var]], .data[[fill.var]]) %>%
    group_by(.data[[x.var]]) %>%
    mutate(
      total_n = sum(n),
      percent = 100 * n / total_n,
      x_with_n = sprintf("%s\n(n=%d)", .data[[x.var]], total_n)
    ) %>%
    filter(total_n > 0) %>%
    ungroup()
  
  if (!is.null(fill.levels)) {
    df[[fill.var]] <- factor(df[[fill.var]], levels = fill.levels)
  }
  
  if (!is.null(x.levels)) {
    present_x <- x.levels[x.levels %in% unique(df[[x.var]])]
    
    order_map <- df %>%
      dplyr::select(all_of(c(x.var, "x_with_n"))) %>%
      distinct() %>%
      dplyr::slice(match(present_x, .data[[x.var]]))
    
    df$x_with_n <- factor(df$x_with_n, levels = order_map$x_with_n)
  }
  
  p <- ggplot(df, aes(x = x_with_n, y = percent, fill = .data[[fill.var]])) +
    geom_col(position = "stack", color = "white", linewidth = 0.2) +
    labs(
      x = xlab %||% x.var,
      y = ylab,
      fill = legend.title %||% fill.var
    ) +
    coord_cartesian(ylim = c(0, ylim.max))
  
  if (add.text) {
    p <- p +
      geom_text(
        aes(label = ifelse(percent < label.threshold, "", sprintf("%.*f%%", label.digits, percent))),
        position = position_stack(vjust = 0.5),
        size = 5
      )
  }
  
  if (!is.null(fill.scale)) {
    p <- p + fill.scale
  }
  
  if (rotate.x) {
    p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  return(p)
}





do.boxplot <- function(
    data,
    x.variable,
    y.variable,
    boxplot.fill.variable = NULL,
    fill.scale = NULL,
    xlab.label = NULL,
    ylab.label = NULL,
    ggtitle.label = NULL,
    show.points = FALSE,
    points.alpha = 0.6,
    outlier.shape = TRUE,
    comparisons = NULL,
    alternative = "two.sided",
    y.adjust = 0,             
    step.increase = 0.1,      
    p.adj.method = "fdr",
    save.figure = TRUE,
    output.filepath = NULL,
    image.width = 10,
    image.height = 6,
    return.plot = FALSE,
    verbose = FALSE,
    facet.by = NULL,
    vars.to.regress = NULL,
    free_y = TRUE,
    facet.nrow = 1
) {
  
  ## -----------------------------
  ## Data preparation
  ## -----------------------------
  data$x <- data[[x.variable]]
  data$y <- data[[y.variable]]
  
  if (!is.null(vars.to.regress)) {
    model <- lm(reformulate(termlabels = vars.to.regress, response = "y"), data = data)
    data$y <- resid(model)
  }
  
  if (!is.null(boxplot.fill.variable)) {
    data$v.fill <- data[[boxplot.fill.variable]]
  } else {
    data$v.fill <- data$x
  }
  
  .filter_comparisons <- function(comparisons, x_values) {
    if (is.null(comparisons)) return(NULL)
    x_values <- unique(as.character(x_values))
    purrr::keep(comparisons, ~ all(.x %in% x_values))
  }
  
  ## -----------------------------
  ## Helper: single-group boxplot
  ## -----------------------------
  .make_single_boxplot <- function(d, local_ylim_min = NULL, local_ylim_max = NULL, current_y_adj = 0) {
    
    d$x <- factor(d$x)
    x_counts <- d %>% dplyr::filter(!is.na(y)) %>% dplyr::count(x)
    x_labels <- setNames(paste0(x_counts$x, "\n(n = ", x_counts$n, ")"), x_counts$x)
    
    valid_comparisons <- .filter_comparisons(comparisons, d$x)
    stat.test <- NULL
    
    if (!is.null(valid_comparisons) && length(valid_comparisons) > 0 && n_distinct(d$x) >= 2) {
      stat.test <- d %>%
        dplyr::filter(!is.na(y)) %>%
        rstatix::wilcox_test(y ~ x, comparisons = valid_comparisons, alternative = alternative) %>%
        rstatix::adjust_pvalue(method = p.adj.method) %>%
        rstatix::add_significance("p") %>%
        rstatix::add_xy_position(x = "x", step.increase = step.increase)
      
      stat.test$p.label <- sapply(stat.test$p.adj, format.p.value)
      stat.test$y.position <- stat.test$y.position + current_y_adj
      
      if (is.null(local_ylim_max)) {
        local_ylim_max <- max(stat.test$y.position) * 1.1
      }
    } else {
      if (is.null(local_ylim_max)) local_ylim_max <- max(d$y, na.rm = TRUE) * 1.1
    }
    
    if (is.null(local_ylim_min)) local_ylim_min <- min(d$y, na.rm = TRUE)
    
    p <- ggplot(d, aes(x, y)) +
      geom_boxplot(aes(fill = v.fill), outlier.shape = outlier.shape) +
      scale_x_discrete(labels = x_labels) +
      coord_cartesian(ylim = c(local_ylim_min, local_ylim_max)) +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    
    if (!is.null(boxplot.fill.variable)) p <- p + fill.scale
    if (show.points) p <- p + geom_point(position = position_jitter(width = 0.15), alpha = points.alpha)
    if (!is.null(stat.test) && nrow(stat.test) > 0) {
      p <- p + ggpubr::stat_pvalue_manual(stat.test, label = "p.label", tip.length = 0.01, size = 6)
    }
    
    return(p)
  }
  
  ## -----------------------------
  ## Plot logic
  ## -----------------------------
  if (!is.null(facet.by)) {
    
    # Compute shared y limits across all facets when free_y = FALSE
    global_ylim_min <- if (!free_y) min(data$y, na.rm = TRUE) else NULL
    
    global_ylim_max <- if (!free_y) {
      all_y_adj <- if (!is.null(names(y.adjust))) {
        sapply(unique(as.character(data[[facet.by]])), function(fv) {
          if (fv %in% names(y.adjust)) y.adjust[[fv]] else y.adjust[1]
        })
      } else {
        rep(y.adjust[1], length(unique(data[[facet.by]])))
      }
      
      if (!is.null(comparisons)) {
        facet_vals <- unique(as.character(data[[facet.by]]))
        bracket_maxes <- purrr::map2_dbl(facet_vals, all_y_adj, function(fv, adj) {
          d <- data[as.character(data[[facet.by]]) == fv, ]
          d$x <- factor(d$x)
          valid_comps <- .filter_comparisons(comparisons, d$x)
          if (!is.null(valid_comps) && length(valid_comps) > 0 && dplyr::n_distinct(d$x) >= 2) {
            st <- tryCatch(
              d %>%
                dplyr::filter(!is.na(y)) %>%
                rstatix::wilcox_test(y ~ x, comparisons = valid_comps, alternative = alternative) %>%
                rstatix::adjust_pvalue(method = p.adj.method) %>%
                rstatix::add_xy_position(x = "x", step.increase = step.increase),
              error = function(e) NULL
            )
            if (!is.null(st)) max(st$y.position + adj) * 1.1 else max(d$y, na.rm = TRUE) * 1.1
          } else {
            max(d$y, na.rm = TRUE) * 1.1
          }
        })
        max(bracket_maxes)
      } else {
        max(data$y, na.rm = TRUE) * 1.1
      }
    } else {
      NULL
    }
    
    plots <- split(data, data[[facet.by]]) %>%
      purrr::imap(function(d, facet_value) {
        
        this_adj <- if (!is.null(names(y.adjust)) && as.character(facet_value) %in% names(y.adjust)) {
          y.adjust[[as.character(facet_value)]]
        } else {
          y.adjust[1] 
        }
        
        facet_title <- if (!is.null(ggtitle.label)) {
          if (!is.null(names(ggtitle.label)) &&
              as.character(facet_value) %in% names(ggtitle.label)) {
            ggtitle.label[[as.character(facet_value)]]
          } else {
            ggtitle.label[1]
          }
        } else {
          paste(facet_value)
        }
        
        p <- .make_single_boxplot(d, local_ylim_min = global_ylim_min, local_ylim_max = global_ylim_max, current_y_adj = this_adj) +
          labs(x = xlab.label, y = ylab.label) +
          ggtitle(facet_title)
        p
      })
    
    p <- cowplot::plot_grid(plotlist = plots, nrow = facet.nrow, align = "h")
    
  } else {
    p <- .make_single_boxplot(data, current_y_adj = y.adjust[1])
    if (!is.null(ggtitle.label)) p <- p + ggtitle(ggtitle.label)
    p <- p + labs(x = xlab.label, y = ylab.label)
  }
  
  ## -----------------------------
  ## Save and Return
  ## -----------------------------
  if (save.figure) {
    pdf(output.filepath, width = image.width, height = image.height)
    print(p)
    dev.off()
  }
  
  if (return.plot) return(p)
}





geneSetOverlapRow <- function(setA_name, setB_name, setA, setB, universe) {
  res <- geneSetOverlapTest(setA, setB, universe)
  
  data.frame(
    setA = setA_name,
    setB = setB_name,
    overlap = res$overlap,
    size.A = res$size.A,
    size.B = res$size.B,
    universe.size = res$universe.size,
    odds.ratio = res$odds.ratio,
    fold.enrichment = res$fold.enrichment,
    p.value = res$p.value,
    stringsAsFactors = FALSE
  )
}


geneSetOverlapTest <- function(setA,
                               setB,
                               universe,
                               alternative = "two.sided") {
  # Restrict to universe
  setA <- intersect(unique(setA), universe)
  setB <- intersect(unique(setB), universe)
  
  # Sizes
  A <- length(setA)
  B <- length(setB)
  k <- length(intersect(setA, setB))
  N <- length(universe)
  
  # Contingency table
  mat <- matrix(
    c(k,
      A - k,
      B - k,
      N - A - B + k),
    nrow = 2
  )
  
  # Fisher's exact test
  ft <- fisher.test(mat, alternative = alternative)
  
  # Fold enrichment
  fold_enrichment <- (k / A) / (B / N)
  
  # Return a tidy list
  list(
    p.value = ft$p.value,
    odds.ratio = unname(ft$estimate),
    fold.enrichment = fold_enrichment,
    overlap = k,
    size.A = A,
    size.B = B,
    universe.size = N,
    contingency.table = mat
  )
}


plot_gene_expression <- function(seurat_obj_metadata, 
                                 features, 
                                 cell_type_col = "CellType",
                                 cell_type_levels = NULL,
                                 color_values = NULL,
                                 color_breaks = NULL,
                                 show_color_legend = FALSE,
                                 by_sample = TRUE, 
                                 num_rows = 1,
                                 strip_text_size = 14,
                                 plot_width = 16,
                                 plot_height = 8,
                                 mean_size = NULL,
                                 output_path = NULL,
                                 facet_wrap_variable = NULL,
                                 min_cells = 10,
                                 alternative = "two.sided"
                                 ) {
  
  # 1. Extract metadata and subset/pivot for chosen features
  df <- seurat_obj_metadata %>%
    pivot_longer(cols = all_of(features), names_to = "Gene", values_to = "Expression")
  
  # Rename cell type column to a standard name for internal use
  df <- df %>% dplyr::rename(CellType = !!sym(cell_type_col))
  
  # Handle initial factor leveling dynamically
  if (!is.null(cell_type_levels)) {
    df$CellType <- factor(df$CellType, levels = cell_type_levels)
  } else {
    df$CellType <- factor(df$CellType)
  }
  df$Gene <- factor(df$Gene, levels = features)
  
  # 2. Statistical Testing and Data Aggregation
  if (by_sample) {
    if (!all(c("Patient", "Sample", "Timepoint") %in% colnames(df))) {
      stop("For 'by_sample = TRUE', metadata must contain 'Patient', 'Sample', and 'Timepoint' columns.")
    }
    
    df <- df %>% 
      group_by(Patient, Sample, Timepoint, CellType, Gene) %>%
      summarize(Expression = mean(Expression),
                count = n(),
                .groups = "drop") %>%
      group_by(Sample, Gene) %>%
      filter(n() == 2 & all(count > min_cells)) %>% 
      ungroup() %>%
      arrange(Sample, CellType)
    
    # Counts unique samples per cell type
    counts_df <- df %>%
      group_by(CellType) %>%
      summarize(n_count = n_distinct(Sample), .groups = "drop")
    
  } else {
    # Counts unique cells (rows) per cell type
    counts_df <- df %>%
      group_by(CellType) %>%
      summarize(n_count = n(), .groups = "drop")
  }
  
  # Map counts back to create new explicit factor levels
  label_mapping <- counts_df %>%
    mutate(new_label = paste0(CellType, "\n(n=", n_count, ")"))
  
  # Update df with the new labels and update factor levels
  df <- df %>%
    left_join(label_mapping, by = "CellType") %>%
    mutate(CellType = factor(new_label, levels = label_mapping$new_label))
  
  # 3. Wilcoxon Tests 
  if (by_sample) {
    stat.test <- df %>%
      group_by(Gene) %>%
      wilcox_test(Expression ~ CellType, paired = TRUE) %>%
      adjust_pvalue(method = "fdr") %>%
      add_significance("p.adj") %>%
      add_xy_position(x = "CellType", step.increase = 0.05) 
  } else {
    stat.test <- df %>%
      group_by(Gene) %>%
      wilcox_test(Expression ~ CellType, paired = FALSE, alternative = alternative) %>%
      adjust_pvalue(method = "fdr") %>%
      add_significance("p.adj") %>%
      add_xy_position(x = "CellType", step.increase = 0.05) 
  }
  
  if (nrow(stat.test) > 0) {
    stat.test$p.label <- sapply(stat.test$p.adj, format.p.value)
  }
  
  # Calculate summary for fallback geom_point
  summary_df <- df %>%
    group_by(Gene, CellType) %>%
    summarise(
      mean_expr = mean(Expression),
      pct_expressing = mean(Expression > 0) * 100,
      .groups = "drop"
    )
  
  # 4. Build the Plot
  p <- ggplot(df, aes(x = CellType, y = Expression, color = CellType)) +
    geom_violin(trim = FALSE, scale = "width", fill = NA) + 
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA)
  
  if (by_sample) {
    p <- p + geom_point(position = position_jitter(width = 0.1), alpha = 0.6)
  } else {
    if (is.null(mean_size)) {
      p <- p + geom_point(
        data = summary_df,
        aes(x = CellType, y = mean_expr, size = pct_expressing),
        shape = 23, fill = "red", color = "black", inherit.aes = FALSE
      ) + labs(size = "% expressing")
    } else if (mean_size > 0) {
      p <- p + geom_point(
        data = summary_df,
        aes(x = CellType, y = mean_expr, size = mean_size),
        shape = 23, fill = "red", color = "black", inherit.aes = FALSE, show.legend = FALSE
      ) + labs(size = NULL)
    }
  }
  
  # Handle Faceting Dynamically
  if (!is.null(facet_wrap_variable)) {
    facet_formula <- as.formula(paste("~", facet_wrap_variable, "+ Gene"))
    p <- p + facet_wrap(facet_formula, scales = "free_y")
  } else {
    p <- p + facet_wrap(~Gene, nrow = num_rows, scales = "free_y")
  }
  
  # Theme and labeling
  p <- p + 
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "right"
    ) +
    labs(x = "Cell type", y = "Expression") +
    theme(strip.text = element_text(size = strip_text_size))
  
  if (nrow(stat.test) > 0) {
    p <- p + stat_pvalue_manual(stat.test, label = "p.label", size = 5)
  }
  
  # 5. Apply color scale using direct arguments
  if (!is.null(color_values) && !is.null(color_breaks)) {
    new_breaks <- label_mapping$new_label[match(color_breaks, label_mapping$CellType)]
    named_values <- setNames(color_values, new_breaks)
    
    p <- p + scale_color_manual(
      breaks = new_breaks,
      values = named_values,
      guide  = if (show_color_legend) "legend" else "none"
    )
  }
  
  # 6. Save and Export
  if (!is.null(output_path)) {
    pdf(output_path, width = plot_width, height = plot_height)
    print(p)
    dev.off()
    message("Plot successfully saved to: ", output_path)
  }
  
  return(p)
}





# For Visium HD analysis
get_cells_within_radius <- function(centroids, center_id, radius, remove_center = TRUE) {
  center <- centroids[centroids$cell_id == center_id, c("x", "y")]
  
  d <- sqrt(
    (centroids$x - center$x)^2 +
      (centroids$y - center$y)^2
  )
  
  neighborhood <- centroids[d <= radius, ]$cell_id
  
  if (remove_center) {
    neighborhood <- setdiff(neighborhood, center_id)
  }
  
  return(neighborhood)
}




make_neighborhood_plot <- function(df, comparisons = NULL, y_var = "reactive.oligo.general", 
                                   y_var_label = NULL, alternative.use = "less", 
                                   facet_var = NULL, filename, width = 10, height = 8,
                                   ref_group = NULL) {
  
  df <- df %>% filter(!is.na(neighborhood.annotation), !is.na(.data[[y_var]]))
  df$neighborhood.annotation <- as.factor(df$neighborhood.annotation)
  
  valid_groups <- df %>%
    group_by(neighborhood.annotation) %>%
    summarise(n = n(), .groups = "drop") %>%
    filter(n >= 3) %>%
    pull(neighborhood.annotation)
  
  df <- df %>% filter(neighborhood.annotation %in% valid_groups)
  df$neighborhood.annotation <- droplevels(df$neighborhood.annotation)
  
  # Build comparisons from ref_group if provided, overriding any passed comparisons
  if (!is.null(ref_group)) {
    if (!ref_group %in% levels(df$neighborhood.annotation)) {
      stop("ref_group '", ref_group, "' not found in neighborhood.annotation after filtering.")
    }
    other_groups <- setdiff(levels(df$neighborhood.annotation), ref_group)
    comparisons <- lapply(other_groups, function(g) c(ref_group, g))
  }
  
  if (is.null(comparisons)) {
    stop("Either 'comparisons' or 'ref_group' must be provided.")
  }
  
  stat.test <- NULL
  formula <- as.formula(paste(y_var, "~ neighborhood.annotation"))
  
  # Helper: filter comparisons to those with sufficient data in a given data frame
  get_valid_comparisons <- function(sub_df, comps) {
    Filter(function(comp) {
      x_vals <- sub_df[[y_var]][sub_df$neighborhood.annotation == comp[1]]
      y_vals <- sub_df[[y_var]][sub_df$neighborhood.annotation == comp[2]]
      x_vals <- x_vals[!is.na(x_vals)]
      y_vals <- y_vals[!is.na(y_vals)]
      length(x_vals) >= 2 && length(y_vals) >= 2
    }, comps)
  }
  
  if (is.null(facet_var)) {
    
    valid_comparisons <- get_valid_comparisons(df, comparisons)
    
    if (length(valid_comparisons) == 0) {
      message("Skipping stat test: insufficient data for any requested comparisons.")
    } else {
      if (length(valid_comparisons) < length(comparisons)) {
        dropped <- setdiff(
          sapply(comparisons, paste, collapse = " vs "),
          sapply(valid_comparisons, paste, collapse = " vs ")
        )
        message("Dropping comparisons due to insufficient data: ", paste(dropped, collapse = "; "))
      }
      
      stat.test <- df %>%
        wilcox_test(formula, 
                    comparisons = valid_comparisons, 
                    alternative = alternative.use) %>%
        add_significance() %>%
        adjust_pvalue(method = "fdr") %>%
        add_xy_position(x = "neighborhood.annotation")
      
      stat.test$p.label <- sapply(stat.test$p.adj, format.p.value)
    }
    
  } else {
    
    run_facet_wilcox <- function(sub_df) {
      valid_comparisons <- get_valid_comparisons(sub_df, comparisons)
      
      if (length(valid_comparisons) == 0) return(NULL)
      
      if (length(valid_comparisons) < length(comparisons)) {
        dropped <- setdiff(
          sapply(comparisons, paste, collapse = " vs "),
          sapply(valid_comparisons, paste, collapse = " vs ")
        )
        facet_val <- unique(sub_df[[facet_var]])
        message("Facet '", facet_val, "': dropping comparisons due to insufficient data: ", paste(dropped, collapse = "; "))
      }
      
      tryCatch(
        sub_df %>%
          wilcox_test(formula, 
                      comparisons = valid_comparisons, 
                      alternative = alternative.use),
        error = function(e) {
          message("Skipping comparison in facet due to error: ", conditionMessage(e))
          NULL
        }
      )
    }
    
    stat.test <- df %>%
      group_by(.data[[facet_var]]) %>%
      nest() %>%
      mutate(res = map(data, run_facet_wilcox)) %>%
      select(-data) %>%
      filter(!map_lgl(res, is.null)) %>%
      filter(map_lgl(res, is.data.frame))
    
    if (nrow(stat.test) == 0) {
      message("Skipping stat test: no facets have sufficient data for comparisons")
      stat.test <- NULL
    } else {
      stat.test <- tidyr::unnest(stat.test, cols = all_of("res"))
      
      stat.test <- stat.test %>%
        add_significance() %>%
        adjust_pvalue(method = "fdr")
      
      y_maxes <- df %>%
        group_by(.data[[facet_var]], neighborhood.annotation) %>%
        summarise(y_max = max(.data[[y_var]], na.rm = TRUE), .groups = "drop") %>%
        group_by(.data[[facet_var]]) %>%
        summarise(y_max = max(y_max, na.rm = TRUE), .groups = "drop")
      
      stat.test <- stat.test %>%
        left_join(y_maxes, by = facet_var) %>%
        group_by(.data[[facet_var]]) %>%
        mutate(
          xmin = match(group1, levels(df$neighborhood.annotation)),
          xmax = match(group2, levels(df$neighborhood.annotation)),
          y.position = y_max + (row_number() * 0.1 * (max(df[[y_var]], na.rm = TRUE) - min(df[[y_var]], na.rm = TRUE)))
        ) %>%
        ungroup()
      
      stat.test$p.label <- sapply(stat.test$p.adj, format.p.value)
    }
  }
  
  # --- Plotting ---
  
  if (is.null(y_var_label)) {
    y_var_label <- y_var
  }
  
  p <- ggplot(df, aes(neighborhood.annotation, .data[[y_var]])) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(position = position_jitter(width = 0.2), alpha = 0.1) +
    xlab("Neighborhood") +
    ylab(y_var_label) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  if (!is.null(stat.test)) {
    p <- p + stat_pvalue_manual(stat.test, label = "p.label", label.size = 6, step.increase = 0.1)
  }
  
  if (!is.null(facet_var)) {
    p <- p + facet_wrap(as.formula(paste("~", facet_var)), scales = "free")
  }
  
  pdf(filename, width = width, height = height)
  print(p)
  dev.off()
}







##### DENOISING METRICS #####

mixing_metric <- function(
    X,                        # numeric matrix: cells x dims (e.g., PCA)
    groups,                   # vector/factor of length n_cells
    k = 30,
    metric = c("seurat", "entropy", "lisi"),
    return_per_cell = TRUE,
    scale_rows = FALSE
){
  metric <- match.arg(metric)
  X <- as.matrix(X)
  stopifnot(nrow(X) == length(groups))
  n <- nrow(X)
  groups <- as.factor(groups)
  G <- nlevels(groups)
  
  if (scale_rows) {
    X <- scale(X)
  }
  
  # kNN search (exclude self if returned)
  kn <- get.knn(X, k = k + 1)
  idx <- kn$nn.index
  dist <- kn$nn.dist
  
  # Drop self neighbors
  is_self_first <- dist[,1] == 0
  if (any(is_self_first)) {
    idx_trim <- matrix(NA_integer_, n, k)
    dist_trim <- matrix(NA_real_, n, k)
    for (i in seq_len(n)) {
      if (is_self_first[i]) {
        idx_trim[i,] <- idx[i, 2:(k+1)]
        dist_trim[i,] <- dist[i, 2:(k+1)]
      } else {
        idx_trim[i,] <- idx[i, 1:k]
        dist_trim[i,] <- dist[i, 1:k]
      }
    }
    idx <- idx_trim
    dist <- dist_trim
  } else {
    idx <- idx[, 1:k, drop = FALSE]
    dist <- dist[, 1:k, drop = FALSE]
  }
  
  neigh_groups <- matrix(groups[idx], nrow = n, ncol = k)
  
  if (metric == "seurat") {
    same <- (neigh_groups == matrix(groups, n, k))
    first_other_rank <- apply(!same, 1, function(v) {
      p <- which(v)[1]
      if (is.na(p)) k + 1 else p
    })
    per_cell <- 1 - (first_other_rank - 1) / k
  }
  
  if (metric == "entropy") {
    per_cell <- numeric(n)
    for (i in seq_len(n)) {
      tab <- table(neigh_groups[i,])
      p <- as.numeric(tab) / sum(tab)
      H <- if (all(p == 0)) 0 else -sum(p * log(p))
      per_cell[i] <- if (G > 1) H / log(G) else 0
    }
  }
  
  if (metric == "lisi") {
    per_cell <- numeric(n)
    for (i in seq_len(n)) {
      tab <- table(neigh_groups[i,])
      p <- as.numeric(tab) / sum(tab)
      lisi <- 1 / sum(p^2)
      per_cell[i] <- if (G > 1) (lisi - 1) / (G - 1) else 0
    }
  }
  
  out <- list(
    metric = metric,
    k = k,
    n_groups = G,
    per_cell = if (return_per_cell) per_cell else NULL,
    mean = mean(per_cell),
    median = median(per_cell),
    sd = sd(per_cell)
  )
  class(out) <- "mixing_metric"
  out
}

mixing_metric_seurat <- function(
    object,
    group.by,
    reduction = "pca",
    dims = 1:30,
    k = 30,
    metric = c("seurat", "entropy", "lisi"),
    return_per_cell = TRUE
){
  metric <- match.arg(metric)
  stopifnot(group.by %in% colnames(object@meta.data))
  
  emb <- Embeddings(object, reduction = reduction)[, dims, drop = FALSE]
  groups <- object@meta.data[[group.by]]
  
  out <- mixing_metric(
    X = emb,
    groups = groups,
    k = k,
    metric = metric,
    return_per_cell = return_per_cell
  )
  
  if (return_per_cell) {
    nm <- paste0("mixing_", metric)
    object[[nm]] <- out$per_cell
    out$column_added <- nm
    out$seurat_object <- object
  }
  
  return(out)
}



get_sample_drivers <- function(object, sample_id, pc = "PC_1", n_genes = 50) {
  # 1. Subset to the specific sample
  s_obj <- subset(object, subset = sample.id == !!sample_id)
  
  # 2. Get global loadings for the chosen PC
  loadings <- Loadings(object, reduction = "pca")[, pc]
  
  # 3. Calculate gene variance within THIS sample only
  # We only look at genes that were part of the PCA
  scaled_data <- GetAssayData(s_obj, slot = "scale.data")
  relevant_genes <- intersect(names(loadings), rownames(scaled_data))
  
  sample_var <- apply(scaled_data[relevant_genes, ], 1, var)
  
  # 4. Weight the loading by the local variance
  # This highlights genes that have high weight AND are varying in this patient
  driver_score <- loadings[relevant_genes] * sample_var
  
  top_pos <- sort(driver_score, decreasing = TRUE)[1:n_genes]
  top_neg <- sort(driver_score, decreasing = FALSE)[1:n_genes]
  
  return(list(
    sample = sample_id,
    positive_drivers = names(top_pos),
    negative_drivers = names(top_neg)
  ))
}


