# rm(list = ls())
 

# Script pour SpliceWiz

# Installation

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("SpliceWiz")

install.packages("DoubleExpSeq")
BiocManager::install(c("DESeq2", "limma", "edgeR"))

# Chargement de SpliceWiz

library(SpliceWiz)

# Interface utilisateur graphique (GUI) de SpliceWiz

if(interactive()) {
  spliceWiz(demo = TRUE)
}


# Construire la référence de SpliceWiz

# La référence SpliceWiz est utilisée pour quantifier l'épissage alternatif dans les fichiers BAM, ainsi que dans le classement en aval, l'analyse différentielle et la visualisation.


ref_path <- file.path(tempdir(), "Reference")
buildRef(
  reference_path = ref_path,
  fasta = "/home/AJ84_SUR/Decouverte/Test_DEXSeq/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
  gtf = "/home/AJ84_SUR/Decouverte/Test_DEXSeq/GRCh38.99.Chr.gtf",
  ontologySpecies = "Homo sapiens"
)


# Traiter les fichiers BAM à l'aide de SpliceWiz

pb_path <- file.path(tempdir(), "pb_output")

bam_path <- "~/Decouverte/Test_DEXSeq/E_MTAB_11609"

bam_files <- list.files(path = bam_path, pattern = "\\.bam$", full.names = TRUE)

Sample_names <- gsub(pattern = ".bam$", replacement = "", x = basename(bam_files))


processBAM(
  bamfiles = bam_files,
  sample_names = Sample_names,
  reference_path = ref_path,
  output_path = pb_path
)


# Rassembler d'expérience

expr <- findSpliceWizOutput(pb_path)


nxtse_path <- file.path(tempdir(), "NxtSE_output")
collateData(
  Experiment = expr,
  reference_path = ref_path,
  output_path = nxtse_path, 
  n_threads = 4
)


# Importer l'epérience

se <- makeSE(nxtse_path)


# Analyse différentielle
# Affectations d'annotation aux échantillons


colData(se)$condition <- rep(c("Parentales", "Vem_Resistant"), each = 4)
colData(se)$batch <- rep(c("Condition1", "Condition2"), each = 4)


# Filtrage des événements de haute confiance

se.filtered <- se[applyFilters(se),]

# Effectuer une analys différentielle

require("edgeR")
res_edgeR <- ASE_edgeR(
  se = se.filtered,
  test_factor = "condition",
  test_nom = "Vem_Resistant",
  test_denom = "Parentales"
)



# Visualisation 
# Volcano plots


library(ggplot2)

ggplot(res_edgeR,
       aes(x = logFC, y = -log10(FDR))) + 
  geom_point() +
  labs(title = "Differential analysis - Resistant vs Parentale",
       x = "Log2-fold change", y = "FDR (-log10)")


ggplot(res_edgeR,
       aes(x = logFC, y = -log10(FDR))) + 
  geom_point() + facet_wrap(vars(EventType)) +
  labs(title = "Differential analysis - Resistant vs Parentale",
       x = "Log2-fold change", y = "FDR (-log10)")


library(stringr)

new_res_edgeR <- res_edgeR

new_res_edgeR$EventLabel <- sapply(str_split(res_edgeR$EventName, ":"), function(x) x[2])



ggplot(new_res_edgeR, aes(x = logFC, y = -log10(FDR))) + 
  geom_point() + facet_wrap(vars(EventType)) +
  geom_text(aes(label = EventLabel), vjust = "inward", hjust = "inward", check_overlap = TRUE, size = 2.5) + 
  labs(title = "Differential analysis - Resistant vs Parentale",
       x = "Log2-fold change", y = "FDR (-log10)")


save.image(file = '/home/AJ84_SUR/test_SpliceWiz/SpliceWiz.RData')


# Scatter plots

ggplot(res_edgeR, aes(x = 100 * AvgPSI_B, y = 100 * AvgPSI_A)) + 
  geom_point() + xlim(0, 100) + ylim(0, 100) +
  labs(title = "PSI values across conditions",
       x = "PSI of condition Resistant", y = "PSI of condition Parentale")


if("res_edgeR" %in% ls() && "AvgPSI_B" %in% names(res_edgeR)) {
  print("La variable est présente et prête à être utilisée.")
} else {
  print("La variable n'est pas présente. Vérifiez votre dataframe ou chargez les données nécessaires.")
}



# save.image(file = '/home/AJ84_SUR/test_SpliceWiz/SpliceWiz.RData')



res_edgeR <- ASE_edgeR(
  se = se.filtered,
  test_factor = "condition",
  test_nom = "Vem_Resistant",
  test_denom = "Parentales"
)
 #> Jan 05 19:25:25 Performing edgeR contrast for included / excluded counts separately
#> Jan 05 19:25:26 Performing edgeR contrast for included / excluded counts together

res_edgeR.filtered <- res_edgeR[res_edgeR$abs_deltaPSI > 0.05,]
res_edgeR.filtered$EventName[1]
#> [1] "NSUN5/ENST00000252594_Intron2/clean"


save.image(file = '/home/AJ84_SUR/test_SpliceWiz/SpliceWiz.RData')


# On veut extraire l'événement pour PICALM 

picalm_events <- res_edgeR.filtered[grep("PICALM", res_edgeR.filtered$EventName),]


dataObj <- getCoverageData(se, Gene = "PICALM", tracks = colnames(se))

plotObj <- getPlotObject(dataObj, Event = picalm_events$EventName[1])


plotView(
  plotObj,
  centerByEvent = TRUE,  # whether the plot should be centered at the 'Event'
  trackList = list(1,2,3,4,5,6,7,8),
  plotJunctions = TRUE
)


# BRAF

braf_events <- res_edgeR.filtered[grep("BRAF", res_edgeR.filtered$EventName),]


dataObj_BRAF <- getCoverageData(se, Gene = "BRAF", tracks = colnames(se))

plotObj_BRAF <- getPlotObject(dataObj_BRAF, Event = braf_events$EventName[1])


plotView(
  plotObj_BRAF,
  centerByEvent = TRUE,  # whether the plot should be centered at the 'Event'
  trackList = list(1,2,3,4,5,6,7,8),
  plotJunctions = TRUE
)




# EVI5L


EVI5L_events <- res_edgeR.filtered[grep("EVI5L", res_edgeR.filtered$EventName),]


dataObj_EVI5L <- getCoverageData(se, Gene = "EVI5L", tracks = colnames(se))

plotObj_EVI5L <- getPlotObject(dataObj_EVI5L, Event = EVI5L_events$EventName[1])


plotView(
  plotObj_EVI5L,
  centerByEvent = TRUE,  # whether the plot should be centered at the 'Event'
  trackList = list(1,2,3,4,5,6,7,8),
  plotJunctions = TRUE
)


# TYR

TYR_events <- res_edgeR.filtered[grep("TYR", res_edgeR.filtered$EventName),]


dataObj_TYR <- getCoverageData(se, Gene = "TYR", tracks = colnames(se))

plotObj_TYR <- getPlotObject(dataObj_TYR, Event = TYR_events$EventName[1])


plotView(
  plotObj_TYR,
  centerByEvent = TRUE,  # whether the plot should be centered at the 'Event'
  trackList = list(1,2,3,4,5,6,7,8),
  plotJunctions = TRUE
)



# CLSTN1


CLSTN1_events <- res_edgeR.filtered[grep("CLSTN1", res_edgeR.filtered$EventName),]


dataObj_CLSTN1 <- getCoverageData(se, Gene = "CLSTN1", tracks = colnames(se))

plotObj_CLSTN1 <- getPlotObject(dataObj_CLSTN1, Event = CLSTN1_events$EventName[1])


plotView(
  plotObj_CLSTN1,
  centerByEvent = TRUE,  # whether the plot should be centered at the 'Event'
  trackList = list(1,2,3,4,5,6,7,8),
  plotJunctions = TRUE
)




# EPB41


EPB41_events <- res_edgeR.filtered[grep("EPB41", res_edgeR.filtered$EventName),]


dataObj_EPB41 <- getCoverageData(se, Gene = "EPB41", tracks = colnames(se))

plotObj_EPB41 <- getPlotObject(dataObj_EPB41, Event = EPB41_events$EventName[1])


plotView(
  plotObj_EPB41,
  centerByEvent = TRUE,  # whether the plot should be centered at the 'Event'
  trackList = list(1,2,3,4,5,6,7,8),
  plotJunctions = TRUE
)




# MARK3

MARK3_events <- res_edgeR.filtered[grep("MARK3", res_edgeR.filtered$EventName),]


dataObj_MARK3 <- getCoverageData(se, Gene = "MARK3", tracks = colnames(se))

plotObj_MARK3 <- getPlotObject(dataObj_MARK3, Event = MARK3_events$EventName[1])


plotView(
  dataObj_MARK3,
  centerByEvent = TRUE,  # whether the plot should be centered at the 'Event'
  trackList = list(1,2,3,4,5,6,7,8),
  plotJunctions = TRUE
)










