# Test de l'outil ASPli 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ASpli")
BiocManager::install("GenomeInfoDb")
BiocManager::install("GenomeInfoDb", force = TRUE)

BiocManager::install("GenomicRanges")



library(ASpli)
library(GenomicFeatures)

BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")

# gtf preprocessing
gtfFileName <- '/home/AJ84_SUR/test_ASPli/GRCh38.99.Chr.gtf'
library(rtracklayer)
library(GenomicFeatures)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GenomicFeatures")


genomeTxDb <- makeTxDbFromGFF( gtfFileName )

library(GenomicFeatures)


saveDb(genomeTxDb,file="gene.sqlite")

# feature extraction ----
features <- binGenome( genomeTxDb )

#bams and target file ----
BAMFiles <- c('/home/AJ84_SUR/Decouverte/Test_DEXSeq/E_MTAB_11609/Parental_1.bam',
              '/home/AJ84_SUR/Decouverte/Test_DEXSeq/E_MTAB_11609/Parental_2.bam',
              '/home/AJ84_SUR/Decouverte/Test_DEXSeq/E_MTAB_11609/Parental_3.bam',
              '/home/AJ84_SUR/Decouverte/Test_DEXSeq/E_MTAB_11609/Parental_4.bam', 
              '/home/AJ84_SUR/Decouverte/Test_DEXSeq/E_MTAB_11609/Vem_Resistant_1.bam',
              '/home/AJ84_SUR/Decouverte/Test_DEXSeq/E_MTAB_11609/Vem_Resistant_2.bam',
              '/home/AJ84_SUR/Decouverte/Test_DEXSeq/E_MTAB_11609/Vem_Resistant_3.bam',
              '/home/AJ84_SUR/Decouverte/Test_DEXSeq/E_MTAB_11609/Vem_Resistant_4.bam')

row_names <- c( 'Parental_1','Parental_2','Parental_3','Parental_4','Vem_Resistant_1','Vem_Resistant_2','Vem_Resistant_3','Vem_Resistant_4')

targets <- data.frame(row.names = row_names,
                      bam = BAMFiles[1:8],
                      f1 = c( 'Parental','Parental','Parental','Parental','Vem_Resistant','Vem_Resistant','Vem_Resistant','Vem_Resistant'),
                      stringsAsFactors = FALSE)


mBAMs <- data.frame( bam = targets$bam[c(1:8)],
                     condition = c( 'Parental','Parental','Parental','Parental','Vem_Resistant','Vem_Resistant','Vem_Resistant','Vem_Resistant'))
# gbcounts: Permet de lire le comptage par rapport aux fonctionnalités annotées
gbcounts <- gbCounts(features=features, targets=targets,
                     minReadLength = 90, maxISize = 50000)


# asd : Comptage de novo et estimation du signal d épissage basé sur les jonctions
# jCounts : Résumer les indices d'inclusion des jonctions PSI, PIR et PJU
asd <- jCounts(counts=gbcounts, features=features, minReadLength=90)

# gb : Expression différentielle des gènes et estimation du signal d’utilisation des bacs
gb <- gbDUreport(gbcounts, contrast = c(-1,1))

# jdur : Analyse de l'utilisation des jonctions différentielles
jdur <- jDUreport(asd, contrast=c(-1,1))

# sr : Intégration des signaux de bac et de jonction
sr <- splicingReport(gb, jdur, counts=gbcounts)

# is : Résumé de l'intégration des signaux d'épissage le long des régions génomiques
is <- integrateSignals(sr,asd)

# Exportation des résultats
exportIntegratedSignals(is,sr=sr,
                        output.dir = "Resultat_ASPli_EMTAB",
                        counts=gbcounts,features=features,asd=asd,
                        mergedBams = mBAMs)

# Exporter des tableaux vers des fichiers texte
writeAS(as=asd, output.dir="Resultat_tableau_ASPli_EMTAB")


 # gsutil -m cp -r ~/Resultat_tableau_ASPli_EMTAB gs://em52-technology-development-prod-src-explo-prod-e88c

# samtools index /home/AJ84_SUR/Decouverte/Test_DEXSeq/E_MTAB_11609/Parental_1.bam
# samtools index /home/AJ84_SUR/Decouverte/Test_DEXSeq/E_MTAB_11609/Parental_2.bam
# samtools index /home/AJ84_SUR/Decouverte/Test_DEXSeq/E_MTAB_11609/Parental_3.bam
# samtools index /home/AJ84_SUR/Decouverte/Test_DEXSeq/E_MTAB_11609/Parental_4.bam
# samtools index /home/AJ84_SUR/Decouverte/Test_DEXSeq/E_MTAB_11609/Vem_Resistant_1.bam
# samtools index /home/AJ84_SUR/Decouverte/Test_DEXSeq/E_MTAB_11609/Vem_Resistant_2.bam
# samtools index /home/AJ84_SUR/Decouverte/Test_DEXSeq/E_MTAB_11609/Vem_Resistant_3.bam
# samtools index /home/AJ84_SUR/Decouverte/Test_DEXSeq/E_MTAB_11609/Vem_Resistant_4.bam





