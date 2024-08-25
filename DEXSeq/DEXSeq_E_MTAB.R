# LEBIB Inès - IdRS - 02/07/2023

# Le but de ce script est de tester l'outil DEXSeq avec le jeu de données E_MTAB_11609

# Le package Bioconductor DEXSeq est une méthode pour quantifier les 
# événements d'épissage alternatif au niveau exoniques.

# Liens utiles pour la documentation de DEXSeq : 
# https://bioconductor.org/packages/release/bioc/html/DEXSeq.html
# https://bioconductor.org/packages/release/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html
#https://rdrr.io/bioc/DEXSeq/f/vignettes/DEXSeq.Rmd


# Lien du jeu de données
# https://www.ncbi.nlm.nih.gov/Traces/study/?acc=E-MTAB-11609&o=acc_s%3Aa

# Le jeu de données que nous allons analyser est E_MTAB_11609
# Il s'agit d'un modèle de résistance acquis au vémurafénib dans 
# le mélanome malin (lignée SK-MEL-239, présentant la mutation BRAF V600E). 

# Séquençage de l’ARN dans des cellules parentales et clones C3 SK-MEL-239

# 4 réplicats biologiques de cellules parentales 
# 4 réplicats biologiques de cellules résistantes - clones C3 SK-MEL-239 
# 30 millions de reads PE de 150 pb par échantillon, 8 échantillons 


install.packages('knitr')

knitr::opts_chunk$set(tidy = FALSE,
                      cache = TRUE,
                      dev = "png",
                      message = FALSE,
                      error = FALSE,
                      warning = TRUE)

 library(reticulate)

 BiocManager::install("BiocStyle")
 BiocStyle::Biocpkg("DEXSeq")
 BiocStyle::Biocpkg("Rsubread")
 BiocStyle::Biocpkg("DESeq2")




#### Etape 1 : 
# Alignement des reads sur un génome de référence

# Installation du packgae HTSeq sur Python 
# https://pypi.org/project/HTSeq/
# /!\ Importer HTSeq dans le même répertoire que votre script R


# Telechargement et chemin des scripts python
pythonScriptsDir =  system.file ( "python_scripts" , package = "DEXSeq" )
list.files (pythonScriptsDir)
system.file ( "python_scripts" , package = "DEXSeq" , mustWork = TRUE )




# Preparation de l'annotation 
# Le script python 'dexseq_prepare_annotation.py' va 'réduire' le 
# fichier GTF en GFF en définissant des bacs de comptage d'exons
# input :  GTF , output : GFF

# /!\ Répertoire actuel doit contenir fichier GTF
# Appelez script python

# Ligne de commande à lancer à partir du shell, pas R 
# exemple : python /path/to/library/DEXSeq/python_scripts/dexseq_prepare_annotation.py Drosophila_melanogaster.BDGP5.72.gtf Dmel.BDGP5.25.62.DEXSeq.chr.gff
getwd()
setwd('~/Decouverte/Test_DEXSeq')

# python /home/AJ84_SUR/R/IMG_R_4.1/4.2.1/DEXSeq/python_scripts/dexseq_prepare_annotation_Gene_name.py GRCh38.99.Chr.gtf GRCh38.99.Chr_Gene_name.gff



#####
# Compter des reads
# input : GFF, BAM 
# output : txt
#####

# python /home/AJ84_SUR/R/IMG_R_4.1/4.2.1/DEXSeq/python_scripts/dexseq_count.py -p yes -r pos -f bam GRCh38.99.Chr_Gene_name.gff E_MTAB_11609/E_MTAB_11609/Parental_1.bam ~/Decouverte/Test_DEXSeq/E_MTAB_Count_Reads_Gene_name/Parental_1.txt
# python /home/AJ84_SUR/R/IMG_R_4.1/4.2.1/DEXSeq/python_scripts/dexseq_count.py -p yes -r pos -f bam GRCh38.99.Chr_Gene_name.gff E_MTAB_11609/E_MTAB_11609/Parental_2.bam ~/Decouverte/Test_DEXSeq/E_MTAB_Count_Reads_Gene_name/Parental_2.txt
# python /home/AJ84_SUR/R/IMG_R_4.1/4.2.1/DEXSeq/python_scripts/dexseq_count.py -p yes -r pos -f bam GRCh38.99.Chr_Gene_name.gff E_MTAB_11609/E_MTAB_11609/Parental_3.bam ~/Decouverte/Test_DEXSeq/E_MTAB_Count_Reads_Gene_name/Parental_3.txt
# python /home/AJ84_SUR/R/IMG_R_4.1/4.2.1/DEXSeq/python_scripts/dexseq_count.py -p yes -r pos -f bam GRCh38.99.Chr_Gene_name.gff E_MTAB_11609/E_MTAB_11609/Parental_4.bam ~/Decouverte/Test_DEXSeq/E_MTAB_Count_Reads_Gene_name/Parental_4.txt
# python /home/AJ84_SUR/R/IMG_R_4.1/4.2.1/DEXSeq/python_scripts/dexseq_count.py -p yes -r pos -f bam GRCh38.99.Chr_Gene_name.gff E_MTAB_11609/E_MTAB_11609/Vem_Resistant_1.bam ~/Decouverte/Test_DEXSeq/E_MTAB_Count_Reads_Gene_name/Vem_Resistant_1.txt
# python /home/AJ84_SUR/R/IMG_R_4.1/4.2.1/DEXSeq/python_scripts/dexseq_count.py -p yes -r pos -f bam GRCh38.99.Chr_Gene_name.gff E_MTAB_11609/E_MTAB_11609/Vem_Resistant_2.bam ~/Decouverte/Test_DEXSeq/E_MTAB_Count_Reads_Gene_name/Vem_Resistant_2.txt
# python /home/AJ84_SUR/R/IMG_R_4.1/4.2.1/DEXSeq/python_scripts/dexseq_count.py -p yes -r pos -f bam GRCh38.99.Chr_Gene_name.gff E_MTAB_11609/E_MTAB_11609/Vem_Resistant_3.bam ~/Decouverte/Test_DEXSeq/E_MTAB_Count_Reads_Gene_name/Vem_Resistant_3.txt


# python /home/AJ84_SUR/R/IMG_R_4.1/4.2.1/DEXSeq/python_scripts/dexseq_count.py -p yes -r pos -f bam GRCh38.99.Chr_Gene_name.gff E_MTAB_11609/E_MTAB_11609/Vem_Resistant_4.bam ~/Decouverte/Test_DEXSeq/E_MTAB_Count_Reads_Gene_name/Vem_Resistant_4.txt

getwd()
# gsutil -m cp -r ~/Decouverte/Test_DEXSeq/E_MTAB_Count_Reads_Gene_name/Parental_1_clean.txt gs://em52-technology-development-prod-d74c/Splice_Alternatif_Tool_Evaluation/DataSet/E_MTAB_11609/Count_Reads

# /!\ Si l'output txt contient des guillemets, il faut les retier avec cette commande : 



# sed 's/\"//g' Parental_1.txt > Parental_1_clean.txt
# sed 's/\"//g' Parental_2.txt > Parental_2_clean.txt
# sed 's/\"//g' Parental_3.txt > Parental_3_clean.txt
# sed 's/\"//g' Parental_4.txt > Parental_4_clean.txt


# sed 's/\"//g' Vem_Resistant_1.txt > Vem_Resistant_1_clean.txt
# sed 's/\"//g' Vem_Resistant_2.txt > Vem_Resistant_2_clean.txt
# sed 's/\"//g' Vem_Resistant_3.txt > Vem_Resistant_3_clean.txt
# sed 's/\"//g' Vem_Resistant_4.txt > Vem_Resistant_4_clean.txt

######
# Memo sur parallelisation
######


# Méthodes de parallélisation avec BiocParallel avec implémentations
# de BPPARAM= avec les fonctions estimateDispersions(), testForDEU() 
# et estimateExonFoldChanges()


# Le reste de l'analyse se fait maintenant sur R
# En utilisant les ouput des deux scripts python


countFiles =  list.files (path = '~/Decouverte/Test_DEXSeq/E_MTAB_Count_Reads', pattern = "*_clean.txt$" , full.names = TRUE )
basename(countFiles)

flattenedFile =  list.files (path = '/home/AJ84_SUR/Decouverte/Test_DEXSeq', pattern = "Chr.gff$" , full.names = TRUE )
basename(flattenedFile)


# Maintenant que l'on a les fichiers gff et txt, on va préparer le tableau
# Ce tableau doit contenir :
#   - une ligne pour chaque bibliotheque
#   - des colonnes contenant des informations pertinantes (nom du fichier, conditions expérimentales, informations techniques ...)

# condition : controle experimentale avec niveau de traitement
# libType : single-end, paired-end ...

sampleTable = data.frame(
 row.names = c( "Parental_1", "Parental_2", "Parental_3", "Parental_4",
                  "Vem_Resistant_1", "Vem_Resistant_2", "Vem_Resistant_3", "Vem_Resistant_4" ),
 
 condition = c( "Parental", "Parental", "Parental", "Parental",
                "Vem_Resistant", "Vem_Resistant", "Vem_Resistant", "Vem_Resistant" ),
 
 libType = c( "paired-end", "paired-end", "paired-end", "paired-end", 
              "paired-end", "paired-end", "paired-end", "paired-end" ) )



#####
# Construction de l'objet DEXSeqDataSet
# avec la fonction DEXSeqDataSetFromHTSeq()
#####

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DEXSeq")
library( "DEXSeq" )


#####
# La fonction DEXSeqDataSetFromHTSeq prend 4 arguments: 
#   - un vecteur avec noms de fichiers count (output de dexseq_count.py)
#   - sample table
#   - formule : ~ sample + exon + condition:exon
# la formule permet de s'interesser aux différences d'epissage des exons
#   - GFF  
#####


# Chaque ligne de la classe DEXSeqDataSet contient dans chaque colonne des
# données de comptage d'un exon donné  (« this ») et les données de comptage
# de la somme des autres exons appartenant au même gène (« others »)
# Toutes ces informations sont rappelées dans la colonne colData()

dxd = DEXSeqDataSetFromHTSeq(
 countFiles,
 sampleData=sampleTable,
 design= ~ sample + exon + condition:exon,
  flattenedfile=flattenedFile )

# Les visualisations de dxd

# colData(dxd)
# head( rowRanges(dxd), 3 )
# head( counts(dxd), 5 )
# split( seq_len(ncol(dxd)), colData(dxd)$exon )
# head( featureCounts(dxd), 5 )
# head( rowRanges(dxd), 3 )
# sampleAnnotation( dxd )

#####
# Normalisation
#####

dxd = estimateSizeFactors( dxd )

# Visualisation de dxd normalisé 

# colData(dxd)
# head( rowRanges(dxd), 3 )
# head( counts(dxd), 5 )
# split( seq_len(ncol(dxd)), colData(dxd)$exon )
# head( featureCounts(dxd), 5 )
# head( rowRanges(dxd), 3 )
# sampleAnnotation( dxd )


# Estimation de la dispersion
# càd estimation de la variabilité des données 
# Permet de distinguer les variations techniques et biologiques (bruit)
# des effets reels d'evenement d'épissage


# Calcul dispersion par exon : 
#   - estimation de la vraisemblance de profil ajustée de Cox-Reid
#   - relation dispersion-moyenne est ajustée aux valeurs de dispersion individuelles
#   - valeurs ajustées sont prises pour réduire les estimations par exon


# Visualisation de la dispersion par exon par rapport 
# au nombre moyen normalisé, aux valeurs ajustées résultantes et aux 
# estimations de dispersion 

# plotDispEsts(dxd)


######
# Test de l'utilisation différentielle des exons
######

# La fonction va effectuer des tests pour chaque exon de chaque gène
# Ajuster un modèle linéaire généralisé et 
# comparer au modèle plus petit (modèle nul)

# ~sample + exon + condition:exon~ sample + exon


# On supprime libType car une seule variable unique : paired - end 
# Variable ne va pas être utile pour calculer le modèle linéaire

# A la place de : 
# formulaFullModel    =  ~ sample + exon + type:exon + condition:exon
# formulaReducedModel =  ~ sample + exon + type:exon 

# On va écrire : 

formulaFullModel    =  ~ sample + exon + condition:exon
formulaReducedModel =  ~ sample + exon 

# Regarder BPPARAM pour paralleliser 
# dans estimateDispersions

dxd = estimateDispersions( dxd, formula = formulaFullModel )


dxd = testForDEU( dxd, 
                  reducedModel = formulaReducedModel, 
                  fullModel = formulaFullModel )


dxr2 = DEXSeqResults( dxd )

table( dxr2$padj < 0.05 )
#  FALSE   TRUE 
# 341182   6179 

table( dxr2$pvalue < 0.05 )

indices <- which( dxr2$pvalue < 0.05 )

dxr2[indices]

DEXSeqHTML( dxr2, FDR=0.05, color=c("#FF000080", "#0000FF80") )

############## point sophie
### Regarder si on peut avoir nom du gène au lieu du gène ID - ENSG
### Se renseigner sur fichier html et métrique que l'on peut ajouter
### DONNER LES LIMITES DE DEXSeq
#############





### Pour avoir les gènes name de tous les gènes ID  

dxr3 <- dxr2
dxr3 <- data.frame(dxr3)

data1  <- merge(genes, dxr3,  by="groupID")
data2 <- unique(data1[, 1:2])


## Controle du FDR au niveau des gènes

numbOfGenes <- sum( perGeneQValue(dxr2) < 0.05)
numbOfGenes

data_2 <- data.frame(perGeneQValue(dxr2) < 0.05)
data_3 <- subset(data_2, perGeneQValue.dxr2....0.05 == 'TRUE')



#######
# Visualisation
#######


# Pour commencer, on reprend les données de l'article, 
# là où on a déjà des résultats
 
 
  # ID Ensembl BRAF
# plotDEXSeq( dxr2, "ENSG00000157764", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
# plotDEXSeq( dxr2, "ENSG00000157764", displayTranscripts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq( dxr4, "ENSG00000157764", expression=FALSE, splicing=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )


###
# REGARDER A QUOI CORRESPOND CES EXONS - DOMAINES Prot ? 
###

  # ID Ensembl EVI5L - perte exon 12 dans article 
# plotDEXSeq( dxr2, "ENSG00000142459", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
# plotDEXSeq( dxr2, "ENSG00000142459", displayTranscripts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq( dxr2, "ENSG00000142459", expression=FALSE, splicing=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )


  # ID Ensembl TYR - perte exon 3 dans article
# plotDEXSeq( dxr2, "ENSG00000077498", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
# plotDEXSeq( dxr2, "ENSG00000077498", displayTranscripts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq( dxr2, "ENSG00000077498", expression=FALSE, splicing=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )


  # ID Ensembl CAPN3 - perte exon 15 dans article 
# plotDEXSeq( dxr2, "ENSG00000092529", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
# No read counts falling in this gene, there is nothing to plot.


  # ID Ensembl CLSTN1 - perte exon 11 et 3 dans article
# plotDEXSeq( dxr2, "ENSG00000171603", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
# plotDEXSeq( dxr2, "ENSG00000171603", displayTranscripts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq( dxr2, "ENSG00000171603", expression=FALSE, splicing=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )


  # ID Ensembl FANCA  
# plotDEXSeq( dxr2, "ENSG00000187741", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
# plotDEXSeq( dxr2, "ENSG00000187741", displayTranscripts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq( dxr2, "ENSG00000187741", expression=FALSE, splicing=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )



  # ID Ensembl EPB41 - perte exon 15 dans article
# plotDEXSeq( dxr2, "ENSG00000159023", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
# plotDEXSeq( dxr2, "ENSG00000159023", displayTranscripts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq( dxr2, "ENSG00000159023", expression=FALSE, splicing=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )


  # ID Ensembl MARK3 - perte exon 16 dans article 
# plotDEXSeq( dxr2, "ENSG00000075413", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
# plotDEXSeq( dxr2, "ENSG00000075413", displayTranscripts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq( dxr2, "ENSG00000075413", expression=FALSE, splicing=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )


  # ID Ensembl PICALM  
# plotDEXSeq( dxr2, "ENSG00000073921", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
# plotDEXSeq( dxr2, "ENSG00000073921", displayTranscripts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq( dxr2, "ENSG00000073921", expression=FALSE, splicing=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )



  # ID Ensembl SYNE2 
# plotDEXSeq( dxr2, "ENSG00000054654", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
# plotDEXSeq( dxr2, "ENSG00000054654", displayTranscripts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq( dxr2, "ENSG00000054654", expression=FALSE, splicing=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )


  # ID Ensembl MPRIP 
# plotDEXSeq( dxr2, "ENSG00000133030", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
# plotDEXSeq( dxr2, "ENSG00000133030", displayTranscripts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq( dxr2, "ENSG00000133030", expression=FALSE, splicing=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )



  # ID Ensembl FAM126A 
# plotDEXSeq( dxr2, "ENSG00000122591", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
# plotDEXSeq( dxr2, "ENSG00000122591", displayTranscripts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq( dxr2, "ENSG00000122591", expression=FALSE, splicing=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )


  # ID Ensembl CHID1 
# plotDEXSeq( dxr2, "ENSG00000177830", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
# plotDEXSeq( dxr2, "ENSG00000177830", displayTranscripts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )
plotDEXSeq( dxr2, "ENSG00000177830", expression=FALSE, splicing=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )







# save(dxd, file = 'dxd_E_MTAB_DEXSeq.RData')
# save(dxr1, file = 'dxr1_E_MTAB_DEXSeq.RData')
# save(dxr2, file = 'dxr2_E_MTAB_DEXSeq.RData')



####
# On suit le complément de materiels et méthodes
#### 

getwd()

setwd('~/Decouverte/Test_DEXSeq/E_MTAB_Count_Reads')

load("~/Decouverte/Test_DEXSeq/E_MTAB_Count_Reads/dxd_E_MTAB_DEXSeq.RData")

dxd = estimateExonFoldChanges( dxd, fitExpToVar="condition")


# extract results and write into table:
dxr1 = DEXSeqResults( dxd )
dfdxr1 =data.frame (dxr1$groupID, dxr1$pvalue, dxr1$padj, dxr1$log2fold_Vem_Resistant_Parental, dxr1$genomicData )

write.table (dfdxr1, file = "DEXSeq-table.txt", append = FALSE, quote = TRUE,
             sep = "\t",
             eol = "\n", na = "NA", dec = ".", row.names = FALSE,
             col.names = TRUE)


setwd('~/Decouverte/Test_DEXSeq/E_MTAB_Count_Reads/DEXSeq_Report_Article')


DEXSeqHTML( dxr1, file="DEXSec-wholegenome.html", FDR=0.01,
            color=c("#FF000080", "#00ff00") )


BiocManager::install("biomaRt")


DEXSeqHTML(dxr1,fitExpToVar="condition",FDR=0.05, 
           color=c("#FF000080", "#0000FF80"),
           mart=ensembl,
           filter="ensembl_gene_id",
           attributes=c("external_gene_id","description"))

















