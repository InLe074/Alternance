##### LEBIB Inès - IDRS #####

### **Introduction** ###

L’épissage alternatif est l’épissage de différents exons à partir de l’ARN pré-messager. Il permet de diversifier les séquences des protéines ainsi que leurs fonctions et leurs régulations. 

Le package R/Bioconductor DEXSeq est une méthode pour quantifier les événements d'épissage alternatif au niveau exonique.


Liens utiles pour la documentation de DEXSeq : 

* https://bioconductor.org/packages/release/bioc/html/DEXSeq.html

* https://bioconductor.org/packages/release/bioc/vignettes/DEXSeq/inst/doc/DEXSeq.html

* https://rdrr.io/bioc/DEXSeq/f/vignettes/DEXSeq.Rmd


### **Packages et librairies à installer** ###

```
install.packages('knitr')

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
  
  
BiocManager::install("BiocStyle")
BiocStyle::Biocpkg("DEXSeq")
BiocStyle::Biocpkg("Rsubread")
BiocStyle::Biocpkg("DESeq2")

library(reticulate)
library( "DEXSeq" )

knitr::opts_chunk$set(tidy = FALSE,
                      cache = TRUE,
                      dev = "png",
                      message = FALSE,
                      error = FALSE,
                      warning = TRUE)
```


### **Etape 1 : Comptage** ###

#### **HTSeq Count** ####

Les premières étapes de l'analyse sont effectuées à l'aide de deux scripts Python fournit avec le package HTSeq.

Installation du package [HTSeq](https://pypi.org/project/HTSeq/) sur Python.

`pip install HTSeq`

**Important** : il faut importer le package HTSeq, **toujours depuis le script python**, dans le **même** répertoire que votre **script R** avec la commande `cd` puis `import HTSeq`

Une fois le package installé, vous pouvez utiliser les deux scripts Python *dexseq_prepare_annotation.py* et *dexseq_count.py*


#### **Téléchargement et chemin des scripts python dans le script R** ####

```
pythonScriptsDir =  system.file ( "python_scripts" , package = "DEXSeq" )
list.files (pythonScriptsDir)
system.file ( "python_scripts" , package = "DEXSeq" , mustWork = TRUE )
```


#### **Preparation de l'annotation** ####
Le script python *dexseq_prepare_annotation.py* va 'réduire' le fichier GTF en fichier GFF. 

input :  GTF , output : GFF

Il faut s'assurer que le répertoire actuel contient le fichier GTF. 

Le but est d'appeler le script python en lançant la ligne de commande ci-dessous à partir du **shell**, et non la console R. Le chemin des scripts python est indiqué dans le `system.file`


`python /path/to/library/DEXSeq/python_scripts/dexseq_prepare_annotation.py Drosophila_melanogaster.BDGP5.72.gtf Dmel.BDGP5.25.62.DEXSeq.chr.gff`


#### **Compter les reads** ####

input : BAM et GFF, output : txt

Pour chaque fichier BAM, les reads qui chevauchent les exons vont être comptés. 

Il faut lancer cette commande en une seule fois, dans le **shell**, pour chaque fichier BAM.

`python /path/to/library/DEXSeq/python_scripts/dexseq_count.py Dmel.BDGP5.25.62.DEXSeq.chr.gff untreated1.bam untreated1fb.txt`

- **Paired-end data : ** Si les données sont paired-end, il faut ajouter l'option `-p yes` et si elles sont triées par position il faut ajouter `-r pos`. 

- **Fichiers SAM et BAM : ** par défaut, le script s'attend à avoir en input un fichier SAM mais il peut aussi lire les fichiers BAM. Il faudra cependant ajouter l'option `-f bam`. 

- **Qualité de l'alignement : ** Il y a la possibilité d'ajouter une option supplémentaire pour spécifier la qualité d'alignement minimale. Tous les reads dont la qualité est inférieur à celle spécifiée seront ignorés `-a-a 10`

Ces options sont à ajouter après l'appel du fichier python *dexseq_count.py*



Important : si l'output contient des guillemets, il faut les retier avec cette commande, toujours dans le **shell**, pour chaque fichier txt : 

`sed 's/\"//g' untreated1fb.txt > untreated1fb_clean.txt`

### **Analyse des données dans un script R** ###

Le reste de l'analyse se fait maintenant avec un script R en utilisant les ouput des deux scripts python. 

```
countFiles = list.files(inDir, pattern="fb.txt$", full.names=TRUE)
basename(countFiles)
flattenedFile = list.files(inDir, pattern="gff$", full.names=TRUE)
basename(flattenedFile)
```

Maintenant que l'on a les fichiers GFF et txt, on va préparer le tableau `sampleTable`.

Ce tableau doit contenir :

   - une ligne pour chaque échantillon
   
   - des colonnes contenant des informations pertinantes (nom du fichier, conditions expérimentales, informations techniques ...)


*condition* : contrôle experimental avec niveau de traitement

*libType* : le type de bibliothèque : single-end, paired-end ...


```
sampleTable = data.frame(
   row.names = c( "treated1", "treated2", "treated3", 
      "untreated1", "untreated2", "untreated3", "untreated4" ),
   condition = c("knockdown", "knockdown", "knockdown",  
      "control", "control", "control", "control" ),
   libType = c( "single-end", "paired-end", "paired-end", 
      "single-end", "single-end", "paired-end", "paired-end" ) )
```


#### **Construction de l'objet DEXSeqDataSet avec la fonction DEXSeqDataSetFromHTSeq()** #### 

La fonction DEXSeqDataSetFromHTSeq prend 4 arguments: 

 - un vecteur avec noms de fichiers count (output de *dexseq_count.py*)
 
 - `sampleTable`
 
 - formule : `~ sample + exon + condition:exon`
 
 - GFF (output de *dexseq_prepare_annotation.py*)
 
La formule permet de s'interesser aux différences d'utilisation des exons

Chaque ligne de la classe DEXSeqDataSet contient dans chaque colonne des données de comptage d'un exon donné  (« *this* ») et les données de comptage de la somme des autres exons appartenant au même gène (« *others* »).

Toutes ces informations sont rappelées dans la colonne `colData()`.


```
dxd = DEXSeqDataSetFromHTSeq(
 countFiles,
 sampleData=sampleTable,
 design= ~ sample + exon + condition:exon,
  flattenedfile=flattenedFile )
```

- Pour accéder aux 5 premières lignes des données de décompte : `head( counts(dxd), 5 )`

- Pour accéder aux 5 premières lignes des données de décompte appartenant aux régions exoniques (« this ») (sans afficher la somme des décomptes du reste des exons du même gène) : `head( featureCounts(dxd), 5 )`


#### **Normalisation des données** ####


`dxd = estimateSizeFactors( dxd )`


#### **Estimation de la variance** ####
L'estimation de la variance va permettre distinguer les variations techniques et biologiques (bruit) des effets réels d'evenement d'épissage.


Calcul de la variance par exon : 

   - Estimation de la vraisemblance de profil ajustée de Cox-Reid
   
   - Relation variance-moyenne est ajustée aux valeurs de la variance individuelle
   
   - Valeurs ajustées sont prises afin de ramener les estimations par exon vers les valeurs ajustées



##### **Test de l'analyse différentielle des exons** #####

La fonction va effectuer des tests pour chaque exon de chaque gène : elle va ajuster un modèle linéaire généralisé et comparer au modèle plus petit (modèle nul) avec la formule : 


Modèle complet : `formulaFullModel    =  ~ sample + exon +  condition:exon`

Modèle réduit : `formulaReducedModel =  ~ sample + exon `


`dxd = estimateDispersions( dxd, formula = formulaFullModel )`


Visualisation de la variance par exon par rapport au nombre moyen normalisé :
`plotDispEsts(dxd)`

La fonction `testForDEU` effectue des tests pour chaque exon de chaque gène et va prendre en compte les deux formules.

```
dxd = testForDEU( dxd, 
                  reducedModel = formulaReducedModel, 
                  fullModel = formulaFullModel )
```

Pour avoir un tableau récapitulatif des résultats : 

`dxr2 = DEXSeqResults( dxd )`


La description de chacune des colonnes de l'objet DEXSeqResults se trouve dans les colonnes de métadonnées : `mcols(dxr2)$description`

- Pour savoir combien de régions exoniques sont significatives avec un taux de fausses découvertes de 10 % : `table ( dxr1$padj < 0.1 )`

- Pour voir combien de gènes subissent un événements d'épissage : `table( tapply( dxr1$padj < 0.1, dxr1$groupID, any ) )`

- Calcule du nombre de gènes, avec un FDR de 10% avec au moins un exon épissé différentiellement : `somme ( perGeneQValue (dxr) <  0,1 )`


### **Visualisation** ###

DEXSeq fournit un moyen de visualiser les résultats d'une analyse avec `plotDEXSeq()`

`plotDEXSeq( dxr2, "FBgn0010909", legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )`

Ce graphique montre les valeurs d’expression ajustées de chacun des exons du gène FBgn0010909, pour chacune des deux conditions, traitées et non traitées.

Il est aussi possible de visualiser les modèles de transcription afin d'analyser les résultats différentiels des épissages d'exon dans les isoformes : 
`plotDEXSeq( dxr2, "FBgn0010909", displayTranscripts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )`


Il existe aussi une visualisation qui permet d'examiner les valeurs de comptages des échantillons individuels.

`plotDEXSeq( dxr2, "FBgn0010909", expression=FALSE, norCounts=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )` 

Pour observer les changements d'épissage des exons sans a priori sur la régulation du gène (à la hausse ou à la baisse) : 

`plotDEXSeq( dxr2, "FBgn0010909", expression=FALSE, splicing=TRUE, legend=TRUE, cex.axis=1.2, cex=1.3, lwd=2 )`


Afin d'avoir une vue d’ensemble détaillée et facilement consultable de tous les résultats d’analyse, le package fournit un générateur de rapports HTML.

Elle permet d'avoir une table détaillée de résultats avec des liens vers des tracés pour les résultats significatifs.

`DEXSeqHTML( dxr2, FDR=0.05, color=c("#FF000080", "#0000FF80") )`



























