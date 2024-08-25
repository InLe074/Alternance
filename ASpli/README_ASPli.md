##### LEBIB Inès - IDRS #####
###### 2024-07-02 ######
#### **Introduction** ####

L’épissage alternatif est l’épissage de différents exons à partir d'un même ARN pré-messager. Il permet de diversifier les séquences protéiques ainsi que leurs fonctions et leurs régulations. 

ASpli, est un outil bioinformatique implémenté en R, qui permet l'identification de changements dans les événements d'épissage alternatif annotés et nouveaux et peut traiter des plans expérimentaux simples, multifactoriels ou appariés. 

Liens utiles pour la documentation de ASpli : 

* https://www.bioconductor.org/packages/devel/bioc/vignettes/ASpli/inst/doc/ASpli.pdf
* https://doi.org/10.1093/bioinformatics/btab141


Toute l'anayse se fera en langage de programmation **R**

### **Packages et librairies à installer** ###

```
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ASpli")
BiocManager::install("GenomeInfoDb")
BiocManager::install("GenomicRanges")
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
BiocManager::install("GenomicFeatures")


library(ASpli)
library(GenomicFeatures)

```


### **Etape 1 : Pré-processing du fichier GTF** ###

Cette étape va permettre d'avoir un fichier au format TxDb (Transcript DataBase) à partir d'un fichier GTF. Il va être utilisé pour l'étude de l'expression génique et l'annotation génomique. 

```
gtfFileName <- '/path/to/GRCh38.99.Chr.gtf'
genomeTxDb <- makeTxDbFromGFF( gtfFileName )
saveDb(genomeTxDb,file="gene.sqlite")
```

### **Etape 2 : Features extraction** ###

Cette ligne de commande va permettre d'extraire les caractéristiques importantes des données génomiques qui va être utilisé pour le comptage des reads par exemple.

```
features <- binGenome( genomeTxDb )
```

### **Etape 3 : Préparation des fichiers BAM et définition des cibles** ###

Création de la liste des fichiers BAM

```
BAMFiles <- c('/path/to/file_1.bam',
              '/path/to/file_2.bam', ...)
              
```

Création d'un data.frame pour les cibles avec les données phénotypiques

```
row_names <- c( "Control_1", "Control_2", "Control_3", 
                "Treatment_1", "Treatment_2", "Treatment_3")

targets <- data.frame(row.names = row_names,
                      bam = BAMFiles[1:6],
                      f1 = c( "Control", "Control", "Control", 
                              "Treatment", "Treatment", "Treatment"),
                      stringsAsFactors = FALSE)


mBAMs <- data.frame( bam = targets$bam[c(1:6)],
                     condition = c( "Control", "Control", "Control", 
                              "Treatment", "Treatment", "Treatment"))
                              
```

Comptage des reads alignés avec la fonction `gbCounts`.

`minReadLength=90` spécifie la longueur de lecture minimale pour qu'une lecture soit considérée dans l'analyse.


```
gbcounts <- gbCounts(features=features, targets=targets,
                     minReadLength = 90, maxISize = 50000)

```

Comptage de novo des jonctions et estimation les signaux d'épissage avec la fonction `jCounts`.


```                                                     
asd <- jCounts(counts=gbcounts, features=features, minReadLength=90)

```

Analyser des événements différentiels des bins de gènes avec la fonction `gbDUreport`.

`contrast = c(1, -1)` définit le contraste pour les comparaisons entre les conditions.

```
gb <- gbDUreport(gbcounts, contrast = c(1, -1))
```

Rapport des événements d'épissage différentiel des jonctions avec la fonction `jDUreport`.

Prend en entrée `asd` de l'étape précédente et applique le même contraste.

```
jdur <- jDUreport(asd, contrast=c(1, -1))

```

Rapports des événements différentielse des bins et des jonctions pour créer un rapport d'épissage global avec la fonction `splicingReport`.

```
sr <- splicingReport(gb, jdur, counts=gbcounts)
```

Intégrons des signaux d'épissage le long des régions génomiques avec la fonction `integrateSignals`.

```
is <- integrateSignals(sr,asd)

```

Enfin, nous exportons les signaux d'épissage intégrés et les rapports associés vers un dossier de résultats.

`exportIntegratedSignals` exporte les données analysées vers un dossier spécifique.

`output.dir` = "ASPli_Resultat"' définit le dossier de sortie pour les fichiers de résultats.

```
exportIntegratedSignals(is,sr=sr,
                        output.dir = "ASPli_Resultat",
                        counts=gbcounts,features=features,asd=asd,
                        mergedBams = mBAMs)
                        
```


#### **Détails des outputs ** ####


#### Structure de l'Objet ASpliCounts

L'objet `ASpliCounts` inclut les données de comptage des gènes, les densités de reads des gènes, les comptes des bins et les densités de reads des bins, ainsi que les comptes des jonctions.

#### Comptes de Gènes (countsg) et Densités de Reads des Gènes (rdsg)

- `row.names`: Nom du gène tel que rapporté dans les données d'annotation.
- `symbol`: Un nom optionnel pour le gène, qui doit être fourni au moment de l'extraction des caractéristiques.
- `locus_overlap`: Loci chevauchants.
- `gene_coordinates`: Format chromosome : start-end.
- `start`, `end`, `length`: Position de début, de fin et longueur du gène.
- `effective_length`: Longueur efficace du gène en utilisant seulement les exons annotés.
- `sample_data`: Données de comptage/densité de reads pour les gènes (une colonne par échantillon).

#### Comptes de Bins (countsb) et Densités de Reads des Bins (rdsb)

- `row.names`: binGenome: Binning du génome.
- `feature`: E pour les bins exons, I pour les bins introns et Io pour les introns avant division.
- `event`: Événement d'épissage alternatif assigné.
- `sample_data`: Données de comptage/densité de reads pour les bins (une colonne par échantillon).

#### Comptes de Jonctions (countsj)

- `row.names`: Nom de la jonction au format chromosome.start.end.
- `junction`: Si la jonction coïncide avec une jonction déduite de l'annotation, le nom est affiché tel qu'il est donné dans la section 6.1:binGenome: Binning du génome, sinon c'est noHit.
- `gene`: Locus qui contient la jonction.
- `strand`: Brin du gène.
- `multiple_hit`: Oui si la jonction couvre plusieurs gènes.
- `bin_spanned`: Noms des bins traversés par la jonction.
- `j_within_bin`: Si la jonction se trouve dans un seul bin, le nom de ce bin est affiché.
- `sample_data`: Données de comptage de reads pour les jonctions (une colonne par échantillon).



#### Données de Jonction pour l'Analyse de l'Épissage Alternatif en RNA-seq

Les objets `countsei1` et `countsei2` semblent contenir des données relatives aux jonctions dans le contexte de l'épissage alternatif analysé par RNA-seq.

#### Structure de l'Objet de Jonction

Les données de jonction contiennent plusieurs composantes clés :

- `row.names` : Le nom de la jonction au format `chromosome.start.end`, `binGenome: Binning the genome`.
  
- `junction` : Indique si la jonction coïncide avec une jonction déduite de l'annotation. Si c'est le cas, le nom est affiché tel que donné : `binGenome: Binning the genome`, sinon il contient `noHit`.
  
- `gene` : Le nom du locus qui contient la jonction.

- `strand` : Le brin sens du gène.

- `multipleHit` : Indique `yes` si la jonction couvre plusieurs gènes.

- `sample data` : Les comptes de jonction (une colonne par échantillon).


#### iPIR (Intron Percent Inclusion Ratio)

Cet indice est calculé pour estimer la fréquence d'inclusion des introns.

- `event`: Type d'événement assigné par ASpli lors du binning.
- `J1`: Liste de jonctions séparées par des points-virgules dont l'extrémité correspond au début de l'intron.
- `J2`: Liste de jonctions séparées par des points-virgules dont l'extrémité correspond à la fin de l'intron.
- `J3`: Liste de toutes les jonctions chevauchant l'intron.

La formule pour calculer le iPIR pour chaque condition est la suivante :

```math
PIR = (J1 + J2) / (J1 + J2 + 2 * J3)
```

#### altPSI (Alternative Percent Spliced In)

Cet indice est calculé pour estimer la fréquence d'épissage alternatif.

- `event`: Type d'événement assigné par ASpli lors du binning.
- `J1(J2)`: Liste de jonctions séparées par des points-virgules avec une extrémité correspondant à l'extrémité alternative 5' ou 3'.
- `J3`: Liste de toutes les jonctions chevauchant l'intron.

La formule pour calculer le altPSI pour chaque condition est la suivante :

```math
PSI = (J1(J2)) / (J1(J2) + J3)
```

où J12 est utilisé comme J1 s'il s'agit d'un événement d'épissage alternatif 5', et comme J2 s'il s'agit d'un événement d'épissage alternatif 3'.


#### esPSI (Exon Skipping Percent Spliced In)

- `event`: Type d'événement assigné par ASpli lors du binning.
- `J1`: Liste de jonctions séparées par des points-virgules avec une extrémité sur l'exon alternatif.
- `J2`: Liste de jonctions séparées par des points-virgules avec une extrémité sur l'exon alternatif.
- `J3`: Liste de jonctions d'exclusion de l'exon alternatif.

Les colonnes entre J1-J2 et J2-J3 représentent les comptes de jonctions dans les échantillons pour chaque bin. Le PSI est calculé pour chaque condition comme suit :

```math
PSI = (J1 + J2) / (J1 + J2 + 2 * J3)
```


#### junctionsPIR (Percent Intron Retention)

- `PIR`: Pour chaque jonction expérimentale, utilisant les comptes `ei1` et `ie2`. La jonction d'exclusion est la jonction elle-même.
- `hitIntron` et `hitIntronEvent`: Si la jonction correspond à un bin, nom et type d'événement assignés par ASpli.

Les colonnes entre représentent les comptes J1, J2, et J3 dans les différents échantillons. Le PIR est calculé pour chaque condition comme suit :

```math
PIR = (J1 + J2) / (J1 + J2 + 2 * J3)
```



#### junctionsPJU (Percent Junction Usage)

Cette section décrit les métriques fournies pour évaluer les événements  des jonctions.

- `Rowname`: Plage de J3.
- `Junction`: Nom de la jonction.
- `Gene`: Gène auquel elle appartient.
- `Strand`: Brin du gène.
- `multipleHit`: Indique si la jonction chevauche plusieurs gènes.
- `symbol`: Symbole du gène.
- `gene_coordinates`: Coordonnées du gène.
- `bin_spanned`: Liste de bins séparés par des points-virgules couverts par cette jonction.
- `j_within_bin`: Autres jonctions dans cette région.
- `StartHit`: Jonctions partageant le début avec cette jonction et le calcul de PJU/1 est `PJU/1 = J3/(J1 + J3)` pour chaque condition.
- `EndHit`: Jonctions partageant la fin avec cette jonction et le calcul de PJU/2 est `PJU/2 = J3/(J2 + J3)` pour chaque condition.

Les colonnes entre contiennent :

- Les comptes de J3 dans les différents échantillons pour chaque région.

ASpli fournit un cadre intégré pour l'analyse des jonctions dans les données de séquençage d'ARN, permettant d'identifier des événements d'épissage alternatif.

#### Calculs supplémentaires

- `PJU/1` et `PJU/2` sont calculés en utilisant les formules susmentionnées pour chaque condition.
- Les plages pour `StartHit` sont J1 et pour `EndHit` sont J2.



#### ASpliJDU

Cet objet fournit des informations sur les jonctions connues ou nouvelles et les bins couverts.

- `spanned bins`: Les bins que la jonction couvre entièrement.
- `exontron`: Indique si la jonction est entièrement incluse dans un bin, ce qui peut signifier que cet événement d'épissage alternatif est un exontron possible.

#### localec

Informations sur les changements statistiquement significatifs dans l'utilisation des jonctions à l'intérieur des clusters de `locale-junction`.

- `size`: Nombre de jonctions appartenant au cluster.
- `cluster.LR`: Rapport de vraisemblance différentiel d'utilisation du cluster.
- `pvalue, FDR`: Valeur p et FDR différentiels d'utilisation du cluster.
- `range`: Coordonnées génomiques du cluster.
- `participation`: Valeur maximale de participation de la jonction à l'intérieur du cluster.
- `dParticipation`: Participation delta de la jonction significative à l'intérieur du cluster.

#### localej

Informations relatives à une jonction spécifique au sein d'un cluster.

- `cluster`: Nom du cluster auquel la jonction appartient.
- `log.mean`: Logarithme de la moyenne des comptes pour cette jonction sur toutes les conditions.
- `logFC`: Logarithme du changement de pli de la jonction à travers les conditions.
- `pvalue and FDR`: Valeur p et FDR de la jonction.
- `annotated`: Indique si la jonction est annotée ou nouvelle.
- `participation`: Valeur de participation maximale observée à travers des conditions contrastées.
- `dParticipation`: Participation delta de la valeur de participation maximale observée à travers des conditions contrastées.

Les colonnes entre indiquent les comptes de jonctions pour tous les échantillons.

#### anchorc

- `cluster.LR`: Rapport de vraisemblance des événements différentiels du cluster.
- `pvalue and FDR`: Valeur p et FDR des événements différentiels du cluster.


#### anchorj

- `log.mean`: Logarithme de la moyenne des comptages pour cette jonction à travers toutes les conditions.
- `logFC`: Changement de pli logarithmique de la jonction à travers les conditions.
- `LR`: Rapport de vraisemblance des événements différentiels de la jonction.
- `pvalue and FDR`: p-value de la jonction et FDR correspondant.
- `J1.pvalue :`: P-value de la jonction J1.
- `J2.pvalue :`: P-value de la jonction J2.
- `NonUniformity`: Si un test de non uniformité a été effectué, les nombres plus proches de zéro signifient l'uniformité et ceux plus proches de un signifient la non uniformité.
- `dPIR`: Delta PIR de la jonction.
- `annotated :`:  Si la jonction est annotée ou nouvelle.
- `J counts:`: Comptages des jonctions pour tous les échantillons.



#### jir

- `J3`: jonctions J3.
- `logFC`: Changement de pli logarithmique de la jonction à travers les conditions.
- `LR`: Rapport de vraisemblance des événements différentiels de la jonction.
- `pvalue and FDR`: P-value de la jonction et FDR correspondant.
- `NonUniformity`: Si un test de non uniformité a été effectué, les nombres plus proches de zéro signifient l'uniformité et ceux plus proches de un signifient la non uniformité.
- `dPIR`: Delta PIR de la jonction.
- `multiplicity :`:  Si plusieurs jonctions traversent la région.
- `J counts:`: Comptages des jonctions pour tous les échantillons.



#### jes

- `event`: Type d'événement attribué par ASpli lors du regroupement.
- `J3`: jonctions J3.
- `logFC`: Changement de pli logarithmique de la jonction à travers les conditions.
- `LR`: Rapport de vraisemblance des événements différentiels de la jonction.
- `pvalue and FDR`: P-value de la jonction et FDR correspondant.
- `NonUniformity`: Si un test de non uniformité a été effectué, les nombres plus proches de zéro signifient l'uniformité et ceux plus proches de un signifient la non uniformité.
- `dPSI`: Delta PSI de la jonction.
- `multiplicity :`:  Si plusieurs jonctions traversent la région.
- `J counts:`: Comptages des jonctions pour tous les échantillons.



##### jalt

- `event`: Type d'événement assigné par ASpli lors du binning.
- `J3`: Jonctions.
- `logFC`: Changement de pli logarithmique de la jonction à travers les conditions.
- `log.mean`: Logarithme de la moyenne des comptes à travers toutes les conditions pour cette jonction.
- `logFC`: Changement de pli logarithmique de la jonction à travers les conditions.
- `LR`: Ratio de vraisemblance des événements différentiels des jonctions.
- `pvalue, FDR`: P-valeur de la jonction et FDR correspondant.
- `dPSI`: Delta PSI de la jonction.
- `multiplicity`: Indique si de multiples jonctions traversent la région.
- `J counts`: Comptes des jonctions pour tous les échantillons.

















