##### LEBIB Inès - IDRS #####

### **Introduction** ###

L’épissage alternatif est l’épissage de différents exons à partir de l’ARN pré-messager. Il permet de diversifier les séquences des protéines ainsi que leurs fonctions et leurs régulations. 



SpliceWiz est une interface graphique pour l'épissage différentiel et la visualisation dans R.


SpliceWiz est offre une gestion rapide des fichiers BAM d'alignement et SpliceWiz offre la possibilité de faire des traitements en multithread


### **Packages et librairies à installer** ###

```

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SpliceWiz")

library(SpliceWiz)

```

### **Construire la référence de SpliceWiz**

#### La référence SpliceWiz est utilisée pour quantifier l'épissage alternatif dans les fichiers BAM, ainsi que dans le classement en aval, l'analyse différentielle et la visualisation.

```
ref_path <- file.path(tempdir(), "Reference")
buildRef(
  reference_path = ref_path,
  fasta = "/home/path/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
  gtf = "/home/path/GRCh38.99.Chr.gtf",
  ontologySpecies = "Homo sapiens"
)
```


### Traiter les fichiers BAM à l'aide de SpliceWiz

```

pb_path <- file.path(tempdir(), "pb_output")

bam_path <- "/home/Directory/to/BAM"

bam_files <- list(bam_files)

Sample_names <- c("")

```

Traitez ces fichiers BAM à l'aide de SpliceWiz

```
processBAM(
  bamfiles = bam_files,
  sample_names = Sample_names,
  reference_path = ref_path,
  output_path = pb_path
)


```


### Rassembler les données 

```
expr <- findSpliceWizOutput(pb_path)


nxtse_path <- file.path(tempdir(), "NxtSE_output")
collateData(
  Experiment = expr,
  reference_path = ref_path,
  output_path = nxtse_path
)

```


### Importer les données 

```
se <- makeSE(nxtse_path)

```

### Analyse différentielle 

### Attribution des annotation aux échantillons

```
colData(se)$condition <- rep(c("DMSO", "E7107"), each = 3)
colData(se)$batch <- rep(c("Condition1", "Condition2"), each = 3)
```

### Filtrage des événements de haute confiance

SpliceWiz propose des filtres par défaut pour identifier et supprimer les événements d'épissage alternatifs à faible confiance.

```
se.filtered <- se[applyFilters(se),]
```

### Effectuer une analyse différentielle












