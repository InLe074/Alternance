##### LEBIB Inès - IDRS #####
###### 2023-11-16 ######
#### **Introduction** ####

L’épissage alternatif est l’épissage de différents exons à partir d'un même ARN pré-messager. Il permet de diversifier les séquences protéiques ainsi que leurs fonctions et leurs régulations. 

rMATS -turbo est un outil développé en Python qui détecte 5 modèles d’épissage : saut d’exon (SE), sites d’épissage alternatifs en 5’ (A5SS), sites d’épissage alternatifs en 3’ (A3SS), exons mutuellement exclusifs (MXE) ou intron retenu (rétention d’intron) (RI). Le nombre de reads de support peut être compté par les reads de jonction uniquement (JC) ou par les reads de jonction et d’exon (JCEC). Le fichier de sortie des différentes méthodes de comptage est indiqué dans le nom du fichier.

Liens utiles pour la documentation de rMATS : 

* https://github.com/Xinglab/rmats-turbo 
* https://rnaseq-mats.sourceforge.io/ 

Toute l'anayse se fera en **ligne de commande** sur **Python** dans un environnement **conda**. 

Pour savoir comment créer un environnement conda sur JupyterLab, je vous invite à consulter le readme **Install_Python_Packages** ou le pdf **RDP-Install Python packages on Jupyter Notebook & JupyterLab Posit Workbench-220523-121139**

#### **Dépendances requises** ####

* Python (3.6.12 or 2.7.15)
* Cython (0.29.21 or 0.29.15 for Python 2)
* BLAS, LAPACK
* GNU Scientific Library (GSL 2.5)
* GCC (>=5.4.0)
* gfortran (Fortran 77)
* CMake (3.15.4)
* PAIRADISE (optional)
* Samtools (optional)
* STAR (optional)

#### **Etape 1 : Cloner le github** ####

La première étape est de cloner le github [rMATS turbo](https://github.com/Xinglab/rmats-turbo#rmats-turbo-v420) : `git clone https://github.com/Xinglab/rmats-turbo.git` en **ligne de commande** sur **Python** dans un environnement **conda**. 


#### **Etape 2 : Construire rMATS** ####

La seconde étape est de construire rMATS avec la ligne de commande (en se mettant dans le même répertoire où a été cloné le github) : 
`./build_rmats`
Puis lancez : 
`python rmats.py {arguments}`

L'utilisation de build_rmats se fait avec un environnement conda :

```
./build_rmats --conda

--conda: create a conda environment for Python and R dependencies
```

**Attention** : Durant toutes les prochaines étapes, lors de l'écriture du chemin d'accès, toujours bien commencer par la **racine du chemin** !

#### **Etape 3 : L'analyse** ####

L'analyse de rMATS comporte deux étapes : la préparation de l'analyse et le résultat.

##### **Etape 3.1 : La préparation** #####

Lors de l'étape de préparation, les fichiers d'entrée sont traités et un résumé est enregistré dans les fichiers `.rmats`. 

Les fichiers `.rmats` contiennent des informations sur les événements d'épissage repérés dans les données analysées tels que les types d'événements d'épissage, les statistiques, les scores associés et des annotations génomiques indiquant la localisation de ces événements dans le génome. Leurs chemins complets sont spécifiés dans le fichier `.txt` donné en entrée.

Supposons que nous ayons 8 BAM et deux machines dotées chacune de 4 threads CPU. Chaque machine peut exécuter l'étape de préparation sur 4 BAM simultanément.

Lors de l'écriture du chemin d'accès, toujours bien commencer par la **racine du chemin** !

Divisez les BAM en deux groupes, suivant leurs conditions. 

* `/path/to/Condition_1.txt` :

`/path/to/1.bam,/path/to/2.bam,/path/to/3.bam,/path/to/4.bam`

* `/path/to/Condition_2.txt` :

`/path/to/5.bam,/path/to/6.bam,/path/to/7.bam,/path/to/8.bam`


L'étape de la préparation s'exécute avec `Condition_1.txt` : 
`python rmats.py --b1 /path/to/Condition_1.txt --gtf /path/to/the.gtf -t paired --readLength 150 --nthread 4 --od /path/to/output --tmp /path/to/tmp_output_Condition_1 --task prep`

Nous répétons la commande en l'adaptant au fichier `Condition_2.txt` :
`python rmats.py --b1 /path/to/Condition_2.txt --gtf /path/to/the.gtf -t paired --readLength 150 --nthread 4 --od /path/to/output --tmp /path/to/tmp_output_Condition_2 --task prep`


##### **Etape 3.2 : L'analyse** #####

Lors de l'étape de l'analyse, les fichiers `.rmats` sont lus et les fichiers de sortie finaux sont créés. Ces fichiers contiendront tous les événements d'épissages trouvés par l'outil rMATS. 

L'argument `--task` permet à l'étape de préparation de rMATS d'être exécuté indépendamment pour différents sous-ensembles de fichiers BAM d'entrée.

L'étape de l'analyse peut être exécutée sur les fichiers `.rmats` générés indépendamment. Cela permet d'exécuter le calcul à des moments différents et/ou sur des machines différentes.

Divisez les BAM en deux groupes. Cette répartition sert à comparer statistiquement les deux groupes.

Lors de l'écriture du chemin d'accès, toujours bien commencer par la **racine du chemin** !

* `/path/to/Condition_1.txt` :

`/path/to/1.bam,/path/to/2.bam,/path/to/3.bam,/path/to/4.bam`

* `/path/to/Condition_2.txt` :

`/path/to/5.bam,/path/to/6.bam,/path/to/7.bam,/path/to/8.bam`


Copiez tous les fichiers .rmats dans un seul répertoire avec le script `cp_with_prefix.py`



```
python cp_with_prefix.py prep_1_ /path/to/output_Condition_1/ /path/to/tmp_output_Condition_1/*.rmats
python cp_with_prefix.py prep_2_ /path/to/output_Condition_2/ /path/to/tmp_output_Condition_2/*.rmats
```

Puis exécuter :
`python rmats.py --b1 /path/to/Condition_1.txt --b2 /path/to/Condition_2.txt --gtf /path/to/the.gtf -t paired --readLength 50 --nthread 4 --od /path/to/output --tmp /path/to/tmp_output_post --task post`

Le répertoire de sortie contiendra les fichiers nécessaires `fromGTF.[AS].txt` et `{JC,JCEC}.raw.input.[AS].txt`. Les fichiers `fromGTF.[AS].txt` peuvent être utilisés tels quels pour toutes les comparaisons impliquant les échantillons, mais les informations pertinentes à une comparaison spécifique sont extraites des fichiers `{JC,JCEC}.raw.input.[AS].txt`.


#### **Output - Fichiers de Sortie** ####

Dans rMATS-turbo, chaque modèle d'épissage alternatif a un ensemble correspondant de fichiers de sortie. Dans les modèles de noms de fichiers ci-dessous, `[AS_Event]` est remplacé par l'un des cinq modèles d'épissage alternatifs de base : exon sauté (`SE`), sites d'épissage alternatifs en 5' (`A5SS`), sites d'épissage alternatifs en 3' (`A3SS`), exons mutuellement exclusifs (`MXE`) ou intron retenu (`RI`). Comme le montre le diagramme, le nombre de lectures de support peut être compté par les lectures de jonction uniquement (`JC`) ou par les lectures de jonction et d'exon (`JCEC`).

#### Répertoire `--od`

`--od` contient les fichiers de sortie finaux de l'étape de post-production :

- `[AS_Event].MATS.JC.txt` : Sortie finale qui contient la liste des événements et les comptes de read. Seuls les reads de jonction d'épissage sont comptés.

- `[AS_Event].MATS.JCEC.txt` : Sortie finale contenant la liste des événements et le comptes de reads. Les reads de jonctions d’épissage et les reads exoniques sont comptés.

- `fromGTF.[AS_Event].txt`: Tous les événements d’épissage alternatif identifiés dérivés du GTF et de l’ARN.

- `fromGTF.novelJunction.[AS_Event].txt` : Événements d'épissage alternatif qui ont été identifiés uniquement après avoir pris en compte l'ARN (par opposition à l'analyse du GTF seul). N'inclut pas les événements avec un site d'épissage non annoté.

- `fromGTF.novelSpliceSite.[AS_Event].txt` : Ce fichier contient uniquement les événements qui incluent un site d'épissage non annoté. Seulement pertinent si l'option `--novelSS` est activée.

- `JC.raw.input.[AS_Event].txt` : Comptes d'événements incluant uniquement les reads qui couvrent les jonctions définies par rMATS.

- `JCEC.raw.input.[AS_Event].txt` : Comptes d'événements incluant à la fois les reads qui couvrent les jonctions définies par rMATS et les reads qui ne traversent pas une frontière d'exon.

#### Colonnes partagées

- **ID**: ID d'événement rMATS
- **GeneID**: Identifiant du gène
- **geneSymbol**: Nom du gène
- **chr**: Chromosome
- **strand**: Brin du gène
- **IJC_SAMPLE_1**: Nombre d'inclusions pour l'échantillon 1. Les réplicats sont séparés par des virgules.
- **SJC_SAMPLE_1**: Comptes à sauter pour l'échantillon 1. Les réplicats sont séparés par des virgules.
- **IJC_SAMPLE_2**: Les inclusions comptent pour l'échantillon 2. Les réplicats sont séparés par des virgules.
- **SJC_SAMPLE_2**: Comptes à sauter pour l'échantillon 2. Les répétitions sont séparées par des virgules.
- **IncFormLen**: Formulaire de longueur d'inclusion, utilisé pour la normalisation.
- **SkipFormLen**: Longueur du formulaire à sauter, utilisée pour la normalisation.
- **PValue**: Importance de la différence d’épissage entre les deux groupes d’échantillons. (Uniquement disponible si le modèle statistique est activé)
- **FDR**: Taux de fausses découvertes calculé à partir de la valeur p. (Uniquement disponible si le modèle statistique est activé)
- **IncLevel1**: Niveau d'inclusion pour l'échantillon 1. Les réplicats sont séparés par des virgules. Calculé à partir de comptes normalisés.
- **IncLevel2**: Niveau d'inclusion pour l'échantillon 2. Les réplicats sont séparés par des virgules. Calculé à partir de comptes normalisés.
- **IncLevelDifference**: moyenne(IncLevel1) - moyenne(IncLevel2).

###  Colonnes spécifiques à l'événement (coordonnées de l'événement)

#####  **SE**
- `exonStart_0base` `exonEnd` `upstreamES` `upstreamEE` `downstreamES` `downstreamEE`
- Le formulaire d'inclusion comprend l'exon cible (`exonStart_0base`, `exonEnd`).

#####  **MXE**
- `1stExonStart_0base` `1stExonEnd` `2ndExonStart_0base` `2ndExonEnd` `upstreamES` `upstreamEE` `downstreamES` `downstreamEE`
- Si le brin est `+`, alors la forme d'inclusion inclut le 1er exon (`1stExonStart_0base`, `1stExonEnd`) et saute le 2ème exon.
- Si le brin est `-`, alors la forme d'inclusion inclut le 2ème exon (`2ndExonStart_0base`, `2ndExonEnd`) et saute le 1er exon.

##### **A3SS, A5SS**
- `longExonStart_0base` `longExonEnd` `shortES` `shortEE` `flankingES` `flankingEE`
- Le formulaire d'inclusion inclut l'exon long (`longExonStart_0base`, `longExonEnd`) au lieu de l'exon court (`shortES` `shortEE`).

##### **IR**
- `riExonStart_0base` `riExonEnd` `upstreamES` `upstreamEE` `downstreamES` `downstreamEE`
- La forme d'inclusion inclut (conserve) l'intron (`upstreamEE`, `downstreamES`).



#####  **summary.txt**
- Bref résumé de tous les types d’événements d’épissage alternatifs. Comprend le nombre total d'événements et le nombre d'événements significatifs. Par défaut, les événements sont considérés comme significatifs si FDR <= 0,05. Le résumé peut être régénéré avec différents critères en exécutant rMATS_P/summary.py.

#####  **--tmp**
- Contient les fichiers intermédiaires générés par l'étape de préparation :
  - `[datetime]_[id].rmats` : Résumé généré à partir du traitement d'un BAM.
  - `[datetime]_bam[sample_num]_[replicate_num]/Aligned.sortedByCoord.out.bam` : Résultat du mappage des fichiers FASTQ d'entrée
  - `[datetime]_read_outcomes_by_bam.txt` : nombre de reads utilisés à partir de chaque BAM ainsi que le nombre de raisons pour lesquelles les reads n'ont pas pu être utilisés



#### Tous les arguments ####

```
python rmats.py -h

usage: rmats.py [options]

optional arguments:
  -h, --help            show this help message and exit
  --version             show program's version number and exit
  --gtf GTF             An annotation of genes and transcripts in GTF format
  --b1 B1               A text file containing a comma separated list of the
                        BAM files for sample_1. (Only if using BAM)
  --b2 B2               A text file containing a comma separated list of the
                        BAM files for sample_2. (Only if using BAM)
  --s1 S1               A text file containing a comma separated list of the
                        FASTQ files for sample_1. If using paired reads the
                        format is ":" to separate pairs and "," to separate
                        replicates. (Only if using fastq)
  --s2 S2               A text file containing a comma separated list of the
                        FASTQ files for sample_2. If using paired reads the
                        format is ":" to separate pairs and "," to separate
                        replicates. (Only if using fastq)
  --od OD               The directory for final output from the post step
  --tmp TMP             The directory for intermediate output such as ".rmats"
                        files from the prep step
  -t {paired,single}    Type of read used in the analysis: either "paired" for
                        paired-end data or "single" for single-end data.
                        Default: paired
  --libType {fr-unstranded,fr-firststrand,fr-secondstrand}
                        Library type. Use fr-firststrand or fr-secondstrand
                        for strand-specific data. Only relevant to the prep
                        step, not the post step. Default: fr-unstranded
  --readLength READLENGTH
                        The length of each read
  --variable-read-length
                        Allow reads with lengths that differ from --readLength
                        to be processed. --readLength will still be used to
                        determine IncFormLen and SkipFormLen
  --anchorLength ANCHORLENGTH
                        The "anchor length" or "overhang length" used when
                        counting the number of reads spanning splice
                        junctions. A minimum number of "anchor length"
                        nucleotides must be mapped to each end of a given
                        junction. The minimum value is 1 and the default value
                        is set to 1 to make use of all possible splice
                        junction reads.
  --tophatAnchor TOPHATANCHOR
                        The "anchor length" or "overhang length" used in the
                        aligner. At least "anchor length" NT must be mapped to
                        each end of a given junction. The default is 1. (Only
                        if using fastq)
  --bi BINDEX           The directory name of the STAR binary indices (name of
                        the directory that contains the SA file). (Only if
                        using fastq)
  --nthread NTHREAD     The number of threads. The optimal number of threads
                        should be equal to the number of CPU cores. Default: 1
  --tstat TSTAT         The number of threads for the statistical model. If
                        not set then the value of --nthread is used
  --cstat CSTAT         The cutoff splicing difference. The cutoff used in the
                        null hypothesis test for differential splicing. The
                        default is 0.0001 for 0.01% difference. Valid: 0 <=
                        cutoff < 1. Does not apply to the paired stats model
  --task {prep,post,both,inte,stat}
                        Specify which step(s) of rMATS to run. Default: both.
                        prep: preprocess BAMs and generate a .rmats file.
                        post: load .rmats file(s) into memory, detect and
                        count alternative splicing events, and calculate P
                        value (if not --statoff). both: prep + post. inte
                        (integrity): check that the BAM filenames recorded by
                        the prep task(s) match the BAM filenames for the
                        current command line. stat: run statistical test on
                        existing output files
  --statoff             Skip the statistical analysis
  --paired-stats        Use the paired stats model
  --novelSS             Enable detection of novel splice sites (unannotated
                        splice sites). Default is no detection of novel splice
                        sites
  --mil MIL             Minimum Intron Length. Only impacts --novelSS
                        behavior. Default: 50
  --mel MEL             Maximum Exon Length. Only impacts --novelSS behavior.
                        Default: 500
  --allow-clipping      Allow alignments with soft or hard clipping to be used
  --fixed-event-set FIXED_EVENT_SET
                        A directory containing fromGTF.[AS].txt files to be
                        used instead of detecting a new set of events

```

