# A Pipeline for Diversity Outbred genome reconstruction using GigaMUGA microarray



[TOC]

---

## 1. Background

### 1.1 About GigaMUGA

GeneSeek Genomic Profiler (GCP) mouse genotyping arrays enable genome wide profiling of Single Nucleotide Polymorphisms (SNPs). [GigaMUGA](https://www.illumina.com/products/by-type/microarray-kits/ggp-mouse.html) (GM), first introduced in 2012, is a third-generation array that extends the capabilities of previous GCP arrays (e.g. MegaMUGA) by adding more than 65,000 new SNPs. GM can profile over More than 143,000 SNPs from 159 inbred mouse lines with an average marker spacing of ~22.5 kb.

### 1.2 Importance of DO genome reconstruction in disease genetics research

Much like humans, Diversity Outbred (DO) mice are a population of high genetic diversity. Do mice are reproduced from an initial cross of eight inbred mouse strains (commonly known as founder strains), then mixed and bred for many generations to generate a high level of genetic diversity. This makes DO mice an extremely valuable tool for studying complex genetic traits and diseases. In disease genetics research, the goal is often to identify genetic variants associated with a particular phenotypic condition (i.e. a disease). Both genotyping and genome reconstruction plays a crucial role in this process. 

Genome reconstruction is a technique for inferring genomic ancestry. Since each DO mouse is a mosaic of eight founder genomes,  given a particular DO mouse, a high-density genotyping can trace each segments of  the DO genome back to one of the eight founders. This can be valuable in a number of aspects including identification of disease associated genetic variants, narrowing down a potential candidate gene list, understanding gene-gene interactions etc.

### 1.3 DO genome reconstruction process

* **Genotyping**: A genome-wide allelic profiling using either high-density microarrays or sequencing techniques.
* **Comparison to founder strain:** The process of comparing identified allelic profiles to the eight founder strains.
* **Haplotype inference:** A statistical approach for for identifying the most likely set of "haplotypes" or ancestry for each segment of a DO genome.
* **Identification of recombination points:** For each DO genome, identification of the genomic locations where the predicted most likely haplotype switch from one of the eight founders to another. 

This results of this process is a map of the DO mouse's genome, that approximate the most likely origins for each genomic segment. This map can then be used in association studies (e.g. GWAS) to identify regions of the genome associated with particular traits or diseases.

## 2. Genotyping dataset

### 2.1 GeneSeek FinalReport.txt

The `FinalReport.txt` represents the primary data collected through a GigaMUGA experiment. This ASCII formatted text file is produced by the Illumina Genome Studio software and is divided in two sections a) `Header` and b) `Data`. The header section provides the experimental metadata (e.g. date of the processing, number of detected SNPs, total number of samples etc). Genotyping information is captured in the data section in a tab delimited format. Each row in this section represents a SNP in a individual sample. 

```text
# snapshot of a FinalReport.txt file headers and column names
[Header]
GSGT Version    2.0.4
Processing Date 10/25/2022 12:27 PM
Content         GigaMuga_11769261_A.bpm
Num SNPs        143259
Total SNPs      143446
Num Samples     57
Total Samples   57
[Data]
SNP Name        Sample ID       Allele1 - Forward       Allele2 - Forward       Allele1 - Top   Allele2 - Top   Allele1 - AB    Allele2 - AB    GC Score        X       Y
```

#### 2.1.1 Description of the data columns

* `SNP Name`: Name of the SNP marker (e.g. AmpR002)
* `Sample ID`: An unique id that identify the sample, often matches with the Mouse id.
* `Allele-1 - Forward` & `Allele2 - Forward`: The diploid alleles (two different versions) for a given SNP on the forward (+) strand. When `Allele-1 - Forward` = `Allele2 - Forward`, the SNP location can be considered homozygous in the sample; else the location is heterozygous.
* `Allele1 - Top` & `Allele2 - Top`:   The allele calls on the top strand for the SNP, where the top strand refers to a consistent reference strand that is independent of the sample genome.
* `Allele1 - AB` & `Allele2 - AB`:  The allele calls translated to into a "A" or "B" format implemented in In Illumina's genotyping arrays. The two possible alleles for a SNP are designated as "A" (Ref) and "B" (Alt), which can correspond to any of the four nucleotides (A, C, G, or T) depending on the SNP.
* `GC Score`:  The "**GenCall**" score, a quality metric that ranges from 0 to 1. A higher score indicates greater confidence in the accuracy of the genotypes assigned to a SNP.
* `X` & `Y`: The normalized intensities for the A and B alleles, respectively.

### 2.2 LocusSummary.csv

The genotyping module included in Illumina GenomeStudio® produce a comprehensive report pertaining to each profiled SNP. This is commonly known as the '**LocusSummary**' file. This file provide quality matrices, genotype frequencies and other data for each location profiled through the genotyping array.  

#### 2.2.1 Description of the columns

* `Row`: The row index number.

* `Locus_Name`: The name of the profiled locus (i.e. the name of the marker or SNP).

* `Illumicode_Name`: The Illumina-specific identifier for the SNP (i.e Locus ID from the manifest file)

* `#No_Calls`: The number of samples for which a genotype could not be determined at this SNP (i.e. number of loci with GenCall scores below the call region threshold).

* `#Calls`: The number of samples for which a genotype was successfully determined at this SNP (i.e. number of loci with GenCall scores above the call region threshold).

* `Call_Freq`: 
  $$
  CallFreq = \frac{Calls}{(NoCalls + Calls)}
  $$
  
* `A/A_Freq`, `A/B_Freq`, `B/B_Freq`: Frequencies of each possible genotype at this SNP; where, `A/A`, `A/B`, and `B/B` represents the homozygous A, heterozygous and homozygous B alleles respectively. 

* `Minor_Freq`: The frequency of the minor allele at this SNP across all samples.

* `Gentrain_Score`: A score that reflects the reliability of the SNP genotyping. Higher values indicate more reliable genotyping.

* `50%_GC_Score`, `10%_GC_Score`: These columns provide information about the distribution of GenCall scores across samples. They represent the GenCall scores below which 50% and 10% of samples fall, respectively.

* `Het_Excess_Freq`: The frequency of heterozygous genotypes that exceed the expectation under Hardy-Weinberg equilibrium.

* `ChiTest_P100`: The p-value from a chi-square test of Hardy-Weinberg equilibrium.

* `Cluster_Sep`: A measure of the separation between clusters of different genotypes in the intensity data. Higher values indicate better separation and more reliable genotyping.

* `AA_T_Mean`, `AA_T_Std`, `AB_T_Mean`, `AB_T_Std`, `BB_T_Mean`, `BB_T_Std`: These columns provide the mean and standard deviation of the normalized theta (T) values for each genotype.

* `AA_R_Mean`, `AA_R_Std`, `AB_R_Mean`, `AB_R_Std`, `BB_R_Mean`, `BB_R_Std`: These columns provide the mean and standard deviation of the normalized R (signal intensity) values for each genotype.

* `Plus/Minus Strand`: This column specifies whether the SNP is found on the plus (+) or minus (-) strand of the DNA.

## 3. Getting started

### 3.1 Setting up project folder

* Clone this repository

  ```bash
  #!/usr/bin/bash
  # git --version
  git clone <repository link>
  ```

* Create required folders

  ```bash
  #!/usr/bin/bash
  cd GigaMUGA_pipeline
  mkdir -p data/samples
  mkdir results
  ```

* Arrange sample files

  ```bash
  #!/usr/bin/bash
  # move sample files to data/samples
  cp -r <path>/<EXP_Name>_FinalReport.txt data/samples
  ... # follow same process if more than one *FinalReport.txt
  cp -r <path>/<EXP_Name>_LocusSummary.csv data/samples
  ... # follow same process if more than one *LocusSummary.csv
  ```

### 3.2 Get required external data 

#### 3.2.1 Collaborative Cross primary dataset

GigaMUGA SNP data for the CC founder lines are made available by the Jackson Laboratories (originally hosted at ftp.jax.org/MUGA). The [zip file](https://figshare.com/articles/dataset/GM_primary_files_zip/5404675/1) contains the following primary .RData files:

* `GM_code`: strain name -> genotype code (115 lines)

* `GM_geno`: genotypes (markers x strains) as single-letter codes N,C,T,A,G,H (143259 x 115)

* `GM_sex`: sex of individuals (M/F)

* `GM_snps`: snps x info (marker, chr, pos, cM, A1, A2, type, ..., tier); chr = 1-19,X,Y,M,P

* `GM_x`, `GM_y`: x and y intensities, markers x strains

The extracted zip file produces a folder `GM_primary_files`. The scripts provided in this repository expects this folder inside the `data/`folder under the project directory.

#### 3.2.2 GigaMUGA founder dataset for R/QTL2

Processed GigaMUGA SNP genotype data on the Collaborative Cross founders, useful for preparing Diversity Outbred mouse genotypes for use with R/qtl2 (http://kbroman.org/qtl2). The [zip file](https://figshare.com/articles/dataset/GM_processed_files_zip/5404759) contains the following datasets:

* `GM_allelecodes.csv`: allele codes used for all SNPs

* `GM_foundergeno*.csv`: consensus genotypes for each of the 8 founder lines

* `GM_gmap*.csv`: Genetic map positions of the SNPs

* `GM_info.csv`: marker, chr, Mbp position, cM position, and tier

* `GM_pmap*.csv`: Physical map positions (in Mbp) of the SNPs

The extracted zip file produces a folder `GM_processed_files`. The scripts provided in this repository expects this folder inside the `data/`folder under the project directory.

## 4. Genome Reconstruction pipeline

### 4.1 Parse genotyping data 

```bash
#!/usr/bin/bash
# example
Rscript src/parse_geneseek_FinalReport.R --help
```

**Command line options:**

| Option              | Description                                                | Default                               |
| ------------------- | ---------------------------------------------------------- | ------------------------------------- |
| `-i` or `--input`   | Path to the folder with the `FinalReport` file             | `data/samples/`                       |
| `-o` or `--out_dir` | Path to the folder where the parsed output will be written | `results/parse_geneseek_FinalReport/` |

### 4.2 Genotyping quality control

```bash
#!/usr/bin/bash
# example
Rscript src/genotype_qc.R --help
```

**Command line options:**

| Option                    | Description                                                  | Default                                                     |
| ------------------------- | ------------------------------------------------------------ | ----------------------------------------------------------- |
| `-i` or `--input`         | Path to the parsed `FinalReport` file from the previous step | `results/parse_geneseek_FinalReport/parsed_FinalReport.txt` |
| `-l` or `--locus_summary` | Path to the folder with the `LocusSummary.csv` file          | `data/samples/`                                             |
| `-o` or `--out_dir`       | Path to the folder where the parsed output will be written   | `results/genotype_qc/`                                      |

### 4.3 Convert genotype data to qtl2 acceptable format

```bash
#!/usr/bin/bash
# example
Rscript src/geneseek2qtl2.R --help
```

**Command line options:**

| ption                                 | Description                                                  | Default                               |
| ------------------------------------- | ------------------------------------------------------------ | ------------------------------------- |
| `-i` or `--input`                     | Path to the genotype data after QC                           | `results/genotype_qc/Genotype_QC.csv` |
| `-g` or  `--gigaMuga_processed_files` | Standard GigaMUGA processed zipped and provided through qtl2 | `data/GM_processed_files/GM/`         |
| `-o` or `--out_dir`                   | Path to the folder where the parsed output will be written   | `results/genotype_qc/`                |

## References

1. GenomeStudio® Genotyping Module v2.0 Software Guide- Illumina; [Link](https://support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/genomestudio/genomestudio-2-0/genomestudio-genotyping-module-v2-user-guide-11319113-01.pdf)
1. Preparing Diversity Outbred (DO) mouse data for R/qtl2 - Karl Broman; [Link](https://kbroman.org/qtl2/pages/prep_do_data.html)
1. 

## Acknowledgement
