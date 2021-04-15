# ISME.Turnham.et al.2021
Repository of supplemetal tables, figures, materials, and scripts for the paper "Mutualistic microalgae co-diversify with reef corals that acquire symbionts during egg development" Kira E. Turnham, Drew C. Wham, Eugenia Sampayo, and Todd C. LaJeunesse. (In review, 2021 ISME-J)

doi: TBD
<br />
## Files overview

### Supplemental Tables
* **Table_S1_samples.pdf** List of *Cladocopium* samples with corresponding species names, ITS2 designations, host *Pocillopora* species, country of
origin, collection site, and approximate latitude/longitude coordinates.

* **Supplemental_Tables_S2-S4.pdf** 
* * S2. PCR primers and associated references for mitochondrial, nuclear, chloroplast genes and non-coding regions
* * S3. *Cladocopium* microsatellite loci used in this study and corresponding allele size ranges, repeat motifs, annealing temperatures, and references
* * S4. Cell size differences between *Cladocopium pacificum* and *C.latusorum* from *Pocillopora* host samples obtained in Palau. The mean length and width of ovate cells are listed based off of individual measurements of 55-110 cells per sample

## Supplemental Figures
* **Supplemental_Figures_S1_S2.pdf** 
* * S1. Genotype accumulation boxplot showing the percent of unique multilocus genotypes recovered for A) *Cladocopium latusorum* and B) *C. pacificum* based on increasing numbers of loci subsampled from the total dataset of 8 loci. Bars represent 95% confidence intervals.
* * S2. Multi-locus genotypes analyzed via STRUCTURE plot showing no overlap between *C. latusorum* and *C. pacificum* at K=3 and K=4. Partitioning of each species to these higher K values shows within-species subdivision. Further analyses, beyond the scope of this work, are needed to more accurately assess within-species population genetic structure.

### Sequence alignments
* **Supplemental_Dataset_S_concatenated_alignment.nex** *Cladocopium* concatenated nexus alignment: ITS2 LSU cp23S cob cox1
* **Supplemental_Dataset_psbA_alignment.nex** *Cladocopium* psbA nexus alignment
* **psba_alignment_mol_clock.nex** *Cladocopium* psbA nexus alignment used in molecular clock analyses (contains *Cladocopium* associated with sibling species *Porites panamensis* from Eastern Pacific and *Porites porites var. colonensi* from Western Atlantic

### MrModelTest output used in downstream Bayesian analyses
* **mr_model_test_output.txt** MrModelTest output containing Akaike Information Criterion

### Multi-locus genotypes and R scripts for microsatellite analyses
Allele freqeuncy, STRUCUTRE, and t-SNE plots were created using PlotSTR, BayesAllele, and t-SNE scripts developed by Drew Wham and available at https://github.com/DrewWham/theclonalescnet (see also Wham et al., 2016).
* **Supplemental_Dataset_microsatellites.docx** *C. pacificum* and *C. latusorum* multi-locus genotypes
* **Supplemental_Dataset_linkage.txt** Linkage disequilibrium analysis between pairs of all microsatellite loci used for multi-locus genotyping performed using Genepop on the Web https://genepop.curtin.edu.au/)
* **PCA_TSNE_and_GMM_script.R** Script for performing principal components analysis (PCA), t-distributed stochastic neighbor embedding (t-SNE), and Gaussian Mixture Model (GMM) clustering analyses
* **bayesallele_frequency.Rmd** Script for allele frequency based on BayesAllele 
* **structureplot_K2.Rmd** Script for STRUCTURE using PlotSTR
