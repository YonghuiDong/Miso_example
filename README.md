# Miso Example


## 1. Description

An R package to fish out isotopically labeled analytes in single, dual or multiple isotope labeling experiment.

**Features**

(1) **Fully automated**: Miso uses xcms object as input, and automatically processes the data by taking advantage of the phenoData from xcms object.

(2) **Flexible**: It enables detecting molecules labeled with various biologically relevant stable isotopes such as hydrogen (2H), carbon (13C), oxygen (18O), nitrogen (15N) and sulfur (34S). It can be applied to single, dual and multiple isotope labeling experiments without modifing any scripts. In addition, different algorithms are used to process the data with and without replicates.

(3) **Efficient**: fast

(4) **User-friendly**: Easy to use. Three different outputs are provided. (i): a full data table list, which includes a full list of isotopologues; (ii) a reduced data table list, in which only the base peaks of each isotopologues are included; (iii) the interactive plot, which allows the user to visually check the result.

## 2. Experiment Design

Figure 1 describes **single**, **dual** and **multiple** stable isotope labeling experiment. Detailed explainations regarding the experiment design and advantages of dual labeling experiments can be found in ([Liron et al., 2009](https://pubs.acs.org/doi/10.1021/ac901495a); [Liron et al., 2018](https://pubs.acs.org/doi/10.1021/acs.analchem.8b01644); [Dong et al., 2019](https://doi.org/10.1093/bioinformatics/btz092))

<p align="center"> 
<img src="Image/workflow.png" width="600">
</p>

<p align="center">
<b> Figure 1.</b>  Experimental set up for single, dual and multiple stable isotope labeling experiments
</p>

## 3. The Example Dataset

The example data (`Data/lcms.rda`) provided in this example is from a dual labeling experiment. It consists of 5 groups: Group A (SP medium, control group), Group B (0.5 mg mL<sup>−1</sup> unlabeled tyrosine), Group C (0.5 mg mL<sup>−1</sup> tyrosine D4, Label I), Group D (0.5 mg <sup>mL−1</sup> tyrosine 13C915N1, Label II), and Group E (0.5 mg mL<sup>−1</sup> Label I + Label II). 

The Groups B, C and D are shown in Figure 1. More details regarding the this dataset and experimental details can be found in ([Liron et al., 2018](https://pubs.acs.org/doi/10.1021/acs.analchem.8b01644)).


## 4. Tutorial

**Attention**: This tutorial is based on the Miso package version 0.2.0, which is slightly different from the one described in the publication ([Dong et al., 2019](https://doi.org/10.1093/bioinformatics/btz092)), which was based on Miso version 0.1.5.

```r

##(1) install stable version of Miso package
install.packages("Miso")
library(Miso)

##(2) Example 1: single labeling experiment.

##(2a) First filtering: fast. 
## This approach is suitable for the the situation that the intra-sample variation is large and/or there are no replicates.
## use ?prefilter() to check all the parameters in this function

explist <- prefilter(lcms)

##(2b) Alternative first filering method: more complete, but slow. 
## This approach is suitable for samples with replicates. 
## use ?prefilter2() to check all the parameters in this function

explist <- prefilter2(lcms)

##(3) Second filtering
## Group C was fed with H2
## Here we are interested in detecting molecules labeled with 4, 3 or 2 H2 (deuterium). 
## n11 = 4, n12 = 2.

exp.B <- explist$exp.B
exp.C <- explist$exp.C
exp.D <- explist$exp.D
iso.C <- diso(iso1 = 'H2', n11 = 4, n12 = 2, exp.base = exp.B, exp.iso = exp.C)

## Group D was fed with C13, and N15
## Here we are interested in detecting molecules labeled with 9, 8, 7 or 6 C13, 
## and 1 or 0 N15 (n11 = 9, n12 = 6 for C13, and n21 = 1, n22 = 0 for N15)

iso.D <- diso(iso1 = 'C13', n11 = 9, n12 = 6, iso2 = 'N15', n21 = 1, n22 = 0,
              exp.base = iso.C[, 1:3], exp.iso = exp.D)

##(4) Generate results
## Two types of results are provided. A Full list and a reduced list which contains only 
## one form of labeled molecules.

full_Result <- Fresult(iso.C, iso.D)
reduced_Result <- Rresult(full_Result, cutint = 4000)

##(5) plot result
## view the first row of Full_result
isoplot(full_Result, 1)
```

## 4. Attention    

1. R memory limit error may appear during data processing especially for high resolution dataset:   

`Error: memory exhausted (limit reached?), Error during wrapup: memory exhausted (limit reached?)` 

This error is due to the following script:

```r
## In group C, we are looking for analytes labeled with 5, 4, or 3 deuteriums (H2).
iso.C <- diso(iso1 = H2, n11 = 4, n12 = 2, exp.base = exp.B, exp.iso = exp.C)
```

To solve this memory limit problem, the above script can be decomposed into 3 sub-scripts, which respectively search for analytes labled with 4, 3, and 2 deuteriums (H2).

```r
iso.C5 <- diso(iso1 = 'H2', n11 = 4, n12 = 4, exp.base = exp.B, exp.iso = exp.C)
iso.C4 <- diso(iso1 = 'H2', n11 = 3, n12 = 3, exp.base = exp.B, exp.iso = exp.C)
iso.C3 <- diso(iso1 = 'H2', n11 = 2, n12 = 2, exp.base = exp.B, exp.iso = exp.C)

## The results are then combined as iso.C:
iso.C <- rbind(iso.C5, iso.C4, iso.C3)
```

The decomposition step is only usually necessasy for iso.C, as the result list has been significantly reduced. we do not have to do it again for iso.D.

2. Optionally, users can also perform deisotoping before running Miso.

```r
##(2) deisotoping and/or deadducting (optional)
library('CAMERA')
an <- xsAnnotate(xset_g_r_g_fill)
an <- groupFWHM(an)
an <- findIsotopes(an, maxcharge = 3)
peaklist <- getPeaklist(an)
peaklist$isotopes <- sub("\\[.*?\\]", "", peaklist$isotopes)
peaklist <- peaklist[peaklist$isotopes == '' | peaklist$isotopes == '[M]+', ]
```

## 5. Session Information

```r
sessionInfo()

R version 3.3.3 (2017-03-06)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: macOS  10.13.6

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] Miso_0.1.3

loaded via a namespace (and not attached):
 [1] magrittr_1.5        magick_1.6          parallel_3.3.3      tools_3.3.3        
 [5] yaml_2.1.14         Rcpp_0.12.14        xml2_1.1.1          stringi_1.1.5      
 [9] S4Vectors_0.12.2    knitr_1.16          stringr_1.2.0       BiocGenerics_0.20.0
[13] stats4_3.3.3   
```
