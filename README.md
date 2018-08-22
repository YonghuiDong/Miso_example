# Miso Example


## 1. Description

An efficient approach to fish out isotopically labeled analyte.

## 2. Usage

```r

##(1) install stable version of Miso package
install.packages("Miso")
library(Miso)

##(2a) First filtering: fast. This approach is suitable for the the situation that the 
## intra-sample variation is large and/or there are no replicates
explist <- prefilter(lcms)

##(2b) Alternative first filering method: more complete, but slow. This approach is suitable 
## for samples with replicates. 
explist <- prefilter2(lcms)

##(3) Second filtering
## Group C was fed with H2
## Here we are interested in detecting molecules labeled with 4, 3 or 2 H2 (deuterium). 
## n11 = 4, n12 = 2.

exp.B <- explist$exp.B[, -2]
exp.C <- explist$exp.C[, -2]
exp.D <- explist$exp.D[, -2]
iso.C <- diso(iso1 = 'H2', n11 = 4, n12 = 2, exp.base = exp.B, exp.iso = exp.C)

## Group D was fed with C13, and N15
## Here we are interested in detecting molecules labeled with 9, 8, 7 or 6 C13, 
## and 1 or 0 N15 (n11 = 9, n12 = 6 for C13, and n21 = 1, n22 = 0 for N15)

iso.D <- diso(iso1 = 'C13', n11 = 9, n12 = 6, iso2 = 'N15', n21 = 1, n22 = 0,
              exp.base = iso.C[,1:2], exp.iso = exp.D)

##(4) Generate results
## Two types of results are provided. A Full list and a reduced list which contains only 
## one form of labeled molecules.

full_Result <- Fresult(iso.C, iso.D)
reduced_Result <- Rresult(full_Result)
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
