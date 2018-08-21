# Miso Example


## 1. Description

An efficient approach to fish out isotopically labeled analyte.

## 2. Download

You can clone or download [this GitHub repositor](https://github.com/YonghuiDong/Miso_example.git), extract it and run the example in your own computer.

## 3. Usage

```r

##(1) install stable version of Miso package
install.packages("Miso")
library(Miso)

##(2) load data files, 2 files are provided in this example, one is an xcms pre-processed data, 
## xset_g_r_g.rda; and the other is a peak table.
load(file = 'data/xset_g_r_g.rda')
load(file = 'data/lcms.rda')

##(3) deisotoping and/or deadducting (optional but recommended)
library('CAMERA')
an <- xsAnnotate(xset_g_r_g)
an <- groupFWHM(an)
an <- findIsotopes(an, maxcharge = 3)
peaklist <- getPeaklist(an)
peaklist$isotopes <- sub("\\[.*?\\]", "", peaklist$isotopes)
peaklist <- peaklist[peaklist$isotopes == '' | peaklist$isotopes == '[M]+', ]

##(4a) First filtering: fast
explist <- prefilter(lcms)

##(4b) Alternative first filering method: more complete, but slow
explist <- prefilter2(lcms)

##(5) Second filtering
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

##(6) Generate results
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
