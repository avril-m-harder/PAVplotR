## PAVplotR: an R package for plotting PAV from a haploid VCF
Accepts as input a haploid VCF that must include the reference haplotype as a sample. For a specified bin size, calculates PAV in an expanded coordinate system to account for insertions relative to the reference and plots those bins with the x-axis (chromosome position) labeled according to the reference haplotype's original coordinate system. Tested on fairly small regions so far (≤ 500 kb, ~50 haplotypes).

## 
### 1. Input formatting
All you need is a haploid VCF (e.g., output from *vg deconstruct* or *SyRI* with a bit of reformatting as necessary). There are a few requirements for this VCF:
* The reference haplotype should be included as a sample with all reference alleles.
* VCF records must not be multiallelic. You can [split multiallelic records](https://samtools.github.io/bcftools/bcftools.html#norm) with `bcftools norm -m -any`.
* The first nucleotide of the REF and ALT alleles within each record should match. A mismatch will result in a warning during calculation of the expanded coordinate system, but calculation will proceed.

Converting multi-haplotype graphs to VCF format can sometimes produce unexpected results. For regions with complex variants, I would recommend comparing PAV output against dotplots of the original input haplotypes to check for inconsistencies.

### 2. Installing and dependencies

Installing `pak` and running `pak::pak('avril-m-harder/PAVplotR')` will install `PAVplotR` and required dependencies if not already installed (`vcfR`, `ggplot2`, `reshape2`, `dplyr`, `ggnewscale`, and `ggtext`).

### 3. Testing it with example data
The below steps will produce 2 plots for a simplified view of the [SHATTERING1 locus](https://www.nature.com/articles/s41586-026-10229-9/figures/3) in 33 sorghum haplotypes[^1].

First, install PAVplotR, load necessary libraries, and set the minimum information required to plot a region of interest, including the bin size for which you'd like to calculate proportional PAV.
```
pak::pak('avril-m-harder/PAVplotR')

library(PAVplotR)
library(vcfR)
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggnewscale)
library(ggtext)

## provide minimum information on the ROI
start.pos <- 12398577
end.pos <- 12404976
ref <- 'BTx623'
roi <- 'SH1'
chr <- 'Chr01'

## set bin size for calculating PAV proportions
bin.size <- 100
```

Using built-in example data, calculate the extended coordinate system to accommodate insertions relative to the reference haplotype with `build_coordsystem()`. Currently, all variants occurring at the same reference coordinate position will be vertically stacked. For example, for 2 distinct insertions and 1 deletion at the same locus, the 2 insertions will be stacked and accounted for with extended coordinates, whereas the deletion will be visualized just following these two insertions where the reference coordinate system picks back up. Calculate PAV bins for each haplotype with `calculate_bins()`.
```
vcf <- PAVplotR::example_vcf

coord.map <- build_coordsystem(vcf,
                               start_pos = start.pos,
                               end_pos = end.pos)

bin.dat <- calculate_bins(vcf,
                          coord_map = coord.map,
                          bin_size = bin.size,
                          start_pos = start.pos,
                          end_pos = end.pos)
```
The first plot produced with `plot_pav()` will plot simple presence-absence for the ROI.
```
plot_pav(coord_map = coord.map,
         presence_matrix = bin.dat$matrix,
         bin_info = bin.dat$bin_info,
         output_fmt = 'tiff',
         roi = roi,
         ref_hap = ref,
         chrom = chr,
         region_start = start.pos,
         region_end = end.pos,
         bin_size = bin.size,
         hap_order = 'refdist',
         width = 8, height = 5)
```
![basic PAV plot](https://github.com/avril-m-harder/PAVplotR/blob/main/man/figures/SH1_BTx623_Chr01_12398577_12404976_100_PAV.tiff)

The second plot produced with `plot_pav_hiliteInsertions()` will plot the same presence-absence information for the ROI, but will highlight insertions relative to the reference. A threshold for highlighting bins containing insertions can be specified with `ins_thresh` (the number of additional bases in a bin, relative to the reference haplotype, to flag the haplotype and bin as containing ≥1 insertion).
```
plot_pav_hiliteInsertions(coord_map = coord.map,
                          presence_matrix = bin.dat$matrix,
                          bin_info = bin.dat$bin_info,
                          output_fmt = 'tiff',
                          roi = roi,
                          ref_hap = ref,
                          chrom = chr,
                          region_start = start.pos,
                          region_end = end.pos,
                          bin_size = bin.size,
                          ins_thresh = 30,
                          color_high = 'darkseagreen3',
                          ins_color_high = 'thistle4',
                          hap_order = 'refdist',
                          width = 8, height = 5)

```



[^1]: Morris, G.P., Harder, A.M., Healey, A.L., McLaughlin C.M., Rifkin, J.L. *et al*. A sorghum pangenome reference improves global crop trait discovery. *Nature* **652**, 1245–1253 (2026). [https://doi.org/10.1038/s41586-026-10229-9]
