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

Running `remotes::install_github('avril-m-harder/PAVplotR')` will install `PAVplotR` and required dependencies if not already installed: `vcfR`, `ggplot2`, `reshape2`, `dplyr`, `ggnewscale`, and `ggtext`.

### 3. Testing it with example data
