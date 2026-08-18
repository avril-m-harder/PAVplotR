## PAVplotR: an R package for plotting PAV from a haploid VCF
Accepts as input a haploid VCF that must include the reference haplotype as a sample. For a specified bin size, calculates PAV in an expanded coordinate system to account for insertions relative to the reference and plots those bins with the x-axis (chromosome position) labeled according to the reference haplotype's original coordinate system. Tested on fairly small regions so far (≤ 500 kb).

## 
### 1. Quick start
#### 1.1 Input formatting
All you need to plot PAV is a haploid VCF (e.g., output from *vg deconstruct* or *SyRI* with a bit of reformatting as necessary). There are a few requirements for this VCF:
* The reference haplotype should be included as a sample with all reference alleles.
* 
