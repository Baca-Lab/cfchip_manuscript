The Rmd file in this repo will generate most of the figures for the cfChIP manuscript

To run this, two types of data are needed that are not included in the repo:
1. RDS files containing granges objects for each sample. These contain the fragment locations. These go in data/granges, eg:
data/granges/Ac_mPC15_3R.RDS
data/granges/Ac_mPC15_all.RDS
data/granges/Ac_mPC16_1R.RDS
data/granges/Ac_mPC17_3R.RDS

2. Peak files containing the MACS2-called peaks for each sample. These go in data/peaks, eg:
data/peaks/Ac_mPC15_3R.rep1_sorted_peaks.narrowPeak.bed
data/peaks/Ac_mPC15_all.rep1_sorted_peaks.narrowPeak.bed
data/peaks/Ac_mPC16_1R.rep1_sorted_peaks.narrowPeak.bed
data/peaks/Ac_mPC17_3R.rep1_sorted_peaks.narrowPeak.bed

