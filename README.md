# Wavelet_analysis

Wavelet_analysis is a repository of python scripts to carry out spectral decomposition
of tiemseries using a set of wavelet basis functions. The methodology is discussed extensively in [Terrence and Compo et al.] (https://psl.noaa.gov/people/gilbert.p.compo/Torrence_compo1998.pdf).


## Files

wavelet_lib.py - script containing functions to compute the wavelet transform, power spectra and significance tests on a timeseries.

wavelet_plot.py - functions to plot wavelet spectra and cross spectra

wavelet_analysis.py - script to load data, call wavelet functions and output spectra

PDO_wavelet_UKESM_SON_DJF_5_yr_window.png - sample output of wavelet power spectra from climate model data
(PDO_wavelet_UKESM_SON_DJF_5_yr_window.png)
cross_coherence_SSWs_PDO_SON_UKESM.png - sample output of wavelet power spectra from climate model data


## Usage

This repository has been used for analysis in [Dimdore-Miles et al. 2020](https://wcd.copernicus.org/preprints/wcd-2020-56/) which analyses long term variability in stratospheric circulation indices
