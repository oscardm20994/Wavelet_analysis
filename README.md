# Wavelet_analysis

Wavelet_analysis is a repository of python scripts to carry out spectral decomposition
of tiemseries using a set of wavelet basis functions. The methodology is 
discussed extensively in [Terrence and Compo et al.](https://psl.noaa.gov/people/gilbert.p.compo/Torrence_compo1998.pdf)

## Mathematical Approach

The wavelet transform of a uniform 1-dimensional time series, x, of length N and timestep &delta;t is given by the convolution between the series and a scaled and translated version of a wavelet function &psi;<sub>0

<img src="https://render.githubusercontent.com/render/math?math=W_n(s) = \sum^{N - 1}_{n' = 0} x_{n'} \psi^* \bigg[(n' - n) \frac{\delta t}{s}\bigg],">

where s is the wavelet scale indicating the frequency of the wavelet. Varying s and translating along the time scale (the index n), W<sub>n</sub> indicates the amplitude of signals at different scales and their variation in time. The scale s is increased in powers of 2 according to 

<img src="https://render.githubusercontent.com/render/math?math=s_j = s_0 2^{j \deltaj},          j = 0, 1, ..., J">

<img src="https://render.githubusercontent.com/render/math?math=J = \delta j^{-1} log_2\bigg(\frac{N \delta t}{s_0}\bigg),">

where s<sub>0</sub> is the shortest resolvable scale of a signal, J corresponds to the longest and &delta;j is the scale resolution. The translated and scaled wavelet has the form

<img src="https://render.githubusercontent.com/render/math?math=\psi^* \bigg[(n' - n) \frac{\delta t}{s}\bigg] = \bigg(\frac{\delta t}{s}\bigg)^{1/2} \psi_0\bigg[(n' - n) \frac{\delta t}{s}\bigg]">

a common form of \psi<sub>0</sub> is as a Morlet wavelet, an oscillatory function enveloped by a Gaussian which is expressed as

<img src="https://render.githubusercontent.com/render/math?math=\psi_0(p) = \pi^{-1/4} e^{i\omega_0 p} e^{\frac{p^2}{2}}.">

It is computationally quicker to compute the wavelet transform in discrete Fourier space. By the convolution theorem, the transform reduces to multiplication

<img src="https://render.githubusercontent.com/render/math?math=W_n(s) = \sum^{N - 1}_{k = 0} \hat{x}_{k} \hat{\psi}^* (s\omega_k) e^{i \omega_k n \delta t},">

where $\hat{x}_{k}$ and $\hat{\psi}$ are the discrete Fourier transforms of the time series $x$ (equation \ref{fourier1}) and the wavelet function (equation \ref{fourier2}) respectively,

<img src="https://render.githubusercontent.com/render/math?math=\hat{x}_k = \frac{1}{N} \sum^{N-1}_{n = 0} x_n e^{\frac{-2\pi i k n}{N}}">

<img src="https://render.githubusercontent.com/render/math?math=\hat{\psi}(s\omega_k) = \bigg(\frac{2 \pi s}{\delta t}\bigg) \pi^{-1/4}H(\omega_k) e^{-(s\omega_k - \omega_0)^2/2}.">

H(&omega<sub>k</sub>;) is the Heavyside function. The square modulus of the wavelet transform gives the wavelet power spectrum which indicates relative strength of signals in the time series as a function of signal period and discretised time. We also define a confidence interval for wavelet power observed at a given period and time for a series by assuming a mean background spectrum corresponding to that of a first order autoregressive (AR1, red noise) process modelled by

<img src="https://render.githubusercontent.com/render/math?math=x_n = \alpha x_{n - 1} + z_n,">

where &alpha: is the lag-1 autocorrelation of the time series and z<sub>n is Gaussian white noise. such a process's wavelet power spectrum is &chi:<sup>2</sup> distributed and therefore can be used to define a 95\% confidence interval for any observed power. Additionally, the cross wavelet spectrum of two time series x and y with associated wavelet spectra Wx<sub>n and Wy<sub>n gives a measure of coincident power (the same period at the same timepoints) between the series. It is given by

<img src="https://render.githubusercontent.com/render/math?math=\vert W^{xy}_n(s)\vert = \vert W^{x*}_n(s) W^{y}_n(s)\vert,">

where Wx<sub>n</sub> is the complex conjugate of the wavelet power spectrum of x. The complex argument of Wxy<sub>n</sub> gives the local phase difference between signals in x and y in frequency-time space. The phase relationship between the two time-series can be represented by a
vector that subtends an angle representing the phase difference: On all plots of cross spectra, arrows to the right (left) denoted signals which are in-phase and correlated (anti-correlated). Vertical arrows indicate a phase relationship of &pi;/2 between the time-series, so that the evolution of
one is correlated with the rate-of-change of the other. As for individual power spectra, we define a confidence interval for which cross power of a larger amplitude is deemed significant (>95% confidence interval) by comparing power exhibited by actual series with a theoretical red noise process. The cross power of two such AR1 processes is theoretically distributed such that the probability of obtaining cross power greater than a set of red-noise processes is

<img src="https://render.githubusercontent.com/render/math?math=D\bigg(\frac{\vert W^{xy}_n(s)\vert}{\sigma_x \sigma_y} < p\bigg) = \frac{Z_\nu(p)}{\nu} \sqrt{P^x_k P^y_k},">


where &sigma; denotes the standard deviation of the time series, Z is the confidence interval defined by p (Z = 3.999 for 95% confidence), &nu; is the degrees of freedom for a real wavelet spectrum (&nu; = 2) and P<sub>xk is the theoretical Fourier spectrum of the AR1 process. For a given wavenumber k, this can be expressed as

<img src="https://render.githubusercontent.com/render/math?math=P_k = \frac{1 - \alpha^2}{\vert 1 - \alpha e^{2i\pi k} \vert^2}.">

## Files

wavelet_lib.py - script containing functions to compute the wavelet transform, power spectra and significance tests on a timeseries.

wavelet_plot.py - functions to plot wavelet spectra and cross spectra

wavelet_analysis.py - script to load data, call wavelet functions and output spectra

[PDO_wavelet_UKESM_SON_DJF_5_yr_window.png](PDO_wavelet_UKESM_SON_DJF_5_yr_window.png) - sample output of wavelet power spectra from climate model data

[cross_coherence_SSWs_PDO_SON_UKESM.png](cross_coherence_SSWs_PDO_SON_UKESM.png) - sample output of wavelet cross power spectra from climate model data.


## Usage

This repository has been used for analysis in [Dimdore-Miles et al. 2020](https://wcd.copernicus.org/preprints/wcd-2020-56/) which analyses long term variability in stratospheric circulation indices.
