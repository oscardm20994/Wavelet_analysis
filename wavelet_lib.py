# !/usr/bin/python
# WAVELET LIBRARY - Based on Wavelet methods outlined in Torrence and Combo (1998)
# Examples of use of this library

import numpy as np
import math
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
import os
import sys
import netCDF4


""" Translating mfiles of the Torrence and Combo to python functions
    1 - wavetest.m
    2 - wave_bases.m
    3 - wave_signif.m
    4 - chisquare_inv.m
    5 - chisquare_solve.m
"""


def nextpow2(i):
    n = 2
    while n < i:
        n = n * 2
    return n


def wave_bases(mother, k, scale, param):
    """Computes the wavelet function as a function of Fourier frequency
    used for the CWT in Fourier space (Torrence and Compo, 1998)
    -- This def is called automatically by def wavelet --

    _____________________________________________________________________
    Inputs:
    mother - a string equal to 'Morlet'
    k      - a vectorm the Fourier frequecies
    scale  - a number, the wavelet scale
    param  - the nondimensional parameter for the wavelet function

    Outputs:
    daughter       - a vector, the wavelet function
    fourier_factor - the ratio os Fourier period to scale
    coi            - a number, the cone-of-influence size at the scale
    dofmin         - a number, degrees of freedom for each point in the
                     wavelet power (Morlet = 2)

    Call function:
    daughter,fourier_factor,coi,dofmin = wave_bases(mother,k,scale,param)
    _____________________________________________________________________
    """
    n = len(k)  # length of Fourier frequencies (came from wavelet.py)
    """CAUTION : default values"""
    if (mother == 'Morlet'):  # choose the wavelet function
        param = 6  # For Morlet this is k0 (wavenumber) default is 6
        k0 = param
        # table 1 Torrence and Compo (1998)
        expnt = -pow(scale * k - k0, 2) / 2 * (k > 0)
        norm = math.sqrt(scale * k[1]) * \
            (pow(math.pi, -0.25)) * math.sqrt(len(k))
        daughter = []  # define daughter as a list

        for ex in expnt:  # for each value scale (equal to next pow of 2)
            daughter.append(norm * math.exp(ex))
        k = np.array(k)  # turn k to array
        daughter = np.array(daughter)  # transform in array
        daughter = daughter * (k > 0)  # Heaviside step function
        # scale --> Fourier
        fourier_factor = (4 * math.pi) / (k0 + math.sqrt(2 + k0 * k0))
        # cone-of- influence
        coi = fourier_factor / math.sqrt(2)
        dofmin = 2  # degrees of freedom
# ---------------------------------------------------------#
    elif (mother == 'DOG'):
        param = 2
        m = param
        expnt = -pow(scale * k, 2) / 2.0
        pws = (pow(scale * k, m))
        pws = np.array(pws)
        """CAUTION gamma(m+0.5) = 1.3293"""
        norm = math.sqrt(scale * k[1] / 1.3293) * math.sqrt(n)
        daughter = []
        for ex in expnt:
            daughter.append(-norm * pow(1j, m) * math.exp(ex))
        daughter = np.array(daughter)
        daughter = daughter[:] * pws
        fourier_factor = (2 * math.pi) / math.sqrt(m + 0.5)
        coi = fourier_factor / math.sqrt(2)
        dofmin = 1
# ---------------------------------------------------------#
    elif (mother == 'PAUL'):  # Paul Wavelet
        param = 4
        m = param
        k = np.array(k)
        expnt = -(scale * k) * (k > 0)
        norm = math.sqrt(scale * k[1]) * \
        (2 ** m / math.sqrt(m * \
                            (math.factorial(2 * m - 1)))) * math.sqrt(n)
        pws = (pow(scale * k, m))
        pws = np.array(pws)
        daughter = []
        for ex in expnt:
            daughter.append(norm * math.exp(ex))
        daughter = np.array(daughter)
        daughter = daughter[:] * pws
        daughter = daughter * (k > 0)     # Heaviside step function
        fourier_factor = 4 * math.pi / (2 * m + 1)
        coi = fourier_factor * math.sqrt(2)
        dofmin = 2
    else:
        print ('Mother must be one of MORLET,PAUL,DOG')

    return daughter, fourier_factor, coi, dofmin


def wavelet(Y, dt, param, dj, s0, j1, mother):
    """Computes the wavelet continuous transform of the vector Y,
       by definition:

    W(a,b) = sum(f(t)*psi[a,b](t) dt)        a dilate/contract
    psi[a,b](t) = 1/sqrt(a) psi(t-b/a)       b displace

    Only Morlet wavelet (k0=6) is used
    The wavelet basis is normalized to have total energy = 1 at all scales

    _____________________________________________________________________
    Input:
    Y - time series
    dt - sampling rate
    mother - the mother wavelet function
    param - the mother wavelet parameter

    Output:
    ondaleta - wavelet bases at scale 10 dt
    wave - wavelet transform of Y
    period - the vector of "Fourier"periods ( in time units) that correspond
             to the scales
    scale - the vector of scale indices, given by S0*2(j*DJ), j =0 ...J1
    coi - cone of influence

    Call function:
    ondaleta, wave, period, scale, coi = wavelet(Y,dt,mother,param)
    _____________________________________________________________________

    """

    n1 = len(Y)  # time series length
    #s0 = 2 * dt  # smallest scale of the wavelet
    # dj = 0.25  # spacing between discrete scales
    # J1 = int(np.floor((np.log10(n1*dt/s0))/np.log10(2)/dj))
    J1 = int(np.floor(np.log2(n1 * dt / s0) / dj))  # J1+1 total os scales
    # print 'Nr of Scales:', J1
    # J1= 60
    # pad if necessary
    x = detrend_mean(Y)  # extract the mean of time series
    pad = 1
    if (pad == 1):
        base2 = nextpow2(n1)  # call det nextpow2
    n = base2
    """CAUTION"""
    # construct wavenumber array used in transform
    # simetric eqn 5
    #k = np.arange(n / 2)
    
    import math
    k_pos, k_neg = [], []
    for i in arange(0, int(n / 2) ):
        k_pos.append(i * ((2 * math.pi) / (n * dt)))  # frequencies as in eqn5
        k_neg = k_pos[::-1]  # inversion vector
        k_neg = [e * (-1) for e in k_neg]  # negative part
        # delete the first value of k_neg = last value of k_pos
        #k_neg = k_neg[1:-1]
    print(len(k_neg),len(k_pos))    
    k = np.concatenate((k_pos, k_neg), axis=0)  # vector of symmetric
    # compute fft of the padded time series
    f = np.fft.fft(x, n)
    scale = []
    for i in range(J1 + 1):
        scale.append(s0 * pow(2, (i) * dj))

    period = scale
    # print period
    wave = np.zeros((J1 + 1, n))  # define wavelet array
    wave = wave + 1j * wave  # make it complex
    # loop through scales and compute transform
    for a1 in range(J1 + 1):
        daughter, fourier_factor, coi, dofmin = wave_bases(
            mother, k, scale[a1], param)  # call wave_bases
        wave[a1, :] = np.fft.ifft(f * daughter)  # wavelet transform
        if a1 == 11:
            ondaleta = daughter
    # ondaleta = daughter
    period = np.array(period)
    period = period[:] * fourier_factor

    # cone-of-influence, differ for uneven len of timeseries:
    if (((n1) / 2.0).is_integer()) is True:
        # create mirrored array)
        mat = np.concatenate(
            (arange(1, int(n1 / 2)), arange(1, int(n1 / 2))[::-1]), axis=0)
        # insert zero at the begining of the array
        mat = np.insert(mat, 0, 0)
        mat = np.append(mat, 0)  # insert zero at the end of the array
    elif (((n1) / 2.0).is_integer()) is False:
        # create mirrored array
        mat = np.concatenate(
            (arange(1, int(n1 / 2) + 1), arange(1, int(n1 / 2))[::-1]), axis=0)
        # insert zero at the begining of the array
        mat = np.insert(mat, 0, 0)
        mat = np.append(mat, 0)  # insert zero at the end of the array
    coi = [coi * dt * m for m in mat]  # create coi matrix
    # problem with first and last entry in coi added next to lines because 
    # log2 of zero is not defined and cannot be plottet later:
    coi[0] = 0.1  # coi[0] is normally 0
    coi[len(coi)-1] = 0.1 # coi[last entry] is normally 0 too
    wave = wave[:, 0:n1]
    return ondaleta, wave, period, scale, coi, f


def wave_signif(Y, dt, scale1, sigtest, lag1, sig1v1, dof, mother, param):
    """CAUTION : default values"""
    import scipy
    from scipy import stats

    n1 = np.size(Y)
    J1 = len(scale1) - 1
    s0 = np.min(scale1)
    dj = np.log10(scale1[1] / scale1[0]) / np.log10(2)
    """CAUTION"""
    if (n1 == 1):
        variance = Y
    else:
        variance = np.var(Y)
    """CAUTION"""
    # sig1v1 = 0.95
    if (mother == 'Morlet'):
        # get the appropriate parameters [see table2]
        param = 6
        k0 = param
        fourier_factor = float(4 * math.pi) / (k0 + np.sqrt(2 + k0 * k0))
        empir = [2, -1, -1, -1]
        if(k0 == 6):
            empir[1:4] = [0.776, 2.32, 0.6]

    if(mother == 'DOG'):
        param = 2
        k0 = param
        m = param
        fourier_factor = float(2 * math.pi / (np.sqrt(m + 0.5)))
        empir = [1, -1, -1, -1]
        if(k0 == 2):
            empir[1:4] = [3.541, 1.43, 1.4]

    if (mother == 'PAUL'):
        param = 4
        m = param
        fourier_factor = float(4 * math.pi / (2 * m + 1))
        empir = [2., -1, -1, -1]
        if (m == 4):
            empir[1:4] = [1.132, 1.17, 1.5]

    period = [e * fourier_factor for e in scale1]
    dofmin = empir[0]  # Degrees of  freedom with no smoothing
    Cdelta = empir[1]  # reconstruction factor
    gamma_fac = empir[2]  # time-decorrelation factor
    dj0 = empir[3]  # scale-decorrelation factor
    freq = [dt / p for p in period]
    fft_theor = [((1 - lag1 * lag1) / (1 - 2 * lag1 *
                                       np.cos(f * 2 * math.pi) + lag1 * lag1))
                 for f in freq]
    fft_theor = [variance * ft for ft in fft_theor]
    signif = fft_theor
    if(dof == -1):
        dof = dofmin
    """CAUTION"""
    if(sigtest == 0):
        dof = dofmin
        chisquare = scipy.special.gammaincinv(dof / 2.0, sig1v1) * 2.0 / dof
        signif = [ft * chisquare for ft in fft_theor]
    elif (sigtest == 1):
        """CAUTION: if len(dof) ==1"""
        dof = np.array(dof)
        truncate = np.where(dof < 1)
        dof[truncate] = np.ones(np.size(truncate))
        for i in range(len(scale1)):
            dof[i] = (
                dofmin * np.sqrt(1 + pow((dof[i] * dt / gamma_fac / scale1[i]),
                                         2)))
        dof = np.array(dof)  # has to be an array to use np.where
        truncate = np.where(dof < dofmin)
        # minimum DOF is dofmin
        dof[truncate] = [dofmin * n for n in np.ones(len(truncate))]
        chisquare, signif = [], []
        for a1 in range(J1 + 1):
            chisquare.append(
                scipy.special.gammaincinv(dof[a1] / 2.0, sig1v1) * 2.0 /
                dof[a1])
            signif.append(fft_theor[a1] * chisquare[a1])
    """CAUTION : missing elif(sigtest ==2)"""
    return signif, fft_theor
   

#function used to normalise timeseries
def normalize(data):
    """
    NORMALIZE FUNCTION by - mean/sqrt(variance)
    """
    variance = np.var(data)
    data = (data - np.mean(data)) / (np.sqrt(variance))
    return data


#continuous wavelet transform function.
def cwt(data, dt, pad, dj, s0, j1, lag1, param, mother, name):
    """
    CONTINUOUS WAVELET TRANSFORM
    pad = 1         # pad the time series with zeroes (recommended)
    dj = 0.25       # this will do 4 sub-octaves per octave
    s0 = 2*dt       # this says start at a scale of 6 months
    j1 = 7/dj       # this says do 7 powers-of-two with dj sub-octaves each
    lag1 = 0.72     # lag-1 autocorrelation for red noise background
    param = 6
    mother = 'Morlet'
    """

    #from cwt.lib_wavelet import wavelet,wave_signif

    variance = np.var(data)
    n = len(data)
    # Wavelet transform
    ondaleta, wave, period, scale, coi, f = wavelet(
        data, dt, param, dj, s0, j1, mother)
    # wave = np.array(wave)
    power = (np.abs(wave) ** 2)
    # Significance levels: (variance=1 for the normalized SST)
    signif, fft_theor = wave_signif(
        1.0, dt, scale, 0, lag1, 0.95, -1, mother, param)
    ones = np.ones((len(signif), n))  # expand signif --> ones (J+1)x(n)
    sig95 = [s * ones[1] for s in signif]   # vstack signif concatenate
    sig95 = power / sig95  # where ratio > 1, power is significant
    # Global wavelet spectrum & significance levels:
    global_ws = variance * (np.sum(power.conj().transpose(), axis=0) / n)
    dof = [n - s for s in scale]
    """CAUTION - DEFAULT VALUES """
    global_signif, fft_theor = wave_signif(
        variance, dt, scale, 1, lag1, 0.95, dof, mother, param)
    # Daughter wavelet

    joint_wavelet = np.concatenate((np.fft.ifft(ondaleta)[int(np.ceil(
        n / 2.)):], np.fft.ifft(ondaleta)[int(np.ceil(n / 2.)):][::-1]), axis=0)
    imag_wavelet = np.concatenate((np.fft.ifft(ondaleta).imag[int(np.ceil(
        n / 2.)):], np.fft.ifft(ondaleta).imag[int(np.ceil(n / 2.)):][::-1]), axis=0)
    nw = np.size(joint_wavelet)  # daughter's number of points
    # admissibility condition
    mean_wavelet = mean(joint_wavelet.real)
    mean_wavelet = np.ones(nw) * mean_wavelet
    result = {'ondaleta': ondaleta, 'wave': wave, 'period': period,
              'scale': scale, 'coi': coi, 'power': power, 'sig95': sig95,
              'global_ws': global_ws, 'global_signif': global_signif,
              'joint_wavelet': joint_wavelet, 'imag_wavelet': imag_wavelet,
              'nw': nw, 'mean_wavelet': mean_wavelet, 'dj': dj, 'j1': j1,
              'dt': dt, 'fft': f, 'mother': mother, 'data': data, 'name': name, 'fft_theor': fft_theor}
    return result

#call like
# result = cwt(data_norm,0.25,1,0.25,2*0.25,7/0.25,0.72,6,'Morlet')


#fast fourier transform of a timeseries 
def fft(data):
    """FFT spectrum
    """
    n = len(data)
    X = np.fft.fft(data)
    sxx = ((X * np.conj(X)) / (n))
    f = -np.fft.fftfreq(n)[int(np.ceil(n / 2.)):]
    sxx = np.abs(sxx)
    sxx = sxx[int(np.ceil(n / 2.)):]
    return f, sxx

# ---------------------------
#           Ploting
# ---------------------------


def levels(result, dtmin):
    """
    Power levels
    """

    dtmax = result['power'].max()
    lev = []
    for i in range(int(log2(dtmax / dtmin))):
        dtmin = dtmin * 2
        lev.append(dtmin)
    return lev



from scipy.stats import chi2
def cross_wavelet(result1, result2, time):
    """ Computes the cross wavelet analysis.
    wave1 = result['wave'] time serie 1
            wave2 = result['wave'] time serie 2
    A normalized time and scale resolved measure for the relationship
    between two time series x1(t) and x2(t) is the wavelet coherency (WCO),
    which is defined as the amplitude of the WCS(wavelet cross spectrum)
    normalized to the two single WPS(wavelet power spectrum) (Maraun and
    Kurths,2004).
    WCOi(s)**2  = ((WPS12)*(WPS12)) / ((WPS1 * WPS2))**2
    _____________________________________________________________________
    Inputs:
    wave1 - wavelet transform of time series x1
            wave2 - wavelet transform of time series x2
    Outputs:
    cohere - wavelet coherency (WCO)
    Call function:
    cohere = cross_wavelet(wave1,wave2)
    """
    wave1 = result1['wave']
    wave2 = result2['wave']

    cross_power = np.abs(wave1 * wave2.conjugate())
    WPS12 = (wave1 * wave2.conjugate())
    WPS1  = (wave1 * wave1.conjugate())
    WPS2  = (wave2 * wave2.conjugate())
    coherence = ((WPS12)*(WPS12)) / ((WPS1 * WPS2))
    coherence = np.real(coherence)
    pot_cohere = np.around(coherence**2,3) # round numbers for digits to be in interval 0 a 1
    phase_angle = np.angle(WPS12)
    


    #significance for cross power, taken from methodology of Grinsted et al. 2004
    #significance depends on product of chai sqaured distributions and theoretical power spectra of each time series
    power1 = np.array(result1['fft_theor'])
    power2 = np.array(result2['fft_theor'])
    power_product = np.sqrt(power1*power2)
    Z = 3.999
    signif = (Z/2) * power_product 
    
    n = len(time)
    ones = np.ones((len(signif), n))  # expand signif --> ones (J+1)x(n)
    sig95 = [s * ones[1] for s in signif]   # vstack signif concatenate
    sig95 = cross_power / sig95  # where ratio > 1, power is significant
    
    
    
    return cross_power, pot_cohere, phase_angle, sig95



    
