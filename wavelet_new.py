import waipy
import iris
import numpy as np
from pylab import *
from waipy import *
from numpy import genfromtxt
import collections
from scipy.stats import norm
from scipy import signal
###stick with UKESM1 data for now, 1000 years
#
def correct_nov_year(SSW_year, SSW_month):
	SSW_year2 = SSW_year[:]
	for i in range(len(SSW_year)):
		if SSW_month[i] == 'Nov':
			SSW_year2[i] += 1

	return SSW_year2

#
def get_coord_from_list(List, coord_name):
	coord = []
	
	for i in range(len(List)):
		coord.append(List[i].coord(coord_name).points[0])
	
	return coord





##northern annular mode
NAM_UKESM = iris.load('/gws/nopw/j04/aopp/odimdore/diagnostics/annular_modes/UKESMNAM*')[0].data

## iris list of SSW events
SSW_UKESM = iris.load('../equator_diagnostics/SSW_UKESM_cube.nc')[0]
SSW_UKESM_DJFM = SSW_UKESM.extract(iris.Constraint(month = ['Dec','Jan','Feb','Mar']))
## get years and months out of list
SSW_UKESM_year = np.sort(SSW_UKESM.coord('year').points)#get_coord_from_list(SSW_UKESM_list, 'season_year')
SSW_UKESM_yearDJFM = np.sort(SSW_UKESM_DJFM.coord('year').points)#get_coord_from_list(SSW_UKESM_list, 'season_year')
SSW_month_UKESM = get_coord_from_list(SSW_UKESM_list, 'month')
SSW_UKESM_yearDJF = get_coord_from_list(SSW_UKESM_listDJF, 'season_year')
#correct November warming years (add 1 to the year)
SSW_UKESM_year2 = correct_nov_year(np.array(SSW_UKESM_year), np.array(SSW_month_UKESM))



#count SSWs in each year
count = collections.Counter(SSW_UKESM_year)
countDJFM = collections.Counter(SSW_UKESM_yearDJFM)


##year list
UKESM_years = np.arange(2060,3061,1)

##make SSW time series
SSW_UKESM,SSW_UKESM_DJFM,SSW_UKESM_DJFM_series = [],[],[]
for i in UKESM_years:
	SSW_UKESM.append(count[i])
	SSW_UKESM_DJFM_series.append(countDJF[i])
SSW_UKESM[SSW_UKESM == 3] = 2



UKESM_DJFM_decadal = []
## decadal smoothing if necessary
for i in range(len(UKESM_years) - 10):
	UKESM_DJFM_decadal.append(mean(SSW_UKESM_DJFM_series[i:i+10]))

	
### load QBO index data
QBO_UKESM = iris.load('QBO_index_20.nc')[0]
QBO_SON = QBO_UKESM.extract(iris.Constraint(clim_season = 'son'))
QBO_SON = QBO_SON.aggregated_by('year', iris.analysis.MEAN).data

##ENSO index data
fname = '/gws/nopw/j04/aopp/strat_clim_data/CMIP6/UKESM1/picontrol/SST_UKESM/UKESM_2060-3060_EN3.4index.nc'
ENSO_UK = iris.load(fname)[0]
ENSO_UK = add_times(ENSO_UK, 'time')
ENSO_SON = ENSO_UK.extract(iris.Constraint(clim_season = 'son'))
ENSO_SON = ENSO_SON.aggregated_by('year', iris.analysis.MEAN).data

### read in AMV index for 
AMV_DJF = iris.load('AMO_UK_DJF.nc')[0].data
AMV = iris.load('AMO_UK.nc')[0].data

##PDO data
PDO_UK_all = iris.load('PDO_UK_NDJF.nc')[0]
PDO_UK_all_SON = iris.load('PDO_UK_SON.nc')[0]

PDO_UK_NDJF = PDO_UK_all.extract(iris.Constraint(month = ['Nov','Dec','Jan','Feb']))
PDO_UK_N = PDO_UK_all.extract(iris.Constraint(month = 'Nov'))
PDO_UK_N = PDO_UK_N.aggregated_by('year', iris.analysis.MEAN).data
PDO_UK = PDO_UK_NDJF.aggregated_by('season_year', iris.analysis.MEAN)[0:-1].data
PDO_UK_SON = PDO_UK_all_SON.extract(iris.Constraint(month = ['Sep' 'Oct', 'Nov']))
PDO_UK_SON = PDO_UK_SON.aggregated_by('year', iris.analysis.MEAN).data

HT_10 = iris.load('HT_series_UKESM_10yrs_correlations.nc')[0].data


##decadal smoothing
PDO_UK_SON_decadal, AMV_UK_DJF_decadal, ENSO_SON_decadal, QBO_SON_decadal = [], [], [], []
for i in range(len(UKESM_years) - 10):
	PDO_UK_SON_decadal.append(np.mean(PDO_UK_SON[i:i+10]))
	AMV_UK_DJF_decadal.append(np.mean(AMV_DJF[i:i+10]))
	ENSO_SON_decadal.append(np.mean(ENSO_SON[i:i+10]))
	QBO_SON_decadal.append(np.mean(QBO_SON[i:i+10]))



##O3 data
O3_ASO_UKESM = iris.load('../wavelet_analysis/ASON_O3_UKESM.nc')[0].data


def get_result(data, lag):
	data_norm = normalize(data)
	var = np.var(data_norm)
	alpha = np.corrcoef(data_norm[0:-1*lag], data_norm[lag:])[0,1]
	print alpha
	result = cwt(data_norm, 1, 1, 0.25, 0.5, 7/0.5, alpha, 6, mother='Morlet', name='x')
	return result, data_norm

#UKESM SSWs per year data
result_SSW_UKESM, data_norm_UKESM = get_result(SSW_UKESM, 1)
result_SSW_UKESM_DJFM, data_norm_UKESM_DJFM = get_result(SSW_UKESM_DJFM_series, 1)
result_SSW_UKESM_DJFM_decadal, data_norm_UKESM_DJFM_decadal = get_result(UKESM_DJFM_decadal, 3)
result_SSW_UKESM_DJF, data_norm_UKESM_DJF = get_result(SSW_UKESM_DJF, 1)

#other data: AMV and QBO index, to add ENSO index later
result_QBO, data_norm_QBO = get_result(QBO_SON, 1)
result_QBO_decadal, data_norm_QBO_decadal  = get_result(QBO_SON_decadal, 1)
result_ENSO_SON, data_norm_ENSO = get_result(ENSO_SON, 1)
result_ENSO_SON_decadal, data_norm_ENSO_decadal = get_result(ENSO_SON_decadal, 1)
result_AMV_DJF, data_norm_AMV = get_result(AMV_DJF, 1)
result_AMV_DJF_decadal, data_norm_AMV_decadal = get_result(AMV_UK_DJF_decadal, 4)
result_AMV, data_norm_AMV = get_result(AMV, 1)
result_PDO_NDJF, data_norm_PDO = get_result(PDO_UK, 1)
result_PDO_N, data_norm_PDO_N = get_result(PDO_UK_N, 1)
result_PDO_SON, data_norm_PDO_SON = get_result(PDO_UK_SON, 1)
result_PDO_SON_decadal, data_norm_PDO_SON_decadal = get_result(PDO_UK_SON_decadal, 3)
result_O3_ASO, data_norm_O3_ASO = get_result(O3_ASO_UKESM, 1)
result_HT_10, data_norm_HT_10 = get_result(HT_10, 4)

#function to create set of surrogates by shuffling SSW time series and calculating power spectra
#outputs points on real power spectrum distinct to the set to a given level
def bootstrap(SSW, n, time):
	real_result = get_result(SSW, 1)[0]
	real_spectra = real_result['power']
	periods = real_result['period']
	nscale = len(real_spectra[:,0])
	ntime = len(real_spectra[0,:])

	set_of_spectra = np.zeros([n,nscale,ntime])

	for i in range(n):
		print i
		shuffled = np.copy(SSW)
		np.random.shuffle(shuffled)
		print shuffled
		shuffled_spectra = get_result(shuffled, 1)[0]['power']
		set_of_spectra[i,:,:] = shuffled_spectra
	
	#find mean and standard deviation of synthetic spectra set
	mean = np.mean(set_of_spectra, axis = 0) 
	stdev = np.std(set_of_spectra, axis = 0)
	
	Z_score = (real_spectra - mean)/stdev
	
	p = norm.cdf(Z_score)
	
	t_sig = time[np.where(p > 0.95)[1]]
	s_sig = periods[np.where(p > 0.95)[0]]
	
	print p	
	#relative_power[relative_power < 1.5] = float('nan') 
	#relative_power[relative_power > 1.5] = 1 

	#t_sig = np.empty(0)
	#s_sig = np.empty(0)
	#for t in range(ntime):
	#	for s in range(nscale):
	#		if real_spectra[s,t] >= mean[s,t] + 1.96*stdev[s,t]:
	#			t_sig = np.append(t_sig, time[t])
	#			s_sig = np.append(s_sig, periods[s])
	
	return p, t_sig, s_sig
	
P_UKESM_DJF, t_sig, s_sig = bootstrap(SSW_UKESM_DJF, 1000, np.arange(1001))
	
				


def wavelet_plot(var, time, data, dtmin, result, axis_name, filename, start_year, data_type):#, t_sig, s_sig):
    """
    PLOT WAVELET TRANSFORM
    var = title name from data
    time  = vector get in load function
    data  = from normalize function
    dtmin = minimum resolution :1 octave
    result = dict from cwt function
    """

    from numpy import log2
    import numpy as np
    #import wavetest
    import matplotlib
    import matplotlib.gridspec as gridspec

    # frequency limit
    # print result['period']
    # lim = np.where(result['period'] == result['period'][-1]/2)[0][0]
    """Plot time series """

    fig = plt.figure(figsize=(12, 8), dpi=80)

    gs1 = gridspec.GridSpec(4, 3)
    gs1.update(left=0.05, right=0.7, wspace=0.5, hspace=0)
    ax1 = plt.subplot(gs1[0, :])
    ax1.xaxis.set_visible(False)
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax2 = plt.subplot(gs1[1:4, :])#, axisbg='#C0C0C0')

    gs2 = gridspec.GridSpec(4, 1)
    gs2.update(left=0.7, right=0.98, hspace=0)
    ax5 = plt.subplot(gs2[1:4, 0], sharey=ax2)
    plt.setp(ax5.get_yticklabels(), visible=False)

    gs3 = gridspec.GridSpec(6, 1)
    gs3.update(left=0.77, top=0.86, right=0.98, hspace=0.6, wspace=0.01)
    ax3 = plt.subplot(gs3[0, 0])

    # ----------------------------------------------------------------------------------------------------------------#
    if data_type == 'SSW':
    	ax1.bar(time, data)
        ax1.set_yticks([0,1,2])
	ax1.set_yticklabels([0,1,2])

    else:
   	ax1.plot(time, data)
    ax1.axis('tight')
    ax1.set_ylabel(axis_name)
    ax1.set_title('%s' % var)
    #ax1.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    #ax1.grid(True)
    ax1.xaxis.set_visible(False)
    ax1.set_xlim(time[0], time[-1]+1)
    # ----------------------------------------------------------------------------------------------------------------#
    ax3.plot(arange(-result['nw'] / 2, result['nw'] / 2),result['joint_wavelet'], 'k', label='Real part')
    ax3.plot(arange(-result['nw'] / 2, result['nw'] / 2),result['imag_wavelet'], '--k', label='Imag part')
    ax3.plot(arange(-result['nw'] / 2, result['nw'] / 2),result['mean_wavelet'], 'g', label='Mean')
    # ax3.axis('tight')
    ax3.set_xlim(-40, 40)
    # ax3.set_ylim(-0.3,0.3)
    # ax3.set_ylim([np.min(result['joint_wavelet']),np.max(result['joint_wavelet'])])
    ax3.set_xlabel('Time (days)', fontsize=10)
    ax3.set_ylabel('Amplitude', fontsize=10)
    ax3.set_title('$\psi$ (t/s) {0} in time domain'.format(result['mother']))
    # ----------------------------------------------------------------------------------------------------------------#
    # ax4.plot(result['ondaleta'],'k')
    # ax4.set_xlabel('Frequency', fontsize=10)
    # ax4.set_ylabel('Amplitude', fontsize=10)
    # ax4.set_title('$\psi^-$  Frequency domain', fontsize=13)
    # ----------------------------------------------------------------------------------------------------------------#
    """ Contour plot wavelet power spectrum """
    #lev = np.arange(15)#np.min(result['power']), np.max(result['power']))
    pc = ax2.contourf(time, np.log2(result['period']), result['power'], cmap = 'Reds')
    print 2060 + time[475:685]
    #ax2.scatter(t_sig, np.log2(s_sig), s = 0.8, marker = '*',color = 'black')
    #ax2.plot([475,475], np.log2(np.array([55,135])), linestyle = '--', color = 'blue')
    #ax2.plot([685,685], np.log2(np.array([55,135])), linestyle = '--', color = 'blue')
    #ax2.plot([475,685], np.log2(np.array([55,55])), linestyle = '--', color = 'blue')
    #ax2.plot([475,685], np.log2(np.array([135,135])), linestyle = '--', color = 'blue')
    # 95% significance contour, levels at -99 (fake) and 1 (95% signif)
    pc2 = ax2.contour(time, np.log2(result['period']), result['sig95'],[-99, 1], linewidths=2)
    #ax2.plot(time, np.log2(result['coi']), 'k')
    # cone-of-influence , anything "below"is dubious
    ax2.fill_between(time, np.log2(result['coi']), int(np.log2(result['period'][-1]) + 1), alpha=0.5, hatch='/')
    position = fig.add_axes([0.07, 0.05, 0.6, 0.01])
    cbar = plt.colorbar(pc, cax=position, orientation='horizontal')
    cbar.set_label('Power')
    yt = range(int(np.log2(result['period'][1])), int(np.log2(result['period'][-1]) + 1))  # create the vector of periods
    
    
    
    #y_labels = np.array([1,2,4,8,16,32,64,128,256])
    #Yticks = np.log2(y_labels)
    if data_type == 'NAM':   
    	Yticks = [float(math.pow(2, p)) for p in yt]  # make 2^periods
    else:
    	Yticks = [float(math.pow(2, p)) for p in yt]
    Yticks = [int(i) for i in Yticks]
    print Yticks
    
    ax2.set_yticks(yt[1:])
    ax2.set_yticklabels(Yticks[1:])
    ax2.set_xticks(np.arange(0, 1100, 100))
    ax2.set_xticklabels(np.arange(0, 1100, 100) + start_year)
    ax2.set_ylim(ymin=(0, np.log2(np.max(result['period']))))
    ax2.set_xlim(time[0],time[-1:])
    ax2.set_ylim(ax2.get_ylim()[::-1])
    ax2.set_xlabel('Year', fontsize=12)
    ax2.set_ylabel('Period (years)', fontsize=12)
    # ----------------------------------------------------------------------------------------------------------------#
    """ Plot global wavelet spectrum """
    #f, sxx = fft(data)
    #ax5.plot(sxx, np.log2(1 / f * result['dt']), 'gray', label='Fourier spectrum')
    ax5.plot(result['global_ws'], np.log2(result['period']), 'b', label='Wavelet spectrum')
    ax5.plot(result['global_signif'], np.log2(result['period']), 'r--', label='95% confidence spectrum')
    ax5.legend(loc=0)
    ax5.set_xlim(0, 1.05 * np.max(result['global_ws']))
    ax5.set_xticks(np.arange(1, int(1.05 * np.max(result['global_ws'])), 5))
    ax5.set_xlabel('Power', fontsize=10)
    ax5.set_title('Global Wavelet Spectrum', fontsize=12)
    # save fig
    plt.show()
    fig.savefig(filename, dpi=180)




def main(data, data_norm, result, title, axis_name, filename, start_year, data_type):#, t_sig, s_sig):
	z = np.linspace(0,len(data_norm), len(data_norm))
	wavelet_plot(title, z, data, 0.03125, result, axis_name, filename, start_year, data_type)#, t_sig, s_sig)
	return result
	
result_UKESM_SSW = main(SSW_UKESM, data_norm_UKESM, result_SSW_UKESM, 'NDJFM SSWs (UKESM)', 'SSWs', 'SSW_wavelet_UKESM_NDJFM.png', 2060, 'SSW')
result_UKESM_SSW = main(SSW_UKESM_DJFM_series, data_norm_UKESM_DJFM, result_SSW_UKESM_DJFM, 'DJFM SSWs (UKESM)', 'SSWs', 'SSW_wavelet_UKESM_DJFM.png', 2060, 'SSW')
result_UKESM_SSW = main(SSW_UKESM_DJF, data_norm_UKESM_DJF, result_SSW_UKESM_DJF, 'DJF SSWs (UKESM)', 'SSWs', 'SSW_wavelet_UKESM_DJF.png', 2060, 'SSW')
result_UKESM_SSW = main(UKESM_DJFM_decadal, data_norm_UKESM_DJFM_decadal, result_SSW_UKESM_DJFM_decadal, 'DJFM SSW rate (UKESM)', 'SSWs/season', 'SSW_wavelet_UKESM_DJFM_decadal.png', 2060, 'other')
result_UKESM_QBO = main(QBO_SON, data_norm_QBO, result_QBO, 'SON mean QBO 20hPa (UKESM)', 'QBO', 'QBO20hPa_wavelet_UKESM.png', 2060, 'NAM')
result_UKESM_ENSO = main(ENSO_SON, data_norm_ENSO, result_ENSO_SON, 'SON mean ENSO index (UKESM)', 'ENSO index', 'ENSO_SON_wavelet_UKESM.png', 2060, 'NAM')
result_UKESM_ENSO = main(ENSO_SON_decadal, data_norm_ENSO_decadal, result_ENSO_SON_decadal, 'SON mean ENSO index (UKESM)', 'ENSO index', 'ENSO_SON_decadal_wavelet_UKESM.png', 2060, 'NAM')
result_AMV_DJF = main(AMV_DJF, data_norm_AMV, result_AMV_DJF, 'DJF AMV (UKESM1)', 'AMV', 'AMV_wavelet_UKESM_DJF.png', 2060, 'NAM')
result_AMV_DJF = main(AMV_UK_DJF_decadal, data_norm_AMV_decadal, result_AMV_DJF_decadal, 'decadal mean DJF AMV (UKESM1)', 'AMV', 'AMV_wavelet_UKESM_DJF_decadal.png', 2060, 'NAM')
result_AMV_DJF = main(AMV, data_norm_AMV, result_AMV, 'AMV (UKESM1)', 'AMV', 'AMV_wavelet_UKESM.png', 2060, 'NAM')
result_PDO_NDJF = main(PDO_UK, data_norm_PDO, result_PDO, 'NDJF PDO (UKESM1)', 'PDO index', 'PDO_NDJFM_wavelet_UKESM.png', 2060, 'NAM')
result_PDO_N = main(PDO_UK_N, data_norm_PDO_N, result_PDO_N, 'November PDO (UKESM1)', 'PDO index', 'PDO_N_wavelet_UKESM.png', 2060, 'NAM')
result_PDO_SON = main(PDO_UK_SON, data_norm_PDO_SON, result_PDO_SON, 'SON PDO (UKESM1)', 'PDO index', 'PDO_SON_wavelet_UKESM.png', 2060, 'NAM')
result_PDO_SON = main(PDO_UK_SON_decadal, data_norm_PDO_SON_decadal, result_PDO_SON_decadal, 'SON decadal mean PDO (UKESM)', 'PDO index', 'PDO_SON_decadal_wavelet_UKESM.png', 2060, 'NAM')
result_O3_SON = main(O3_ASO_UKESM, data_norm_O3_ASO, result_O3_ASO, 'ASO 10-20hPa Ozone mass mixing ratio (UKESM1)', 'O$_3$ mass mixing ratio', 'O3_ASO_10-20hPa_wavelet_UKESM.png', 2060, 'NAM')
result_HT_10 = main(HT_10, data_norm_HT_10, result_HT_10, 'decadal mean Holton-Tan strenght (UKESM)', 'correlation', 'HT_10_wavelet_UKESM.png', 2060, 'NAM')




	


	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	


def plot_cross(var, cross_power, phase_angle, time, result, result1, sig95, figname):
    """ PLOT CROSS POWER
    cross_power = from cross_wavelet function
    coherence   = from cross_wavelet function
    """
    import matplotlib.gridspec as gridspec
    fig = plt.figure(figsize=(12, 11), dpi=80)
    # set plot grid
    gs1 = gridspec.GridSpec(7, 3)
    gs1.update(left=0.05)
    ax1 = plt.subplot(gs1[0, :])
    ax1 = pylab.gca()
    ax1.xaxis.set_visible(False)
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax2 = plt.subplot(gs1[1:4, :])# axisbg='#C0C0C0')
    # plot timeseries
    ax1.plot(time, UKESM_DJFM_decadal, color = 'b')
    ax1.set_ylabel('SSWs/winter', fontsize=13)
    ax1.axis('tight')
    ax3 = ax1.twinx()
    ax3.plot(time, result1['data'], color='r')
    ax3.set_ylabel('AMV index')
    ax3.axis('tight')
    ax1.set_title('%s' % var, fontsize=15)
    ax1.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    #ax1.grid(True)
    ax1.xaxis.set_visible(False)
    ax1.set_yticks([0,0.5,1])
    ax1.set_xlim(2060,3060)

    
    phs_dt = round(len(time) / 20)
    tidx = np.arange(np.max(np.floor(phs_dt / 2)), len(time), phs_dt)
    tidx = [int(i) for i in tidx]
    tidx = np.array(tidx)
    phs_dp = round(len(result['period']) / 20)
    pidx = np.arange(
        np.max(np.floor(phs_dp / 2)), len(result['period']), phs_dp)
    pidx = [int(i) for i in pidx]
    pidx = np.array(pidx)
    X, Y = meshgrid(
        time.astype(np.int64)[tidx], np.log2(result['period'][pidx]))

    #Arrows indicate in phase when pointing to the right and out of phase when pointing left.
    phase_angle1 = phase_angle[:, tidx]
    phase_angle1 = phase_angle1[pidx, :]
    cA = np.exp(1j * phase_angle1)
    U = np.real(cA)
    V = np.imag(cA)
    ax4 = ax2.twiny()
    ax4.xaxis.set_visible(False)
    # ax4.set_xlim(0.9,4.4)
    CS = ax2.contourf(time, np.log2(result['period']), np.abs(cross_power), cmap = 'Reds')
    pc2 = ax2.contour(time, np.log2(result['period']), sig95,[-99, 1], linewidths=2)
    # cone-of-influence , anything "below"is dubious
    ax2.plot(time, np.log2(result['coi']), 'k')
    ax2.fill_between(time, np.log2(result['coi']), int(
        np.log2(result['period'][-1]) + 1), alpha=0.5, hatch='/')
    position = fig.add_axes([0.15, 0.40, 0.6, 0.01])
    # ,norm=normal)#,shrink=0.5,pad=0.08)
    cbar = plt.colorbar(CS, cax=position, orientation='horizontal')
    cbar.set_label('Power')
    Q = ax4.quiver(X.astype(np.int64), Y, U, V, linewidth=0.1)
    ax4.axis('tight')
    yt = range(int(np.log2(result['period'][0])), int(
        np.log2(result['period'][-1]) + 1))  # create the vector of periods
    Yticks = [float(math.pow(2, p)) for p in yt]  # make 2^periods
    Yticks = [int(i) for i in Yticks]
    ax2.set_yticks(yt)
    ax2.set_yticklabels(Yticks)
    ax2.set_ylim(ymin=(np.log2(result['period'][0])), ymax=(
        np.log2(result['period'][-1])))
    ax2.invert_yaxis()
    ax2.set_xlabel('Year', fontsize=13)
    ax2.set_ylabel('Period', fontsize=13)
    ax2.axhline(y=10.5, xmin=0, xmax=1, linewidth=2, color='k')
    ax2.axhline(y=13.3, xmin=0, xmax=1, linewidth=2, color='k')
	
    #ax4 = plt.subplot(gs1[5:, :])
    #plt.title('PDO and SSWs per season fourier filtered keeping periods 60-120 years', fontsize = 15)
    #ax4.plot(np.arange(1000) + 2060, low_PDO_SON, label = 'PDO')
    #ax4.plot(np.arange(1000) + 2060, low_SSW, label = 'SSWs')
    #ax4.plot([2500,2500],[-0.5,0.5], linestyle = '--', color = 'black')
    #ax4.plot([2750,2750],[-0.5,0.5], linestyle = '--', color = 'black')
    #ax4.set_ylim(-0.5,0.5)
    #ax4.set_xlim(2060,3060)
    #ax4.legend()
    fig.tight_layout()
    plt.savefig('%s'%figname, dpi=200)
    return

#Coherence plots
time = 2060 + np.arange(1001-10)
cross_power, coherence, phase_angle, sig95 = cross_wavelet(result_SSW_UKESM_DJFM_decadal, result_AMV_DJF_decadal, time)
plot_cross('Cross Power Spectrum, SSWs per winter vs DJF AMV (UKESM)', cross_power, phase_angle, time, result_SSW_UKESM_DJFM_decadal, result_AMV_DJF_decadal, sig95, 'cross_coherence_decadal_AMV_DJF_UKESM.png')



# define fourier filter function to isolate signals in data with high cross power
def fourier_filter(GPH, lowf, highf):
#time step in days
	time_step = 1
	fft = np.fft.rfft(GPH)
	print fft.shape
	freqs = np.fft.fftfreq(GPH.size, time_step)[0:len(fft)]
	cut_f_signal = fft.copy()
	cut_f_signal[(freqs > highf)] = 0
	cut_f_signal[(freqs < lowf)] = 0

	inverse = np.fft.irfft(cut_f_signal)
	return inverse


low_PDO = fourier_filter(PDO_UK_SON_decadal, lowf, highf)
low_PDO_N = fourier_filter(PDO_UK_N, lowf, highf)
low_PDO_SON = fourier_filter(np.array(PDO_UK_SON_decadal), lowf, highf)
low_AMV_DJF = fourier_filter(np.array(AMV_UK_DJF_decadal), lowf, highf)
low_SSW = fourier_filter(np.array(UKESM_DJFM_decadal), lowf, highf)
#low_NAO = fourier_filter(NAO.data, lowf, highf)

#define time array
time = np.arange(len(shift)) + 2060 - 10

N = 1001-10 # number of smaples
T = 1.0 # sample spacing


lowf = 1./80.
highf = 1./65.
low_PDO_SON = fourier_filter(np.array(PDO_UK_SON_decadal), lowf, highf)
low_SSW = fourier_filter(np.array(UKESM_DJFM_decadal), lowf, highf)
low_AMV_DJF = fourier_filter(np.array(AMV_UK_DJF_decadal), lowf, highf)

#y1 = low_PDO_SON
y1 = low_AMV_DJF
y2 = low_SSW


al1 = np.angle(hilbert(y1),deg=False)
al2 = np.angle(hilbert(y2),deg=False)
fig,ax = plt.subplots(3,1,figsize=(11,7),sharex=True, dpi = 90)
ax[0].plot(y1,color='r',label='AMV')
ax[0].plot(y2,color='b',label='SSWs')
ax[0].legend(bbox_to_anchor=(0., 1.02, 1., .102),ncol=2)
ax[0].set(xlim=[0,N], title='Filtered SSW and AMV timeseries')
ax[1].plot(al1,color='r')
ax[1].plot(al2,color='b')
ax[1].set(title='Angle at each Timepoint')
ax[1].set(yticks = [-3.1,0,3.1], yticklabels = ['$-\pi$','0','$\pi$']) 
phase_synchrony = 1-np.sin(np.abs(al1-al2)/2)
ax[2].plot(phase_synchrony)
ax[2].set(xlim=[0,N],title='Instantaneous Phase Synchrony',xlabel='Year Number',ylabel='Phase Synchrony')
plt.tight_layout()
fig.savefig('instant_phase_shift_SSWs_AMV_filtered_65-100.png', dpi = 200)
plt.show()

cor = [np.corrcoef(phase_synchrony, cross_power[i,0:-1])[0,1] for i in range(44)]




import matplotlib.gridspec as gridspec
fig = plt.figure(figsize=(15, 10))#, dpi = 400)
gs1 = gridspec.GridSpec(2,1)
ax4 = plt.subplot(gs1[0,0])
ax4.set_title('Fourier Filetered November PDO and SSWs per Season (64-120 years retained)', fontsize = 18)
#ax4.scatter(max_time_SSW, np.zeros(len(max_time_SSW)), color = 'orange')
#ax4.scatter(max_time_PDO, np.zeros(len(max_time_PDO)), color = 'blue') 

ax4.plot(np.arange(1000-10) + 2060, low_PDO_SON, label = 'PDO')
ax4.plot(np.arange(1000-10) + 2060, low_SSW, label = 'SSWs')
ax4.plot([2500,2500],[-0.5,0.5], linestyle = '--', color = 'black')
ax4.plot([2750,2750],[-0.5,0.5], linestyle = '--', color = 'black')
ax4.set_ylim(-0.5,0.5)
ax4.set_xlim(2060,3060)
plt.legend()

ax3 = plt.subplot(gs1[1])
ax3.set_title('Relative Phase Shift (' + str(window) + 'year window)' , fontsize = 18) 
ax3.plot(np.arange(len(shift_SON)) + 2060, shift_SON)
ax3.plot([2500,2500],[np.min(shift),np.max(shift)], linestyle = '--', color = 'black')
ax3.plot([2750,2750],[np.min(shift),np.max(shift)], linestyle = '--', color = 'black')
#ax3.set_ylim(-1.5,1.5)
ax3.set_xlim(2060,3060)
ax3.set_ylabel('phase difference (years)', fontsize = 18)
plt.legend()

#ax3 = plt.subplot(gs1[2])
#ax3.set_title('Time Rate of Change in Phase Shift') 
#ax3.plot(time, ROC_shift)
#ax3.plot([2500,2500],[np.min(ROC_shift),np.max(ROC_shift)], linestyle = '--', color = 'black')
#ax3.plot([2750,2750],[np.min(ROC_shift),np.max(ROC_shift)], linestyle = '--', color = 'black')
##ax3.set_ylim(-1.5,1.5)
#ax3.set_xlim(2060,3060)
#ax3.set_ylabel('ROC of phase difference (year/year)')
plt.tight_layout()
plt.show()
fig.savefig('pass_filtered_SON_PDO_vs_SSWs_UKESM.png', dpi = 300)


