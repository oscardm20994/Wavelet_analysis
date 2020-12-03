import iris
import numpy as np
from pylab import *
from wavelet_lib import *
from wavelet_plot import *
from numpy import genfromtxt
import collections
from scipy.stats import norm
from scipy import signal

##year list, 1000 years of data produced by a MetOffice Climate model
years = np.arange(2060,3061,1)

### read in Pacific Decadal Oscillation index fro climate model output
### which tracks sea surface temperatures over 
### the northern Pacific region. Load nc file using iris package
PDO = iris.load('PDO.nc')[0]
PDO = PDO.extract(iris.Constraint(month = ['Sep' 'Oct', 'Nov'])).data


##decadal smoothing of the series to
## examine long timescale variability
PDO_decadal = []
for i in range(len(years) - 10):
	PDO.append(np.mean(PDO[i:i+10]))


#function to construct wavelet object from timeseries
#lag1 autocorrelation included for the significance 
def get_result(data, lag):
	data_norm = normalize(data)
	var = np.var(data_norm)
	alpha = np.corrcoef(data_norm[0:-1*lag], data_norm[lag:])[0,1]
	print alpha
	result = cwt(data_norm, 1, 1, 0.25, 0.5, 7/0.5, alpha, 6, mother='Morlet', name='x')
	return result, data_norm


#call get_result function above for PDO example
result_PDO, data_norm_PDO = get_result(PDO, 1)



#function to call plotter with appropriate time array
def main(data, data_norm, result, title, axis_name, filename, start_year, data_type):
	z = np.linspace(0,len(data_norm), len(data_norm))
	wavelet_plot(title, z, data, 0.03125, result, axis_name, filename, start_year, data_type)
	return

#plot PDO wavelet spectrum
main(PDO, data_norm_PDO, result_PDO, 'PDO', 'K', 'PDO_wavelet.png', 2060, 'SSW')


#exit gracefully
sys.exit()



	


	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	


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


