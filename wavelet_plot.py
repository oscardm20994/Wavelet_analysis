import cartopy.crs as ccrs
import matplotlib.gridspec as gridspec
import iris.plot as iplt
import mpl_toolkits
import iris.analysis.stats
from iris.experimental.equalise_cubes import equalise_attributes
from mpl_toolkits.basemap import Basemap
import scipy.stats
import matplotlib


#script containing functions to plot wavelet power spectra and cross power spectra. 

def wavelet_plot(var, time, data, dtmin, result, axis_name, filename, start_year, data_type)
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
    gs1.update(left=0.05, top = 0.9, right=0.7, wspace=0.3, hspace=0)
    ax1 = plt.subplot(gs1[0, :])
    ax1.xaxis.set_visible(False)
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax2 = plt.subplot(gs1[1:4, :])#, axisbg='#C0C0C0')

    gs2 = gridspec.GridSpec(4, 1)
    gs2.update(left=0.7, right=0.95, hspace=0)
    ax5 = plt.subplot(gs2[1:4, 0], sharey=ax2)
    plt.setp(ax5.get_yticklabels(), visible=True)

    gs3 = gridspec.GridSpec(6, 1)
    gs3.update(left=0.77, top=0.88, right=0.98, hspace=0.3, wspace=0.01)
    ax3 = plt.subplot(gs3[0, 0])

    # ----------------------------------------------------------------------------------------------------------------#
    if data_type == 'SSW':
        ax1.bar(time, data)
        ax1.set_yticks([0,1,2])
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
    pc = ax2.contourf(time, np.log2(result['period']), result['power']*1.2, cmap = 'Reds')
    pc2 = ax2.contour(time, np.log2(result['period']), result['sig95']*1.2,[-99, 1], linewidths=2)
    #pc3 = ax2.contour(time, np.log2(result['period']), result_SSW_UKESM_DJFM_decadal['sig95'],[-99, 1], linewidths=2, colors = 'greenyellow')

    extra_sig = result_SSW_UKESM_DJFM_decadal['sig95']
    #ax2.contour(time, np.log2(result['period']), extra_sig,[-99, 1], linewidths=2, colors = 'greenyellow')
    #ax2.plot(time, np.log2(result['coi']), 'k')
    # cone-of-influence , anything "below"is dubious
    ax2.fill_between(time, np.log2(result['coi']), int(np.log2(result['period'][-1]) + 1), alpha=0.5, hatch='/')
    position = fig.add_axes([0.07, 0.05, 0.6, 0.01])
    cbar = plt.colorbar(pc, cax=position, orientation='horizontal')
    cbar.set_label('Power')
    yt = (np.log2(result['period']))
    yt = range(int(np.log2(result['period'][0])), int(
        np.log2(result['period'][-1]) + 1))# create the vector of periods
    Yticks = [float(math.pow(2, p)) for p in yt]# make 2^periods
    Yticks = [int(i) for i in Yticks]
    ax2.set_yticks(yt)
    ax2.set_yticklabels(Yticks)
    ax2.set_ylim(ymin=(np.log2(result['period'][i1])), ymax=(
        np.log2(result['period'][i2])))
    ax2.invert_yaxis()
    ax2.set_xlabel('Year Number', fontsize=13)
    ax2.set_ylabel('Period (years)', fontsize=13)
    # ----------------------------------------------------------------------------------------------------------------#
    """ Plot global wavelet spectrum """
    #f, sxx = fft(data)
    #ax5.plot(sxx, np.log2(1 / f * result['dt']), 'gray', label='Fourier spectrum')
    ax5.plot(result['global_ws'], np.log2(result['period']), 'b', label='Wavelet spectrum')
    ax5.plot(result['global_signif'], np.log2(result['period']), 'r--', label='95% confidence spectrum')
    ax5.legend(loc=0)
    #ax5.set_yticks(yt[1:])
    ax5.yaxis.tick_right()
    ax5.yaxis.set_label_position("right")

    ax5.set_xlim(0, 1.05 * np.max(result['global_ws']))
    ax5.set_xticks(np.arange(1, int(1.05 * np.max(result['global_ws'])), 1))
    ax5.set_xlabel('Power', fontsize=10)
    ax5.set_ylabel('Period (years)', fontsize=10)

    ax5.set_title('Global Wavelet Spectrum', fontsize=12)
    #ax5.set_xlim(0,11)
    # save fig
    plt.tight_layout()
    fig.subplots_adjust(left = 0.5)
    plt.show()
    fig.savefig(filename, dpi=180)
    return
    
def plot_cross(var, cross_power, phase_angle, time, result, result1, data1, data2, sig95, figname, extra_sig = []):
    """ PLOT CROSS POWER
    cross_power = from cross_wavelet function
    coherence   = from cross_wavelet function
    """
    import matplotlib.gridspec as gridspec
    fig = plt.figure(figsize=(10, 6), dpi = 300)
    # set plot grid
    gs1 = gridspec.GridSpec(7, 1)
    ax1 = plt.subplot(gs1[0:2, 0])
    ax1 = pylab.gca()
    ax1.xaxis.set_visible(False)
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax2 = plt.subplot(gs1[2:, 0])# axisbg='#C0C0C0')
    # plot timeseries
    ax1.plot(time, data2, color = 'b')
    ax1.set_ylabel('ENSO3.4 (K)', fontsize=11)
    ax1.axis('tight')
    ax3 = ax1.twinx()
    ax3.plot(time, data1, color='r', alpha = 0.6)
    ax3.set_ylabel('amplitude ($ms^{-1}$)', fontsize=11)
    #ax3.set_ylabel()
    ax3.axis('tight')
    ax1.set_title('%s' % var, fontsize=15)
    ax1.yaxis.set_major_locator(MaxNLocator(prune='lower'))
    #ax1.grid(True)
    ax1.xaxis.set_visible(False)
    #ax1.set_yticks([0,0.5,1])
    ax3.yaxis.label.set_color('r')
    ax1.yaxis.label.set_color('b')
    ax1.tick_params(axis='y', colors='b')
    ax3.tick_params(axis='y', colors='r')

    ax1.set_xlim(0,1000)
    
    phs_dt = round(len(time) / 12)
    tidx = np.arange(np.max(np.floor(phs_dt / 2)), len(time), phs_dt)
    tidx = [int(i) for i in tidx]
    tidx = np.array(tidx)
    phs_dp = round(len(result['period']) / 35)
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
    if len(extra_sig) == 0:
        p = 1
    else:
        pc3 = ax2.contour(time, np.log2(result['period']), extra_sig,[-99, 1], linewidths=2, colors = 'greenyellow')
    # cone-of-influence , anything "below"is dubious
    ax2.plot(time, np.log2(result['coi']), 'k')
    ax2.fill_between(time, np.log2(result['coi']), int(
        np.log2(result['period'][-1]) + 1), alpha=0.5, hatch='/')
    position = fig.add_axes([0.94, 0.2, 0.005, 0.4])
    # ,norm=normal)#,shrink=0.5,pad=0.08)
    cbar = plt.colorbar(CS, cax = position, orientation='vertical')
    cbar.set_label('Power', fontsize = 8)
    Q = ax4.quiver(X.astype(np.int64), Y, U, V, linewidth=0.1)
    yt = (np.log2(result['period']))
    yt = range(int(np.log2(result['period'][0])), int(
        np.log2(result['period'][-1]) + 1))# create the vector of periods
    Yticks = [float(math.pow(2, p)) for p in yt]# make 2^periods
    Yticks = [int(i) for i in Yticks]
    ax2.set_yticks(yt)
    ax2.set_yticklabels(Yticks)
    ax2.set_ylim(ymin=(np.log2(result['period'][i1])), ymax=(
        np.log2(result['period'][i2])))
    ax2.set_xlabel('Year Number', fontsize=13)
    ax2.set_ylabel('Period (years)', fontsize=13)
    ax2.invert_yaxis()
    plt.tight_layout()
    plt.savefig(figname, dpi= 300)
    plt.show()
    return
