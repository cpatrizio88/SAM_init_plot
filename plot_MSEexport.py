import numpy as np
import matplotlib.pyplot as plt
import glob
import matplotlib

matplotlib.rcParams.update({'font.size': 28})
matplotlib.rcParams.update({'figure.figsize': (18, 10)})
matplotlib.rcParams.update({'lines.linewidth': 3})
matplotlib.rcParams.update({'legend.fontsize': 22})
matplotlib.rcParams.update({'mathtext.fontset': 'cm'})

def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / float(N)
    
fdata='/Users/cpatrizio/data/SAM/'
fout='/Volumes/GoogleDrive/My Drive/MS/MS figures/SST302/varydomsize_SAM_aggr_64vert_ubarzero_STAT/'

domsizes = [768, 1536, 3072, 6144]
domsizes = [768]
colors = ['k', 'r', 'g', 'm']
colors = ['k']

N = 10*4

divh_cs = []
tdays = []
t3Ds = []
A_cs = []

t3Ds_plot = np.arange(0, 350, 0.25)
tdays_plot = np.arange(0, 351.0, 1.0)

for j, domsize in enumerate(domsizes):
    
    
    divh_cs_plot = np.ones(t3Ds_plot.shape)*np.nan
    divh_ds_plot = np.ones(t3Ds_plot.shape)*np.nan
    hhat_plot = np.ones(t3Ds_plot.shape)*np.nan
    dhhatdt_plot = np.ones(t3Ds_plot.shape)*np.nan
    Enet_plot = np.ones(t3Ds_plot.shape)*np.nan
    qnet_s_cs_plot = np.ones(t3Ds_plot.shape)*np.nan
    qnet_TOA_cs_plot = np.ones(t3Ds_plot.shape)*np.nan
    dhhatdt_cs_plot = np.ones(t3Ds_plot.shape)*np.nan
    dhhatdt_ds_plot = np.ones(t3Ds_plot.shape)*np.nan
    A_cs_plot = np.ones(t3Ds_plot.shape)*np.nan
    A_ds_plot = np.ones(t3Ds_plot.shape)*np.nan
    divh_cs_norm_plot = np.ones(t3Ds_plot.shape)*np.nan
    divh_sums_norm_plot = np.ones(t3Ds_plot.shape)*np.nan
    divh_cs_norm_plot_smooth = np.ones(t3Ds_plot.shape)*np.nan
    divh_sums_norm_plot_smooth = np.ones(t3Ds_plot.shape)*np.nan
    dhhatdt_cs_plot_smooth =  np.ones(t3Ds_plot.shape)*np.nan
    dhhatdt_sums_plot_smooth =  np.ones(t3Ds_plot.shape)*np.nan
    divh_sums_plot_smooth = np.ones(t3Ds_plot.shape)*np.nan
    divh_cs_norm_daily_plot = np.ones(tdays_plot.shape)*np.nan
    
    fnames = glob.glob(fdata + '{:d}km*div_h_c_*0.0mm.npz'.format(domsize))
    fnamesd = glob.glob(fdata + '{:d}km*div_h_d_*0.0mm.npz'.format(domsize))
    divh_cs = []
    divh_ds = []
    t3Ds = []
    A_cs = []
    qnet_s_cs = []
    qnet_TOA_cs = []
    dhhatdt_cs = []
    dhhatdt_ds = []
    Enets = []
    hhats = []
    dhhatdts = []
    
    print fnames
    
    end = len(fnames) - 1
    
    for k, fname in enumerate(fnames):

        
        divh_cf = np.load(fname)
        divh_df = np.load(fnamesd[k])
        
        dhhatdt = divh_cf['dhhatdtbar'][:]
        hhat = divh_cf['hhatbar'][:]
        Enet = divh_cf['Enetbar'][:]
        
        divh_c = divh_cf['divh'][:]
        divh_d = divh_df['divh'][:]
        qnet_s_c = divh_cf['qnets'][:]
        qnet_TOA_c = divh_cf['qnetTOA'][:]
        dhhatdt_c = divh_cf['dhhatdt'][:]
        dhhatdt_d = divh_df['dhhatdt'][:]
        
        dhhatdts = np.hstack([dhhatdts, dhhatdt])
        hhats = np.hstack([hhats, hhat])
        Enets = np.hstack([Enets, Enet])
        divh_cs = np.hstack([divh_cs, divh_c])
        divh_ds = np.hstack([divh_ds, divh_d])
        qnet_s_cs = np.hstack([qnet_s_cs, qnet_s_c])
        qnet_TOA_cs = np.hstack([qnet_TOA_cs, qnet_TOA_c])
        dhhatdt_cs = np.hstack([dhhatdt_cs, dhhatdt_c])
        dhhatdt_ds = np.hstack([dhhatdt_ds, dhhatdt_d])
        t3D = divh_cf['time'][:]
        t3Ds = np.hstack([t3Ds, t3D])
        tday = np.unique(np.round(t3D))
        tdays = np.hstack([tdays, tday])
        A_c = divh_cf['A'][:]
        A_cs = np.hstack([A_cs, A_c])
        
        fnamelarge = glob.glob(fdata + '*day340to350*.npz')
        if fnamelarge:
            if (fname == fnamelarge) and (k == end):
                A_cs = A_cs[:-20]
                t3Ds = t3Ds[:-20]
                divh_cs = divh_cs[:-20]
                divh_ds = divh_ds[:-20]
                qnet_s_cs = qnet_s_cs[:-20]
                qnet_TOA_cs = qnet_TOA_cs[:-20]
                dhhatdt_cs = dhhatdt_cs[:-20]
                dhhatdt_ds = dhhatdt_ds[:-20]
        
        
    A_ds = (domsize*1e3)**2 - A_cs    
    smooth_divh_c = running_mean(divh_c, N)
    
    toplot = np.in1d(t3Ds_plot, t3Ds)
    toplot_days = np.in1d(tdays_plot, tdays)
    
    Enet_plot[toplot] = Enets
    dhhatdt_plot[toplot] = dhhatdts
    hhat_plot[toplot] = hhats
    
    divh_cs_plot[toplot] = divh_cs
    divh_ds_plot[toplot] = divh_ds
    qnet_s_cs_plot[toplot] = qnet_s_cs
    qnet_TOA_cs_plot[toplot] = qnet_TOA_cs
    dhhatdt_cs_plot[toplot] = dhhatdt_cs
    dhhatdt_ds_plot[toplot] = dhhatdt_ds
    A_cs_plot[toplot] = A_cs
    A_ds_plot[toplot] = A_ds 
    
    

    plt.figure(1)
    plt.plot(t3Ds_plot, divh_cs_plot, color=colors[j], linewidth=2)
    plt.axhline(0, color='b')
    plt.xlabel('time (days)')
    plt.ylabel('MSE export from convective region (J/s)')
    plt.savefig(fout + 'MSEexport_conv.pdf')
    
    plt.figure(2)
    plt.plot(t3Ds_plot, qnet_s_cs_plot, color=colors[j], linewidth=2)
    plt.axhline(0, color='b')
    plt.xlabel('time (days)')
    plt.ylabel(r'$Q_{net,surf}$ (W/m$^{2}$)')
    plt.title(r'$Q_{net,surf}$ ')
    plt.savefig(fout + 'Qnets_conv.pdf')
    
    plt.figure(3)
    plt.plot(t3Ds_plot, qnet_TOA_cs_plot, color=colors[j], linewidth=2)
    plt.axhline(0, color='b')
    plt.xlabel('time (days)')
    plt.ylabel(r'$Q_{net,TOA}$ (W/m$^{2}$)')
    plt.title(r'$Q_{net,TOA}$')
    plt.savefig(fout + 'QnetTOA_conv.pdf')
    
    plt.figure(4)
    plt.plot(t3Ds_plot, dhhatdt_cs_plot, color=colors[j], linewidth=2)
    plt.axhline(0, color='b')
    plt.xlabel('time (days)')
    plt.ylabel(r'<$\partial h/ \partial t$> (W/m$^{2}$)')
    plt.title(r'<$\partial h/ \partial t$>')
    plt.savefig(fout + 'dhhatdt_conv.pdf')
    
    plt.figure(5)
    plt.plot(t3Ds_plot, dhhatdt_ds_plot, color=colors[j], linewidth=2)
    plt.axhline(0, color='b')
    plt.xlabel('time (days)')
    plt.ylabel(r'<$\partial h/ \partial t$> (W/m$^{2}$)')
    plt.title(r'<$\partial h/ \partial t$>')
    plt.savefig(fout + 'dhhatdt_dry.pdf')
    
    plt.figure(6)
    plt.plot(t3Ds_plot, qnet_TOA_cs_plot - qnet_s_cs_plot, color=colors[j], linewidth=2)
    plt.axhline(0, color='b')
    plt.xlabel('time (days)')
    plt.ylabel(r'$Q_{net}$ (W/m$^{2}$)')
    plt.title(r'$Q_{net}$')
    plt.savefig(fout + 'Qnet_conv.pdf')

    
    plt.figure(7)
    plt.plot(t3Ds_plot, dhhatdt_cs_plot*A_cs_plot + dhhatdt_ds_plot*A_ds_plot, color=colors[j], linewidth=2)
    plt.axhline(0, color='b')
    plt.xlabel('time (days)')
    plt.ylabel(r'<$\partial h/ \partial t$> (J/s)')
    plt.title(r'<$\partial h/ \partial t$>')
    plt.savefig(fout + 'dhhatdt_sum.pdf')
    
    plt.figure(8)
    plt.plot(t3Ds_plot, divh_cs_plot + divh_ds_plot, color=colors[j], linewidth=2)
    plt.axhline(0, color='b')
    plt.xlabel('time (days)')
    plt.ylabel('total MSE export (J/s)')
    plt.savefig(fout + 'MSEexport_sum.pdf')
    
    
    plt.figure(9)
    plt.plot(t3Ds_plot, hhat_plot, color=colors[j], linewidth=2)
    #plt.axhline(0, color='b')
    plt.xlabel('time (days)')
    plt.ylabel('<$h$> (J/m$^{2}$)')
    plt.savefig(fout + 'hhat.pdf')
    
    plt.figure(10)
    plt.plot(t3Ds_plot, Enet_plot, color=colors[j], linewidth=2)
    plt.axhline(0, color='b')
    plt.xlabel('time (days)')
    plt.ylabel('$Q_{net}$ + LHF + SHF (W/m$^{2}$)')
    plt.savefig(fout + 'Enet.pdf')
    
    plt.figure(11)
    plt.plot(t3Ds_plot, dhhatdt_plot, color=colors[j], linewidth=2)
    plt.axhline(0, color='b')
    plt.xlabel('time (days)')
    plt.ylabel(r'<$\partial h/ \partial t$> (W/m$^{2}$)')
    plt.savefig(fout + 'dhhatdt.pdf')


    divh_cs_norm = divh_cs/A_cs
    
    A_ds = (domsize*1e3)**2 - A_cs 
    
    divh_sums_norm = (divh_cs + divh_ds)/((domsize*1e3)**2)
    
    divh_cs_norm_padded = np.pad(divh_cs_norm, (N//2, N-1-N//2), mode='edge')
    divh_sums_norm_padded = np.pad(divh_sums_norm, (N//2, N-1-N//2), mode='edge')
    divh_sums_padded = np.pad(divh_cs + divh_ds, (N//2, N-1-N//2), mode='edge')
    
    
    dhhatdt_cs_padded = np.pad(dhhatdt_cs, (N//2, N-1-N//2), mode='edge')
    dhhatdt_sums_padded = np.pad(dhhatdt_cs*A_cs + dhhatdt_ds*A_ds, (N//2, N-1-N//2), mode='edge')
    
    smooth_dhhatdt_cs = np.convolve(dhhatdt_cs_padded, np.ones((N,))/N, mode='valid') 
    smooth_dhhatdt_sums = np.convolve(dhhatdt_sums_padded, np.ones((N,))/N, mode='valid') 
    smooth_divh_cs_norm = np.convolve(divh_cs_norm_padded, np.ones((N,))/N, mode='valid') 
    smooth_divh_sums_norm = np.convolve(divh_sums_norm_padded, np.ones((N,))/N, mode='valid') 
    smooth_divh_sums = np.convolve(divh_sums_padded, np.ones((N,))/N, mode='valid') 
    
    
    nt = len(divh_cs_norm)
    ntrunc = nt % 4
    divh_cs_norm_temp = divh_cs_norm[ntrunc:]
    
    divh_cs_norm_temp = divh_cs_norm_temp.reshape(nt/4, 4)
    divh_cs_norm_daily = np.mean(divh_cs_norm_temp, axis=1)
    
    #tdailyplot = np.arange(0, 250, len(divh_cs_norm_daily))

    


    #smooth_divh_cs_norm = running_mean(divh_cs_norm, N)
    
    #divh_cs_norm_daily_plot[toplot_days[:-4]] = divh_cs_norm_daily
    

    divh_cs_norm_plot[toplot] = divh_cs_norm
    
    divh_sums_norm_plot[toplot] = divh_sums_norm
    
    
    
    divh_sums_plot_smooth[toplot] =  smooth_divh_sums
    
    divh_cs_norm_plot_smooth[toplot] = smooth_divh_cs_norm
    
    divh_sums_norm_plot_smooth[toplot] = smooth_divh_sums_norm
    
    dhhatdt_cs_plot_smooth[toplot] = smooth_dhhatdt_cs
    
    dhhatdt_sums_plot_smooth[toplot] = smooth_dhhatdt_sums
    
    #divh_cs_norm_plot_smooth[toplot_smooth] = smooth_divh_cs_norm

    plt.figure(12)
    plt.plot(t3Ds_plot, divh_cs_norm_plot-5, color=colors[j], linewidth=1, alpha=0.5)
    plt.plot(t3Ds_plot,  divh_cs_norm_plot_smooth-5, color=colors[j], label='{:d} km'.format(domsize))
    #plt.plot(tdays_plot, divh_cs_norm_daily_plot, color=colors[j], label='{:d} km'.format(domsize))
    #plt.plot(t3Ds_plot,  divh_cs_norm_plot_smooth, color=colors[j])
    plt.axhline(0, color='b')
    plt.xlabel('time (days)')
    plt.ylabel(r'$M$ (W m$^{-2}$)')
    plt.title(r'$M$')
    plt.ylim(-100, 100)
    plt.legend()
    plt.savefig(fout + 'MSEexportnorm_conv_smooth.pdf')
    
    
    
    plt.figure(13)
    plt.plot(t3Ds_plot, dhhatdt_cs_plot, color=colors[j], linewidth=1, alpha=0.5)
    plt.plot(t3Ds_plot,  dhhatdt_cs_plot_smooth, color=colors[j], label='{:d} km'.format(domsize))
    plt.axhline(0, color='b')
    plt.xlabel('time (days)')
    plt.ylabel(r'<$\partial h/ \partial t$> (W/m$^{2}$)')
    plt.title(r'<$\partial h/ \partial t$>')
    plt.ylim(-100, 100)
    plt.legend()
    plt.savefig(fout + 'dhhatdt_conv_smooth.pdf')
    
    
    plt.figure(14)
    plt.plot(t3Ds_plot, divh_sums_norm_plot, color=colors[j], linewidth=1, alpha=0.5)
    plt.plot(t3Ds_plot,  divh_sums_norm_plot_smooth, color=colors[j], label='{:d} km'.format(domsize))
    plt.axhline(0, color='b')
    plt.xlabel('time (days)')
    plt.ylabel(r'$M$ (W m$^{-2}$)')
    plt.title(r'$M$')
    #plt.ylim(-100, 180)
    plt.legend()
    plt.savefig(fout + 'MSEexportnorm_smooth.pdf')
    
    # plt.figure(15)
    # #plt.plot(t3Ds_plot, divh_cs_norm_plot, color=colors[j], linewidth=2, alpha=0.5)
    # #plt.plot(t3Ds_plot,  divh_cs_norm_plot_smooth, color=colors[j], label='{:d} km'.format(domsize))
    # plt.plot(tdays_plot, divh_cs_norm_daily_plot, color=colors[j], label='{:d} km'.format(domsize))
    # #plt.plot(t3Ds_plot,  divh_cs_norm_plot_smooth, color=colors[j])
    # plt.axhline(0, color='b')
    # plt.xlabel('time (days)')
    # plt.ylabel(r'$M$ (W m$^{-2}$)')
    # plt.title(r'$M$')
    # #plt.ylim(-100, 180)
    # plt.legend()
    # plt.savefig(fout + 'MSEexportnorm_conv_daily.pdf')
    
    # plt.figure(16)
    # plt.plot(t3Ds_plot, divh_cs_norm_plot, color=colors[j], linewidth=2, alpha=0.5)
    # plt.plot(t3Ds_plot,  divh_sums_plot_smooth, color=colors[j], label='{:d} km'.format(domsize))
    # #plt.plot(tdays_plot, divh_cs_norm_daily_plot, color=colors[j], label='{:d} km'.format(domsize))
    # #plt.plot(t3Ds_plot,  divh_cs_norm_plot_smooth, color=colors[j])
    # plt.axhline(0, color='b')
    # plt.xlabel('time (days)')
    # plt.ylabel(r'total MSE export (J/s)')
    # plt.title(r'total MSE export')
    # plt.ylim(-100, 100)
    # plt.legend()
    # plt.savefig(fout + 'MSEexport_sum_smooth.pdf')
    # 
    
    # plt.figure(17)
    # plt.plot(t3Ds_plot, dhhatdt_plot, color=colors[j], linewidth=2, alpha=0.5)
    # plt.plot(t3Ds_plot,  dhhatdt_sums_plot_smooth, color=colors[j], label='{:d} km'.format(domsize))
    # #plt.plot(tdays_plot, divh_cs_norm_daily_plot, color=colors[j], label='{:d} km'.format(domsize))
    # #plt.plot(t3Ds_plot,  divh_cs_norm_plot_smooth, color=colors[j])
    # plt.axhline(0, color='b')
    # plt.xlabel('time (days)')
    # plt.ylabel(r'<$\partial h/ \partial t$> (J/s)')
    # plt.title(r'<$\partial h/ \partial t$>')
    # plt.savefig(fout + 'dhhatdt_sum_smooth.pdf')
    
    # plt.figure(18)
    # #plt.plot(t3Ds_plot, divh_cs_norm_plot, color=colors[j], linewidth=2, alpha=0.5)
    # plt.plot(t3Ds_plot,  (divh_sums_plot_smooth - dhhatdt_sums_plot_smooth)/(1e3*domsize)**2, color=colors[j], label='{:d} km'.format(domsize))
    # #plt.plot(tdays_plot, divh_cs_norm_daily_plot, color=colors[j], label='{:d} km'.format(domsize))
    # #plt.plot(t3Ds_plot,  divh_cs_norm_plot_smooth, color=colors[j])
    # plt.axhline(0, color='b')
    # plt.xlabel('time (days)')
    # plt.ylabel(r'total MSE export (W/m$^2$)')
    # plt.title(r'total MSE export')
    # #plt.ylim(-200, 200)
    # plt.legend()
    # plt.savefig(fout + 'MSEexportnorm_sum_smooth.pdf')
    




plt.close('all')


