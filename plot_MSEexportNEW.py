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
#domsizes = [768]
#colors = ['k', 'r', 'g']
colors = ['k', 'r', 'g', 'm']
#colors = ['k', 'r']
#colors=['m']

#domsizes = [768, 1536, 3072]
#colors = ['r', 'g']

N = 5*4

divh_cs = []
tdays = []
t3Ds = []
A_cs = []

t3Ds_plot = np.arange(0, 350, 0.25)
tdays_plot = np.arange(0, 351.0, 1.0)

for j, domsize in enumerate(domsizes):
    
    

    
    dhhatdt_plot = np.ones(t3Ds_plot.shape)*np.nan
    qnet_plot = np.ones(t3Ds_plot.shape)*np.nan
    LHF_plot = np.ones(t3Ds_plot.shape)*np.nan
    SHF_plot = np.ones(t3Ds_plot.shape)*np.nan
    divh_plot = np.ones(t3Ds_plot.shape)*np.nan
    THF_plot = np.ones(t3Ds_plot.shape)*np.nan
    hprimesq_plot = np.ones(t3Ds_plot.shape)*np.nan
    
    divhs_norm_plot = np.ones(t3Ds_plot.shape)*np.nan
    THFs_norm_plot = np.ones(t3Ds_plot.shape)*np.nan
    qnets_norm_plot = np.ones(t3Ds_plot.shape)*np.nan
    divhs_norm_plot_smooth = np.ones(t3Ds_plot.shape)*np.nan
    divhs_plot_smooth = np.ones(t3Ds_plot.shape)*np.nan
    #P_c = 0.0

    
    fnames = glob.glob(fdata + '{:d}km*hprimedivh_*.npz'.format(domsize))
    #fnamesd = glob.glob(fdata + '{:d}km*div_h_d_*_{:1.1f}mm.npz'.format(domsize, P_c))
    divhs = []
    t3Ds = []
    #A_cs = []
    qnets = []
    LHFs = []
    SHFs = []
    dhhatdts = []
    hprimesqs = []

    
    print fnames
    
    end = len(fnames) - 1
    
    for k, fname in enumerate(fnames):

        
        divh_f = np.load(fname)
        #divh_df = np.load(fnamesd[k])
        
        dhhatdt = divh_f['hdhhatdt'][:]
        #hhat = divh_cf['hhatbar'][:]
        qnet = divh_f['hqnet'][:] 
        LHF = divh_f['hLHF'][:]
        SHF = divh_f['hSHF'][:]
        divh = divh_f['hdivh'][:]
        hprimesq = divh_f['hprimesq'][:]
        t3D = divh_f['time'][:]
        t3Ds = np.hstack([t3Ds, t3D])
        
     
        qnets = np.hstack([qnets, qnet])
        LHFs = np.hstack([LHFs, LHF])
        SHFs = np.hstack([SHFs, SHF])
        divhs = np.hstack([divhs, divh])
        dhhatdts = np.hstack([dhhatdts, dhhatdt])
        hprimesqs = np.hstack([hprimesqs, hprimesq])
        
        t3D = divh_f['time'][:]
        tday = np.unique(np.round(t3D))
        tdays = np.hstack([tdays, tday])

        
        fnamelarge = glob.glob(fdata + '*hprimedivh_day340to350*.npz')
        if fnamelarge:
            if (fname == fnamelarge[0]) and (k == end):
                clip=20
                A_cs = A_cs[:-clip]
                t3Ds = t3Ds[:-clip]
                qnets = qnets[:-clip]
                divhs = divhs[:-clip]
                LHFs = LHFs[:-clip]
                SHFs = SHFs[:-clip]
                #dhhatdt = dhhatdts[:-clip]
                hprimesqs = hprimesqs[:-clip]

        
        
    #A_ds = (domsize*1e3)**2 - A_cs    
    #smooth_divh = running_mean(divh_c, N)
    
    toplot = np.in1d(t3Ds_plot, t3Ds)
    toplot_days = np.in1d(tdays_plot, tdays)
    
    #qnet_cs = qnet_TOA_cs - qnet_s_cs
    
    qnet_plot[toplot] = qnets
    divh_plot[toplot] = divhs
    #dhhatdt_plot[toplot] = dhhatdts
    hprimesq_plot[toplot] = hprimesqs
    LHF_plot[toplot] = LHFs
    SHF_plot[toplot] = SHFs
    THF_plot[toplot] = LHFs + SHFs
    

    
    

    # plt.figure(1)
    # plt.plot(t3Ds_plot, divh_cs_plot, color=colors[j], linewidth=2)
    # plt.axhline(0, color='b')
    # plt.xlabel('time (days)')
    # plt.ylabel('MSE export from convective region (J/s)')
    # plt.savefig(fout + 'MSEexport_conv.pdf')
    
    # plt.figure(2)
    # plt.plot(t3Ds_plot, qnet_s_cs_plot, color=colors[j], linewidth=2)
    # plt.axhline(0, color='b')
    # plt.xlabel('time (days)')
    # plt.ylabel(r'$Q_{net,surf}$ (W/m$^{2}$)')
    # plt.title(r'$Q_{net,surf}$ ')
    # plt.savefig(fout + 'Qnets_conv.pdf')
    # 
    # plt.figure(3)
    # plt.plot(t3Ds_plot, qnet_TOA_cs_plot, color=colors[j], linewidth=2)
    # plt.axhline(0, color='b')
    # plt.xlabel('time (days)')
    # plt.ylabel(r'$Q_{net,TOA}$ (W/m$^{2}$)')
    # plt.title(r'$Q_{net,TOA}$')
    # plt.savefig(fout + 'QnetTOA_conv.pdf')
    
    # plt.figure(4)
    # plt.plot(t3Ds_plot, dhhatdt_cs_plot, color=colors[j], linewidth=2)
    # plt.axhline(0, color='b')
    # plt.xlabel('time (days)')
    # plt.ylabel(r'<$\partial h/ \partial t$> (W/m$^{2}$)')
    # plt.title(r'<$\partial h/ \partial t$>')
    # plt.savefig(fout + 'dhhatdt_conv.pdf')
    
    # plt.figure(5)
    # plt.plot(t3Ds_plot, dhhatdt_ds_plot, color=colors[j], linewidth=2)
    # plt.axhline(0, color='b')
    # plt.xlabel('time (days)')
    # plt.ylabel(r'<$\partial h/ \partial t$> (W/m$^{2}$)')
    # plt.title(r'<$\partial h/ \partial t$>')
    # plt.savefig(fout + 'dhhatdt_dry.pdf')
    
    # plt.figure(6)
    # plt.plot(t3Ds_plot, qnet_TOA_cs_plot - qnet_s_cs_plot, color=colors[j], linewidth=2)
    # plt.axhline(0, color='b')
    # plt.xlabel('time (days)')
    # plt.ylabel(r'$Q_{net}$ (W/m$^{2}$)')
    # plt.title(r'$Q_{net}$')
    # plt.savefig(fout + 'Qnet_conv.pdf')

    
    # plt.figure(7)
    # plt.plot(t3Ds_plot, dhhatdt_cs_plot*A_cs_plot + dhhatdt_ds_plot*A_ds_plot, color=colors[j], linewidth=2)
    # plt.axhline(0, color='b')
    # plt.xlabel('time (days)')
    # plt.ylabel(r'<$\partial h/ \partial t$> (J/s)')
    # plt.title(r'<$\partial h/ \partial t$>')
    # plt.savefig(fout + 'dhhatdt_sum.pdf')
    
    # plt.figure(8)
    # plt.plot(t3Ds_plot, divh_cs_plot + divh_ds_plot, color=colors[j], linewidth=2)
    # plt.axhline(0, color='b')
    # plt.xlabel('time (days)')
    # plt.ylabel('total MSE export (J/s)')
    # plt.savefig(fout + 'MSEexport_sum.pdf')
    
    
    # plt.figure(9)
    # plt.plot(t3Ds_plot, hhat_plot, color=colors[j], linewidth=2)
    # #plt.axhline(0, color='b')
    # plt.xlabel('time (days)')
    # plt.ylabel('<$h$> (J/m$^{2}$)')
    # plt.savefig(fout + 'hhat.pdf')
    # 
    # plt.figure(10)
    # plt.plot(t3Ds_plot, Enet_plot, color=colors[j], linewidth=2)
    # plt.axhline(0, color='b')
    # plt.xlabel('time (days)')
    # plt.ylabel('$Q_{net}$ + LHF + SHF (W/m$^{2}$)')
    # plt.savefig(fout + 'Enet.pdf')
    # 
    # plt.figure(11)
    # plt.plot(t3Ds_plot, dhhatdt_plot, color=colors[j], linewidth=2)
    # plt.axhline(0, color='b')
    # plt.xlabel('time (days)')
    # plt.ylabel(r'<$\partial h/ \partial t$> (W/m$^{2}$)')
    # plt.savefig(fout + 'dhhatdt.pdf')


    divhs_norm = divhs/hprimesqs
    qnets_norm = qnets/hprimesqs
    LHFs_norm = LHFs/hprimesqs
    SHFs_norm = SHFs/hprimesqs
    #dhhatdt_norm = dhhatdts/hprimesqs
    
   # A_ds = (domsize*1e3)**2 - A_cs 
    
    
    #divh_sums_norm = (divh_cs + divh_ds)/((domsize*1e3)**2)
    

    
    divh_norm_padded = np.pad(divhs_norm, (N//2, N-1-N//2), mode='edge')
    divh_padded = np.pad(divhs, (N//2, N-1-N//2), mode='edge')
    
    THFs_norm = LHFs_norm + SHFs_norm
    #THFs_norm_padded = np.pad(THFs_norm, (N//2, N-1-N//2), mode='edge')
    #qnet_sums_norm_padded = np.pad(divh_sums_norm, (N//2, N-1-N//2), mode='edge')
    #divh_sums_padded = np.pad(divh_cs + divh_ds, (N//2, N-1-N//2), mode='edge')
    
    #Enet_cs_padded = np.pad(Enet_cs, (N//2, N-1-N//2), mode='edge')
    #qnet_cs_padded = np.pad(qnet_cs, (N//2, N-1-N//2), mode='edge')
    
    #dhhatdt_cs_padded = np.pad(dhhatdt_cs, (N//2, N-1-N//2), mode='edge')
    #dhhatdt_sums_padded = np.pad(dhhatdt_cs*A_cs + dhhatdt_ds*A_ds, (N//2, N-1-N//2), mode='edge')
    
    #smooth_Enet_cs = np.convolve(Enet_cs_padded, np.ones((N,))/N, mode='valid') 
    #smooth_qnet_cs = np.convolve(qnet_cs_padded, np.ones((N,))/N, mode='valid') 
    #smooth_dhhatdt_cs = np.convolve(dhhatdt_cs_padded, np.ones((N,))/N, mode='valid') 
    #smooth_dhhatdt_sums = np.convolve(dhhatdt_sums_padded, np.ones((N,))/N, mode='valid') 
    #smooth_divh_cs_norm = np.convolve(divh_cs_norm_padded, np.ones((N,))/N, mode='valid') 
    #smooth_divh_sums_norm = np.convolve(divh_sums_norm_padded, np.ones((N,))/N, mode='valid') 
    smooth_divh_norms = np.convolve(divh_norm_padded, np.ones((N,))/N, mode='valid') 
    smooth_divhs = np.convolve(divh_padded, np.ones((N,))/N, mode='valid') 
    
    
    
    if k == end:
        print 'domsize', domsize
        print '{3.0f} km average normalized $M$ (20-day ave):', np.nanmean(-3600*24*smooth_divh_norms[-20*4:]) 
        print '{3.0f} km average normalized $Q_{net}$ (20-day ave):', np.nanmean(-3600*24*qnets_norm[-20*4:]) 
        print '{3.0f} km average normalized $LHF$ (20-day ave):', np.nanmean(-3600*24*LHFs_norm[-20*4:]) 
    
    
    nt = len(divhs_norm)
    # ntrunc = nt % 4
    # divh_cs_norm_temp = divh_cs_norm[ntrunc:]
    # 
    # divh_cs_norm_temp = divh_cs_norm_temp.reshape(nt/4, 4)
    # divh_cs_norm_daily = np.mean(divh_cs_norm_temp, axis=1)
    
    #tdailyplot = np.arange(0, 250, len(divh_cs_norm_daily))

    


    #smooth_divh_cs_norm = running_mean(divh_cs_norm, N)
    
    #divh_cs_norm_daily_plot[toplot_days[:-4]] = divh_cs_norm_daily
    

    divhs_norm_plot[toplot] = divhs_norm
    divhs_norm_plot_smooth[toplot] = smooth_divh_norms
    divhs_plot_smooth[toplot] = smooth_divhs
    
    THFs_norm_plot[toplot] = THFs_norm
    qnets_norm_plot[toplot] = qnets_norm
    
    # 
    # divh_sums_plot_smooth[toplot] =  smooth_divh_sums
    # 
    # divh_cs_norm_plot_smooth[toplot] = smooth_divh_cs_norm
    # 
    # divh_sums_norm_plot_smooth[toplot] = smooth_divh_sums_norm
    # 
    # dhhatdt_cs_plot_smooth[toplot] = smooth_dhhatdt_cs
    # 
    # Enet_cs_plot_smooth[toplot] = smooth_Enet_cs
    # 
    # qnet_cs_plot_smooth[toplot] = smooth_qnet_cs
    # 
    # dhhatdt_sums_plot_smooth[toplot] = smooth_dhhatdt_sums
    
    #divh_cs_norm_plot_smooth[toplot_smooth] = smooth_divh_cs_norm

    plt.figure(1)
    plt.plot(t3Ds_plot, -divh_plot, color=colors[j], linewidth=1, alpha=0.5)
    plt.plot(t3Ds_plot,  -divhs_plot_smooth, color=colors[j],  linewidth=3,  label='{:d} km'.format(domsize))
    plt.axhline(0, color='black', linewidth=1)
    plt.xlabel('time (days)')
    plt.ylabel(r'$M$ (J$^2$ m$^{-4}$ s$^{-1}$)')
    plt.title(r'$M$')
    #plt.ylim(-150, 150)
    plt.legend()
    plt.savefig(fout + 'hprimedivh_smooth.pdf')
    
    plt.figure(2)
    plt.plot(t3Ds_plot, -(3600*24)*divhs_norm_plot, color=colors[j], linewidth=1, alpha=0.5)
    plt.plot(t3Ds_plot,  -(3600*24)*divhs_norm_plot_smooth, color=colors[j],  linewidth=3, label='{:d} km'.format(domsize))
    plt.axhline(0, color='black', linewidth=1)
    plt.xlabel('time (days)')
    plt.ylabel(r'$M$ (day$^{-1}$)')
    plt.title(r'$M$')
    plt.ylim(-0.4, 0.4)
    # plt.xlim(90,220)
    # plt.axvline(100, color='black', linewidth=1, alpha=0.5)
    # plt.axvline(120, color='black', linewidth=1, alpha=0.5)
    # plt.axvline(140, color='black', linewidth=1, alpha=0.5)
    # plt.axvline(160, color='black', linewidth=1, alpha=0.5)
    # plt.axvline(180, color='black', linewidth=1, alpha=0.5)
    # plt.axvline(200, color='black', linewidth=1, alpha=0.5)
    plt.legend()
    plt.savefig(fout + 'hprimedivh_norm_smooth.pdf')
    
    plt.figure(3)
    plt.plot(t3Ds_plot, hprimesq_plot, color=colors[j], linewidth=2, label='{:d} km'.format(domsize))
    #plt.plot(t3Ds_plot,  divhs_norm_plot_smooth, color=colors[j], label='{:d} km'.format(domsize))
    #plt.axhline(0, color='black', linewidth=1)
    plt.xlabel('time (days)')
    plt.ylabel(r'$\langle h^{\prime} \rangle^2$ (J$^2$ m$^{-4}$)')
    plt.title(r'$\langle h^{\prime} \rangle^2$')
    #plt.ylim(-1, 1)
    plt.legend()
    plt.savefig(fout + 'hprimesq.pdf')
    
    plt.figure(4)
    #plt.plot(t3Ds_plot, -(3600*24)*divhs_norm_plot, color=colors[j], linewidth=1, alpha=0.5)
    plt.plot(t3Ds_plot, (3600*24)*qnets_norm_plot, color=colors[j], linewidth=2, label='{:d} km'.format(domsize))
    #plt.plot(t3Ds_plot, (3600*24)*THFs_norm_plot, label = r'$\langle h^{\prime} \rangle THF^{\prime}$')               
    #plt.plot(t3Ds_plot,  -(3600*24)*divhs_norm_plot_smooth, color=colors[j], label= r'$M$')
    plt.axhline(0, color='black', linewidth=1)
    plt.xlabel('time (days)')
    plt.ylabel(r'$\langle h^{\prime} \rangle Q^{\prime}_{net}$ (day$^{-1}$)')
    plt.title(r'$\langle h^{\prime} \rangle Q^{\prime}_{net}$')
    plt.ylim(-0.4, 0.4)
    plt.legend()
    plt.savefig(fout + 'hprimeQnet_norm.pdf')
    
    plt.figure(5)
    #plt.plot(t3Ds_plot, -(3600*24)*divhs_norm_plot, color=colors[j], linewidth=1, alpha=0.5)
    plt.plot(t3Ds_plot, (3600*24)*THFs_norm_plot, color=colors[j], linewidth=2, label='{:d} km'.format(domsize))
    #plt.plot(t3Ds_plot, (3600*24)*THFs_norm_plot, label = r'$\langle h^{\prime} \rangle THF^{\prime}$')               
    #plt.plot(t3Ds_plot,  -(3600*24)*divhs_norm_plot_smooth, color=colors[j], label= r'$M$')
    plt.axhline(0, color='black', linewidth=1)
    plt.xlabel('time (days)')
    plt.ylabel(r'$\langle h^{\prime} \rangle THF^{\prime}$ (day$^{-1}$)')
    plt.title(r'$\langle h^{\prime} \rangle THF^{\prime}$')
    plt.ylim(-0.4, 0.4)
    plt.legend()
    plt.savefig(fout + 'hprimeTHF_norm.pdf')
    
    
plt.close('all')


