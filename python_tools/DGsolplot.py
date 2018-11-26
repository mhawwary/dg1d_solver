from matplotlib import pyplot ,ticker   #and the useful plotting library
from numpy import sin,cos,pi,linspace,ones,zeros,abs,min,max,exp, shape, empty_like , size, loadtxt, arange , log
from decimal import Decimal
import csv
from matplotlib import pyplot as plt

import numpy as np
from scipy.optimize import leastsq, least_squares
import pylab as plt

from fft_toolbox_python import load_data, compute_fft, compute_Etotal\
                , savez_spectrum_linearWave, loadz_spectrum_linearWave, plot_KEnerg\
                , compute_Dt_dissipation_rate, compute_dEkdt, plot_E_vs_time\
                , plot_dissip_rates

pyplot.rc('legend',**{'loc':'best'});
pyplot.rcParams[u'legend.fontsize'] = 18
pyplot.rcParams[u'legend.edgecolor']='black'
pyplot.rcParams[u'font.weight']='normal'
#pyplot.rcParams['font.serif']='false'
pyplot.rcParams[u'xtick.labelsize']=14
pyplot.rcParams[u'ytick.labelsize']=14
pyplot.rcParams[u'axes.titlesize']=18
pyplot.rcParams[u'axes.labelsize']=15
pyplot.rcParams[u'axes.spines.right']='true';
pyplot.rcParams[u'axes.spines.top']='true';
pyplot.rcParams[u'lines.linewidth'] = 1.5;
pyplot.rcParams[u'lines.markersize'] = 8;

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
def compute_fitted_wave(Np, Nelem, x_disc, u_disc):
    N_disc = size(x_disc)  
    N_fit = N_disc-2*Nelem 
    x_fit = zeros(N_fit)
    u_fit = zeros(N_fit)

    j=0;
    for i in range(0, N_disc, Np): # plot all the elements curves except the last one to be able to define a label
        xx = x_disc[i:i + Np];
        uu = u_disc[i:i + Np];
        x_fit[j:Np-2+j] = xx[1:Np-1]
        u_fit[j:Np-2+j] = uu[1:Np-1]
        j=j+Np-2;
        
    guess_mean = np.mean(u_fit)
    guess_std = 3*np.std(u_fit)/(2**0.5)/(2**0.5)
    guess_phase = 0
    guess_freq = 6*pi

    guess_amp = 0.13
    
    x_fit_temp =  x_fit #np.linspace(min(x_cont_exact),max(x_cont_exact),200)
    data_first_guess = guess_std*np.sin(guess_freq*x_fit_temp+guess_phase) + guess_mean
    param = np.zeros(4);
    param = [guess_amp, guess_freq, guess_phase, guess_mean];
    
    optimize_func = lambda x: x[0]*np.sin((x[1]*x_fit)+x[2]) + x[3] - u_fit
    est_amp, est_freq, est_phase, est_mean = leastsq(optimize_func, [guess_amp, guess_freq, guess_phase, guess_mean])[0]
    data_fit = est_amp*np.sin(est_freq*x_fit_temp+est_phase) + est_mean
    
    print('Fitting data:\n')
    print('est_amp: ',est_amp,'\test_freq: ',est_freq,'\test_phase: ',est_phase,'\test_mean: ',est_mean)
    
    return x_fit, u_fit, x_fit_temp, data_fit, data_first_guess
    
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#                                  D I F F U S I ON
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
def plot_initial_proj_diffus(dir1, mode, DG, RK, CFL, Nelem, N_disc_ppt, tt_, dt_ \
                , Epsilon, gamma_, diffus_scheme, cont_num, disc_num, \
                T, K_den,wave_form_):
    nelem = int(Nelem)
    Np1 = int(int(DG)+1)
    dt = float(dt_)

    if  (mode == 'test') | (mode == 'dt_const'):
        mm_name = str('_dt') + dt_
        m_name = str('dt=') + dt_
    elif(mode == 'CFL_const') | (mode == 'normal'):
        mm_name = str('_CFL')+ str(CFL)
        m_name = str('CFL=')+ str(CFL)

    #===================================================
    # Reading continuous data
    #===================================================
    #Numerical:
    fname = dir1 + cont_num + str("_N") + Nelem \
            + mm_name + str("_Eps") + str(Epsilon)\
            + str('_')+ str(tt_) + str("t.dat")
    data = loadtxt(fname)
    x_cont = data[:, 0]
    u_cont = data[:, 1]
    del fname, data
    k_freq, u_amp, KE = compute_fft(u_cont[0:-1])
    E_tot = compute_Etotal(k_freq,KE)
    # computing the initial energy spectrum:
    fname = dir1 + cont_num + str("_N") + Nelem \
            + mm_name + str("_Eps") + str(Epsilon)\
            + str('_0.0000t.dat')
    data = loadtxt(fname)
    u_cont0 = data[:, 1]
    k_freq0, u_amp0, KE0 = compute_fft(u_cont0[0:-1])
    E_tot0 = compute_Etotal(k_freq0,KE0)
    G_num = u_amp/u_amp0

    #Exact:
    fname = dir1 + 'time_data/u_cont_exact_'+ str(tt_) + str("t.dat");
    data = loadtxt(fname);  # continuous exact nodal solution
    x_cont_exact = data[:, 0];
    u_cont_exact = data[:, 1];
    del fname, data
    k_freq_exact, u_amp_exact, KE_exact = compute_fft(u_cont_exact[0:-1])
    E_ex = compute_Etotal(k_freq_exact,KE_exact)
    # computing the initial energy spectrum:
    fname = dir1 + 'time_data/u_cont_exact_0.0000t.dat';
    data = loadtxt(fname);  # continuous exact nodal solution
    u_cont_exact0 = data[:, 1];
    del fname, data
    k_freq_exact0, u_amp_exact0, KE_exact0 = compute_fft(u_cont_exact0[0:-1])
    E_ex0 = compute_Etotal(k_freq_exact0,KE_exact0)
    G_ex = u_amp_exact/u_amp_exact0

    print('\nfft_u_amp_num: ',np.sqrt(E_tot), '\nfft_u_amp_ex: ',np.sqrt(E_ex))

    #===================================================
    # Reading discontinuous data
    #===================================================
    #Numerical:
    fname = dir1 + disc_num + str("_N") + Nelem \
            + mm_name + str("_Eps") + str(Epsilon)\
            + str('_')+ str(tt_) + str("t.dat")
    data = loadtxt(fname)
    x_disc = data[:, 0]
    u_disc = data[:, 1]
    N_disc = size(x_disc)
    Np = N_disc_ppt
    del fname, data
    #Exact:
    fname = dir1 + 'time_data/u_disc_exact_N'+Nelem+"_"+ str(tt_) + str("t.dat");
    data = loadtxt(fname);  # continuous exact nodal solution
    x_disc_exact = data[:, 0];
    u_disc_exact = data[:, 1];
    del fname, data

    print('u_cont_max: ',max(u_cont))
    print('u_disc_max: ',max(u_disc))
    print('u_cont_exact_max: ',max(u_cont_exact))
    print('u_disc_exact_max: ',max(u_disc_exact))
    if wave_form_=="Trigonometric":
        u_jump0 = -u_disc[0]
        u_jump1 = 0.5*(u_disc[Np] - u_disc[Np-1])
        u_jump = 0.5*(u_jump0+u_jump1)
        print('u_jump:     ',u_jump)
        print('est_amp:    ',max(u_disc)-u_jump)
        print('error: ',abs(1-max(u_disc)))
        print('error_cont: ',abs(max(u_cont_exact)-max(u_cont)))
        print('error_disc: ',abs(max(u_disc_exact)-max(u_disc)))
    print('N_disc:',N_disc,', Np:',Np)

    #============ Preparing data for Plotting ========================#
    L = abs(x_cont_exact[-1]-x_cont_exact[0])
    tau_p0 = float(gamma_)*float(tt_)*(Np1*nelem/L)**2
    tau_p0 = Decimal(tau_p0)
    tau_p0 =  Decimal(tau_p0.quantize(Decimal('.01')))

    if wave_form_=="Trigonometric":
        spectrum_fname = dir1 + 'tempdata/spectrum_numer_'+diffus_scheme+'p'+str(DG)+'RK'+RK+\
                         '_eta'+str(Epsilon)+'_N'+Nelem+mm_name+str('_')+\
                           str(tau_p0)+str('tau_pi_')+str(K_den)
        spectrum_fname_exact = dir1 + 'tempdata/spectrum_exact_'+diffus_scheme+'p'+str(DG)+'RK'+RK+\
                         '_eta'+str(Epsilon)+'_N'+Nelem+mm_name+str('_')+\
                           str(tau_p0)+str('tau_pi_')+str(K_den)
    elif wave_form_=="Gaussian":
        spectrum_fname = dir1 + 'tempdata/spectrum_numer_'+diffus_scheme+'p'+str(DG)+'RK'+RK+\
                         '_eta'+str(Epsilon)+'_N'+Nelem+mm_name+str('_')+\
                           str(tau_p0)+str('tau')
        spectrum_fname_exact = dir1 + 'tempdata/spectrum_exact_'+diffus_scheme+'p'+str(DG)+'RK'+RK+\
                         '_eta'+str(Epsilon)+'_N'+Nelem+mm_name+str('_')+\
                           str(tau_p0)+str('tau')

        u_amp = 0.5*u_amp # fix me there is a problem in fft_toolbox dealing with Gaussian
        KE = 0.25*KE
        E_tot = 0.25*E_tot
        u_amp_exact = 0.5 * u_amp_exact
        KE_exact = 0.25* KE_exact
        E_ex = 0.25*E_ex
        u_amp_exact0 = 0.5 * u_amp_exact0
        KE_exact0 = 0.25 * KE_exact0
        E_ex0 = 0.25*E_ex0
        u_amp0 = 0.5*u_amp0
        KE0 = 0.25 * KE0
        E_tot0 = 0.25*E_tot0
    else:
        print('unrecognized_waveform')

    # Saving FFT data:
    savez_spectrum_linearWave(k_freq,u_amp,KE,G_num,spectrum_fname)
    savez_spectrum_linearWave(k_freq_exact, u_amp_exact, KE_exact,G_ex\
                                             ,spectrum_fname_exact)
    np.savetxt(spectrum_fname+'.dat',np.transpose([k_freq,u_amp,KE,G_num]))
    np.savetxt(spectrum_fname_exact+'.dat'\
                 ,np.transpose([k_freq_exact, u_amp_exact, KE_exact,G_ex]))

    #=========================== PLOTTING Initial Solution ===============================#
    ylim_0 = list();
    ylim_1 = list();

    markers_on=[0 , 6 , 15, 24, Np-1]
    fig, (ax1,ax2) = plt.subplots(1,2)
    #plotting continuous data:
    ax1.plot(x_cont_exact, u_cont_exact,'-k',lw=0.7\
             ,label='Gaussian wave');
    #plotting discontinuous exact data:
    N_disc_ex = size(x_disc_exact)
    Np = N_disc_ppt
    label_disc_ex = 'Projected solution' ; # discontinuous exact label
    for i in range(0, N_disc_ex - Np, Np):
        xx = x_disc_exact[i:i + Np];
        uu = u_disc_exact[i:i + Np];
        if(int(i/Np)==15):
            ax1.plot(xx, uu, '--ob',lw=0.6,markevery=markers_on);
        else:
            ax1.plot(xx, uu, '--b',lw=1.5);
    i=i+Np
    xx = x_disc_exact[i:i + Np];
    uu = u_disc_exact[i:i + Np];
    ax1.plot(xx, uu, '--ob',markevery=[],label=label_disc_ex,lw=1.5);

    ylim_0.append(min(u_disc_exact));
    ylim_1.append(max(u_disc_exact));
    ylim_0.append(min(u_disc));
    ylim_1.append(max(u_disc));
    ylim_0.append(min(u_cont_exact));
    ylim_1.append(max(u_cont_exact));

    ax1.legend(fontsize=18,loc='upper right');
    ax1.set_title('Physical space',fontsize=20)
    ax1.set_xlabel(r'$x$', labelpad=3,fontsize=18);
    ax1.set_ylabel(r'$u(x,0)$', labelpad=3);
    n_divisions = 8;
    xtick_dx = (x_cont_exact[-1] - x_cont_exact[0] ) / n_divisions;
    xlabels = arange(x_cont_exact[0], x_cont_exact[-1]+xtick_dx,xtick_dx);
    xlocs = xlabels;
    plt.xticks(xlocs, xlabels);
    ax1.grid()
    ax1.set_xlim(min(x_cont_exact), max(x_cont_exact))
    ax1.set_ylim(min(ylim_0), max(ylim_1)*1.05)

    #================== Plot FFT ==========================#
    label_dg = diffus_scheme+'p'+ str(DG) + r'-$\eta$' + str(Epsilon) \
                +"_RK" + str(RK) # discontinuous label
    k_max_ = int(k_freq[-1]);
    if nelem%2==0:
        normalize_factor = Np1*nelem
    elif nelem%2==1:
        normalize_factor = Np1*nelem-1
    AA = 2*pi/(normalize_factor)  # to transform from an integer wavenumber to K = kh/(P+1)\
    #, non-dimensional one, multiplying by 2pi is because the integers k_freq is actually k_freq/2
    k_Nyquist = int(normalize_factor/2)

    #plotting:
    ax2.plot(AA*k_freq_exact,u_amp_exact,'-k',lw=1.0,label=r'Gaussian wave')
    ax2.plot(AA*k_freq,u_amp, '--b',lw=1.0,label=r'Projected solution')
    ax2.plot(-AA*k_freq_exact[::-1],u_amp_exact[::-1],'-k',lw=1.0)
    ax2.plot(-AA*k_freq[::-1],u_amp[::-1], '--b',lw=1.0)

    xlabels = [r'$-\pi$',r'$-3\pi/4$',r'$-\pi/2$',r'$-\pi/4$'\
              ,'0',r'$\pi/4$',r'$\pi/2$',r'$3\pi/4$',r'$\pi$'];
    xlocs = [-pi,-3*pi/4,-pi/2,-pi/4,0,pi/4,pi/2,3*pi/4,pi];
    plt.xticks(xlocs, xlabels);
    ax2.set_xlim(-pi,pi)
    ax2.set_ylim(0.001,0.008)
    ax2.legend();
    ax2.set_xlabel(r'$K$', labelpad=3);
    ax2.set_ylabel(r'$\left| \hat{u}(K;0) \right| $', labelpad=4);
    ax2.grid()
    ax2.set_title('Frequency space')

    fig.set_size_inches(16.0, 5.5, forward=True)
    fig.tight_layout(pad=0.0, w_pad=3.0, h_pad=10.0,rect=(0.0,0.0,1.0,1.0))

    figname_str='tempfig/gaussian_init_p'+str(DG)+'_N'\
               +Nelem+str('_0_00tau')
    figname = dir1+figname_str+'.png'
    plt.savefig(figname,format='png')
    figname = dir1+figname_str+'.eps'
    plt.savefig(figname,format='eps')

    plt.show()
    exit()

    return 'true'

def plot_diffus(dir1, mode, DG, RK, CFL, Nelem, N_disc_ppt, tt_, dt_ \
                , Epsilon, gamma_, diffus_scheme, cont_num, disc_num, \
                T, K_den,wave_form_):
    nelem = int(Nelem)
    Np1 = int(int(DG)+1)
    dt = float(dt_)
    
    if  (mode == 'test') | (mode == 'dt_const'):
        mm_name = str('_dt') + dt_
        m_name = str('dt=') + dt_
    elif(mode == 'CFL_const') | (mode == 'normal'):
        mm_name = str('_CFL')+ str(CFL)     
        m_name = str('CFL=')+ str(CFL)  
        
    #===================================================
    # Reading continuous data
    #===================================================
    #Numerical:
    fname = dir1 + cont_num + str("_N") + Nelem \
            + mm_name + str("_Eps") + str(Epsilon)\
            + str('_')+ str(tt_) + str("t.dat")
    data = loadtxt(fname)
    x_cont = data[:, 0]
    u_cont = data[:, 1]
    del fname, data
    k_freq, u_amp, KE = compute_fft(u_cont[0:-1])
    E_tot = compute_Etotal(k_freq,KE)
    # computing the initial energy spectrum:
    fname = dir1 + cont_num + str("_N") + Nelem \
            + mm_name + str("_Eps") + str(Epsilon)\
            + str('_0.0000t.dat')
    data = loadtxt(fname)
    u_cont0 = data[:, 1]
    k_freq0, u_amp0, KE0 = compute_fft(u_cont0[0:-1])
    E_tot0 = compute_Etotal(k_freq0,KE0)
    G_num = u_amp/u_amp0
    
    #Exact:
    fname = dir1 + 'time_data/u_cont_exact_'+ str(tt_) + str("t.dat");
    data = loadtxt(fname);  # continuous exact nodal solution
    x_cont_exact = data[:, 0];
    u_cont_exact = data[:, 1];
    del fname, data
    k_freq_exact, u_amp_exact, KE_exact = compute_fft(u_cont_exact[0:-1])
    E_ex = compute_Etotal(k_freq_exact,KE_exact)
    # computing the initial energy spectrum:
    fname = dir1 + 'time_data/u_cont_exact_0.0000t.dat';
    data = loadtxt(fname);  # continuous exact nodal solution
    u_cont_exact0 = data[:, 1];
    del fname, data
    k_freq_exact0, u_amp_exact0, KE_exact0 = compute_fft(u_cont_exact0[0:-1])
    E_ex0 = compute_Etotal(k_freq_exact0,KE_exact0)
    G_ex = u_amp_exact/u_amp_exact0
    
    print('\nfft_u_amp_num: ',np.sqrt(E_tot), '\nfft_u_amp_ex: ',np.sqrt(E_ex))

    #===================================================
    # Reading discontinuous data
    #===================================================
    #Numerical:
    fname = dir1 + disc_num + str("_N") + Nelem \
            + mm_name + str("_Eps") + str(Epsilon)\
            + str('_')+ str(tt_) + str("t.dat")
    data = loadtxt(fname)
    x_disc = data[:, 0]
    u_disc = data[:, 1]
    N_disc = size(x_disc)
    Np = N_disc_ppt
    del fname, data
    #Exact:
    fname = dir1 + 'time_data/u_disc_exact_N'+Nelem+"_"+ str(tt_) + str("t.dat");
    data = loadtxt(fname);  # continuous exact nodal solution
    x_disc_exact = data[:, 0];
    u_disc_exact = data[:, 1];
    del fname, data
    
    print('u_cont_max: ',max(u_cont))
    print('u_disc_max: ',max(u_disc))
    print('u_cont_exact_max: ',max(u_cont_exact))
    print('u_disc_exact_max: ',max(u_disc_exact))
    if wave_form_=="Trigonometric":
        u_jump0 = -u_disc[0]
        u_jump1 = 0.5*(u_disc[Np] - u_disc[Np-1])
        u_jump = 0.5*(u_jump0+u_jump1) 
        print('u_jump:     ',u_jump)
        print('est_amp:    ',max(u_disc)-u_jump)
        print('error: ',abs(1-max(u_disc)))
        print('error_cont: ',abs(max(u_cont_exact)-max(u_cont)))
        print('error_disc: ',abs(max(u_disc_exact)-max(u_disc)))
    print('N_disc:',N_disc,', Np:',Np)

    #============ Preparing data for Plotting ========================#
    L = abs(x_cont_exact[-1]-x_cont_exact[0])
    tau_p0 = float(gamma_)*float(tt_)*(Np1*nelem/L)**2
    tau_p0 = Decimal(tau_p0)
    tau_p0 =  Decimal(tau_p0.quantize(Decimal('.01')))
    
    if wave_form_=="Trigonometric":   
        spectrum_fname = dir1 + 'tempdata/spectrum_numer_'+diffus_scheme+'p'+str(DG)+'RK'+RK+\
                         '_eta'+str(Epsilon)+'_N'+Nelem+mm_name+str('_')+\
                           str(tau_p0)+str('tau_pi_')+str(K_den)
        spectrum_fname_exact = dir1 + 'tempdata/spectrum_exact_'+diffus_scheme+'p'+str(DG)+'RK'+RK+\
                         '_eta'+str(Epsilon)+'_N'+Nelem+mm_name+str('_')+\
                           str(tau_p0)+str('tau_pi_')+str(K_den)  
                           
        casename0 = diffus_scheme+'p'+str(DG)+'RK'+RK+\
                    '_eta'+str(Epsilon)+'_N'+Nelem+mm_name+str('_')+\
                    str(tau_p0)+str('tau_pi_')+str(K_den)
        casename1 = diffus_scheme+'p'+str(DG)+'RK'+RK+\
                    '_eta'+str(Epsilon)+'_N'+Nelem+mm_name+str('_')+\
                    str(T)+str('T_pi_')+str(K_den) 
        fig_casetitle_name0 = str(r'K$=\pi/')+str(K_den)+'$'+str(r', $\tau_{p}=$') + str(tau_p0)
        fig_casetitle_name1 = str(r'K$=\pi/')+str(K_den)+'$'+str(', T=') + str(T)
    elif wave_form_=="Gaussian":
        spectrum_fname = dir1 + 'tempdata/spectrum_numer_'+diffus_scheme+'p'+str(DG)+'RK'+RK+\
                         '_eta'+str(Epsilon)+'_N'+Nelem+mm_name+str('_')+\
                           str(tau_p0)+str('tau')
        spectrum_fname_exact = dir1 + 'tempdata/spectrum_exact_'+diffus_scheme+'p'+str(DG)+'RK'+RK+\
                         '_eta'+str(Epsilon)+'_N'+Nelem+mm_name+str('_')+\
                           str(tau_p0)+str('tau')
        casename0 = diffus_scheme+'p'+str(DG)+'RK'+RK+\
                    '_eta'+str(Epsilon)+'_N'+Nelem+mm_name+str('_')+\
                    str(tau_p0)+str('tau')
        casename1 = diffus_scheme+'p'+str(DG)+'RK'+RK+\
                    '_eta'+str(Epsilon)+'_N'+Nelem+mm_name+str('_')+\
                    str(T)+str('T')
        fig_casetitle_name0 = str(r'$\tau_{p}=$') + str(tau_p0) 
        fig_casetitle_name1 = str('T=') + str(T)  
        
        u_amp = 0.5* u_amp # fix me there is a problem in fft_toolbox dealing with Gaussian
        KE = 0.25*KE 
        E_tot = 0.25*E_tot
        u_amp_exact = 0.5 * u_amp_exact
        KE_exact = 0.25* KE_exact
        E_ex = 0.25*E_ex
        u_amp_exact0 = 0.5 * u_amp_exact0
        KE_exact0 = 0.25 * KE_exact0
        E_ex0 = 0.25*E_ex0
        u_amp0 = 0.5*u_amp0
        KE0 = 0.25 * KE0
        E_tot0 = 0.25*E_tot0
    else:
        print('unrecognized_waveform')  
    
    # Saving FFT data:
    savez_spectrum_linearWave(k_freq,u_amp,KE,G_num,spectrum_fname) 
    # Saving FFT data:
    savez_spectrum_linearWave(k_freq_exact, u_amp_exact, KE_exact,G_ex,spectrum_fname_exact) 

    np.savetxt(spectrum_fname+'.dat',np.transpose([k_freq,u_amp,KE,G_num]))
    np.savetxt(spectrum_fname_exact+'.dat'\
                 ,np.transpose([k_freq_exact, u_amp_exact, KE_exact,G_ex]))
    

    #=========================== PLOTTING Initial Solution ===============================#
    ylim_0 = list();
    ylim_1 = list();
    fig = plt.figure();
    #plotting continuous data:
    plt.plot(x_cont_exact, u_cont_exact,'-k',label='Continuous exact solution');

    #plotting discontinuous exact data:
    N_disc_ex = size(x_disc_exact)
    Np = N_disc_ppt
    label_disc_ex = "Projected exact solution" ; # discontinuous exact label
    for i in range(0, N_disc_ex - Np, Np):
        xx = x_disc_exact[i:i + Np];
        uu = u_disc_exact[i:i + Np];
        if(int(i/Np)==15):
            pyplot.plot(xx, uu, '--ob',lw=0.6);
        else:
            pyplot.plot(xx, uu, '--b',lw=0.6);
    i=i+Np
    xx = x_disc_exact[i:i + Np];
    uu = u_disc_exact[i:i + Np];
    pyplot.plot(xx, uu, '--b', label=label_disc_ex,lw=0.6);

    ylim_0.append(min(u_disc_exact));
    ylim_1.append(max(u_disc_exact));
    ylim_0.append(min(u_disc));
    ylim_1.append(max(u_disc));
    ylim_0.append(min(u_cont_exact));
    ylim_1.append(max(u_cont_exact));

    plt.legend();
    plt.title(fig_casetitle_name0)
    plt.xlabel('x', labelpad=10);
    plt.ylabel('u(x)', labelpad=10);
    n_divisions = 8;
    xtick_dx = (x_cont_exact[-1] - x_cont_exact[0] ) / n_divisions;
    xlabels = arange(x_cont_exact[0], x_cont_exact[-1]+xtick_dx,xtick_dx);
    xlocs = xlabels;
    plt.xticks(xlocs, xlabels);
    plt.grid()
    plt.xlim(min(x_cont_exact), max(x_cont_exact))
    plt.ylim(min(ylim_0), max(ylim_1)*1.3)
    fig.set_size_inches(12.0, 10.0, forward=True)

    fig.tight_layout()
    figname = dir1 + 'tempfig/u_init_'+casename0+'.png'
    plt.savefig(figname,format='png')
    figname = dir1 + 'tempfig/u_init_'+casename0+'.eps'
    plt.savefig(figname,format='eps')
    plt.show()

    #=========================== PLOTTING Solution(1) ============================#
    ylim_0 = list();
    ylim_1 = list();
    fig = plt.figure();
    #plotting continuous data:
    plt.plot(x_cont_exact, u_cont_exact, '-k', label='Exact solution');
    label_cont = diffus_scheme+'p'+ str(DG) + r'-$\eta$' + str(Epsilon) \
                +"_RK" + str(RK) +", CFL="+str(CFL) # discontinuous label
    plt.plot(x_cont, u_cont, '--m', label=label_cont,lw=1.2);
    
    #plotting discontinuous numerical data:
    markers_on=[0,int(Np/6),int(Np/4),int(Np/2),int(Np/2)+int(Np/4), int(Np/2)+int(Np/6), Np-1]
    #markers_on=1;
    for i in range(0, N_disc-Np, Np): # plot all the elements curves except the last one to be able to define a label
        xx = x_disc[i:i + Np];
        uu = u_disc[i:i + Np];
        plt.plot(xx, uu, 'om',lw=0.7,markevery=markers_on)           
        
    i=i+Np
    xx = x_disc[i:i+Np];
    uu = u_disc[i:i+Np];
    plt.plot(xx, uu, 'om', markevery=markers_on);
    
    ylim_0.append(min(u_cont_exact));
    ylim_1.append(max(u_cont_exact));
    ylim_0.append(min(u_disc));
    ylim_1.append(max(u_disc));

    plt.legend();
    plt.title(fig_casetitle_name0)    
    plt.xlabel('x', labelpad=10);
    plt.ylabel('u(x)', labelpad=10);
    n_divisions = 8;
    xtick_dx = (x_cont_exact[-1] - x_cont_exact[0] ) / n_divisions;
    xlabels = arange(x_cont_exact[0], x_cont_exact[-1]+xtick_dx,xtick_dx);
    xlocs = xlabels;
    plt.xticks(xlocs, xlabels);
    plt.grid()
    plt.xlim(min(x_cont_exact), max(x_cont_exact))
    plt.ylim(min(ylim_0), max(ylim_1)*1.3)

    fig.tight_layout()
    figname = dir1 + 'tempfig/u_cont_'+casename0+'.png'
    fig.set_size_inches(12.0, 10.0, forward=True)
    plt.savefig(figname)
    
    #=========================== Finding the best fitted sine wave ==================# 
    #if wave_form_=="Trigonometric":
    #    x_nonfit, u_nonfit, x_fitted, u_fitted, u_first_guess = compute_fitted_wave(Np,int(Nelem),x_disc,u_disc)
       
    #    fig = pyplot.figure();
    #    pyplot.plot(x_cont_exact, u_cont_exact, '-k',label='exact solution')
    #    pyplot.plot(x_nonfit, u_nonfit, '.m')
    #    pyplot.plot(x_fit, u_first_guess,'--', label='first guess')
    #    pyplot.plot(x_fit, u_fitted, '-.b',label='after fitting')
    #    pyplot.legend()
    #    pyplot.show()
    
    #=========================== PLOTTING Solution(2) ============================#
    fig = pyplot.figure();
    
    ylim_0 = list();
    ylim_1 = list();
    ylim_0.append(min(u_disc_exact));
    ylim_1.append(max(u_disc_exact));
    ylim_0.append(min(u_disc));
    ylim_1.append(max(u_disc));

    #plotting discontinuous numerical data:
    label_disc = diffus_scheme+'p'+ str(DG) + r'-$\eta$' + str(Epsilon) \
    +"_RK" + str(RK) +", CFL="+str(CFL) # discontinuous label
    for i in range(0, N_disc-Np, Np):
        xx = x_disc[i:i + Np];
        uu = u_disc[i:i + Np];
        pyplot.plot(xx, uu, '.m');
    i=i+Np
    xx = x_disc[i:i + Np];
    uu = u_disc[i:i + Np];
    pyplot.plot(xx, uu, '.m', label=label_disc);
    
    #plotting discontinuous exact data:
    N_disc_ex = size(x_disc_exact)
    Np = N_disc_ppt
    label_disc_ex = "exact discontinuous" ; # discontinuous exact label
    for i in range(0, N_disc_ex - Np, Np):
        xx = x_disc_exact[i:i + Np];
        uu = u_disc_exact[i:i + Np];
        pyplot.plot(xx, uu, '-k',lw=0.6);
    i=i+Np
    xx = x_disc_exact[i:i + Np];
    uu = u_disc_exact[i:i + Np];
    pyplot.plot(xx, uu, '-k', label=label_disc_ex,lw=0.6);
    #pyplot.plot(x_disc,u_disc,'.b',label=label_disc,\
    #    markerfacecolor='None',markevery=10,markersize=10)
    pyplot.legend();
    
    pyplot.title(str("Plotting of discontinuous solutions, ")+fig_casetitle_name0)
    pyplot.xlabel('x', labelpad=10);
    pyplot.ylabel('u(x)', labelpad=10);

    pyplot.xlim(min(x_cont_exact), max(x_cont_exact));
    pyplot.ylim(min(ylim_0) * 1.05, max(ylim_1) * 1.05);
    #pyplot.ylim(-.14, 0.14);

    n_divisions = 8;
    xtick_dx = (x_cont_exact[-1] - x_cont_exact[0] ) / n_divisions;
    xlabels = arange(x_cont_exact[0], x_cont_exact[-1]+xtick_dx,xtick_dx);
    #xlabels = [0,10,20,30,40,50,60,70,80];
    xlocs = xlabels;
    pyplot.xticks(xlocs, xlabels);
    pyplot.grid()
    pyplot.xlim(min(x_cont_exact), max(x_cont_exact))
    pyplot.ylim(min(ylim_0), max(ylim_1)*1.3)

    fig.tight_layout()
    figname = dir1 + 'tempfig/u_disc_'+casename0+'.png'
    fig.set_size_inches(12.0, 10.0, forward=True)
    pyplot.savefig(figname)
    
    #================== Plot FFT ==========================#
    label_dg = diffus_scheme+'p'+ str(DG) + r'-$\eta$' + str(Epsilon) \
                +"_RK" + str(RK) # discontinuous label
    k_max_ = int(k_freq[-1]);
    if nelem%2==0:
        normalize_factor = Np1*nelem
    elif nelem%2==1:
        normalize_factor = Np1*nelem-1
    AA = 2*pi/(normalize_factor)  # to transform from an integer wavenumber to K = kh/(P+1)\
    #, non-dimensional one, multiplying by 2pi is because the integers k_freq is actually k_freq/2
    k_Nyquist = int(normalize_factor/2)
    
    # Plot(1) vs k integer
    fig, ax = plt.subplots()
    plt.plot(k_freq0, u_amp0, '-.^c',markevery=1, \
        label=r'Initial, E$_{tot}$='+str(np.round(E_tot0,8)))
    plt.plot(k_freq_exact, u_amp_exact, '-k',markevery=1\
        , label=r'Exact, E$_{tot}$= '+str(np.round(E_ex,8)))
    plt.plot(k_freq, u_amp, ':om',markevery=1, \
        label=label_dg+r', E$_{tot}$='+str(np.round(E_tot,8)))
    
    #xlabels = ['0',r'$\pi$/6',r'$\pi$/4',r'$\pi$/3',r'$\pi$/2',r'3$\pi$/4',r'$\pi$'];
    #xlocs = arange
    #plt.xticks(xlocs, xlabels);
    
    plt.xlim(0,k_Nyquist)
    plt.legend();
    plt.xlabel('k', labelpad=2);
    plt.ylabel('|u|', labelpad=2);
    plt.grid()
    plt.title(str(r'$\tau_{p}=') + str(tau_p0) + str('$'))

    fig.set_size_inches(10.0, 8.0, forward=True)
    fig.tight_layout(pad=0.0, w_pad=10.0, h_pad=10.0)
    
    #figname = dir1 + 'tempfig/fftint_'+casename0+'.png'
    #fig.set_size_inches(15.0, 9.0, forward=True)
    #plt.savefig(figname,format='png')
    #figname = dir1 + 'tempfig/fft_'+casename0+'.eps'
    #plt.savefig(figname,format='eps')
    
    # Plot(2) vs K=kh/(P+1)
    fig, ax = plt.subplots()
    plt.plot(AA*k_freq_exact, u_amp_exact, '-k',markevery=1\
        , label=r'Exact, E$_{tot}$= '+str(np.round(E_ex,8)))
    plt.plot(AA*k_freq, u_amp, ':om',markevery=1, \
        label=label_dg+r', E$_{tot}$='+str(np.round(E_tot,8)))
    
    xlabels = ['0',r'$\pi$/6',r'$\pi$/4',r'$\pi$/3',r'$\pi$/2'\
    ,r'3$\pi$/4',r'$\pi$'];
    xlocs = [0,pi/6,pi/4,pi/3,pi/2,3*pi/4,pi];
    plt.xticks(xlocs, xlabels);
    plt.xlim(0,pi)
    plt.legend();
    plt.xlabel('K', labelpad=2);
    plt.ylabel('|u|', labelpad=2);
    plt.grid()
    plt.title(str(r'$\tau_{p}=') + str(tau_p0) + str('$'))

    fig.set_size_inches(10.0, 8.0, forward=True)
    fig.tight_layout(pad=0.0, w_pad=10.0, h_pad=10.0)
    
    figname = dir1 + 'tempfig/fft_'+casename0+'.png'
    fig.set_size_inches(10.0, 8.0, forward=True)
    plt.savefig(figname,format='png')
    figname = dir1 + 'tempfig/fft_'+casename0+'.eps'
    plt.savefig(figname,format='eps')
    
    # Plot(3) , G
    fig, ax = plt.subplots()
    plt.plot(AA*k_freq_exact, G_ex\
        ,'-k',markevery=1, label=r'Exact, $e^{-K^2 \tau_{p}}$')
    plt.plot(AA*k_freq, G_num\
        ,':om',markevery=1,label=label_dg)
    
    xlabels = ['0',r'$\pi$/4',r'$\pi$/2',r'3$\pi$/4',r'$\pi$'];
    xlocs = [0,pi/4,pi/2,3*pi/4,pi];
    plt.xticks(xlocs, xlabels);
    plt.xlim(0,pi)
    if tau_p0 < 0.1:
        plt.ylim(0.85,1.01)
    else:
        plt.ylim(0.0,1.1)
    plt.legend();
    plt.xlabel('K', labelpad=2);
    plt.ylabel(r'G(K;$\tau_{p}$)', labelpad=2);
    plt.grid()
    plt.title(str(r'$\tau_{p}=') + str(tau_p0)+str('$'))

    fig.set_size_inches(8.0, 6.0, forward=True)
    fig.tight_layout(pad=0.0, w_pad=10.0, h_pad=10.0)
    
    figname = dir1 + 'tempfig/G_'+casename0+'.png'
    fig.set_size_inches(8.0, 6.0, forward=True)
    plt.savefig(figname,format='png')
    figname = dir1 + 'tempfig/G_'+casename0+'.eps'
    plt.savefig(figname,format='eps')

    print('Kpi : ',AA*k_freq[k_Nyquist])
    print('G_pi: ',G_num[k_Nyquist])
    
    #=========================== PLOTTING Errors Evolution ============================#
    # Read error data:
    res_dir = './Results/'
    fname = dir1 + 'errors/errors_N'+str(Nelem)+'_CFL'+str(CFL)+'_Eps'+str(Epsilon)+'_'+str(T)+'T.dat'
    data = loadtxt(fname);  # continuous exact nodal solution
    time  = data[:, 0];
    if wave_form_=="Trigonometric":
        tau_p = float(gamma_)*time*(Np1*nelem/L)**2
    else:
        tau_p = time; 
    L1err = data[:, 1];
    L2err = data[:, 2];
    L1err_nodal = data[:, 3];
    L2err_nodal = data[:, 4];

    #plotting
    fig = pyplot.figure();
    pyplot.plot(tau_p,L1err,'-.sb',label=r'$L_{1}(u(\xi))$')
    pyplot.plot(tau_p, L2err,'-ok',label=r'$L_{2}(u(\xi))$')
    pyplot.plot(tau_p,L1err_nodal,':vm',label=r'$L_{1}$, nodal')
    pyplot.plot(tau_p, L2err_nodal,'--c',label=r'$L_{2}$, nodal')
    
    #plt.xlim(0,2.0)
    #plt.ylim(0,0.025)
    
    pyplot.xlabel('time')
    pyplot.ylabel('errors')
    pyplot.legend();
    plt.title(fig_casetitle_name1)
    
    figname = dir1 + 'tempfig/errors_'+casename1+'.png'
    fig.set_size_inches(10.0, 8.0, forward=True)
    pyplot.savefig(figname)
    
    #------------------------
    pyplot.show()
    exit()
    return 'true'

def plot_error_vs_time(K_den):
    T=0.347;
    mode='dt_const'
    eta=Decimal(1.00)
    eta_br1 = '1.33'
    eta_br2 = '2.00'
    eta_ldg = '0.00'
    P=2
    RK=3
    Nelem=8
    N_disc_ppt=7
    CFL=0.0001
    CFL_str='CFL0.0001'
    dt_='1.000e-4'
    dt = float(dt_);
    gamma=0.01;


    if  (mode == 'test') | (mode == 'dt_const'):
        mm_name = str('_dt') + dt_
        m_name = str('dt=') + dt_
    elif(mode == 'CFL_const') | (mode == 'normal'):
        mm_name = str('_CFL')+ str(CFL)
        m_name = str('CFL=')+ str(CFL)

    # Read BR1 data:
    res_dir = './BR1p2RK3/'
    fname = res_dir + 'errors/errors_N'+str(Nelem)+'_'+CFL_str+'_Eps'+str(eta_br1)+'_'+str(T)+'T.dat'
    data = loadtxt(fname);  # continuous exact nodal solution
    time_br1  = data[:, 0];
    L2err_br1 = data[:, 2];
    del fname, data

    # Read BR2 data:
    res_dir = './BR2p2RK3/'
    fname = res_dir + 'errors/errors_N'+str(Nelem)+'_'+CFL_str+'_Eps'+str(eta_br2)+'_'+str(T)+'T.dat'
    data = loadtxt(fname);  # continuous exact nodal solution
    time_br2  = data[:, 0];
    L2err_br2 = data[:, 2];
    del fname, data

    # Read LDG data:
    res_dir = './LDGp2RK3/'
    fname = res_dir + 'errors/errors_N'+str(Nelem)+'_'+CFL_str+'_Eps'+str(eta_ldg)+'_'+str(T)+'T.dat'
    data = loadtxt(fname);  # continuous exact nodal solution
    time_ldg = data[:, 0];
    L2err_ldg = data[:, 2];
    del fname, data

    tau_p = float(gamma)*time_ldg*((int(P)+1)*int(Nelem))**2

    fig, ax = pyplot.subplots()

    pyplot.plot(tau_p, L2err_br1, '--b', label=r'BR1p$'+ str(P) + '$' + r'-$\eta'+str(eta_br1)+'$' ,lw=0.7);
    pyplot.plot(tau_p, L2err_br2, '-c', label=r'BR2p$'+ str(P) + '$' + r'-$\eta'+str(eta_br2)+'$' ,lw=0.7);
    pyplot.plot(tau_p, L2err_ldg, '-.m', label=r'LDGp$'+ str(P) + '$' + r'-$\eta'+str(eta_ldg)+'$' ,lw=0.7);
    #pyplot.plot(time_fdc/2, L2err_fdc, ':om',markevery=5, label=r'FD$6$-central',linewidth=0.70,markersize=8.0)

    #pyplot.legend();
    ax.legend(loc='lower right', edgecolor='black',fontsize=20)
    ax.spines['left'].set_visible(True)
    ax.spines['bottom'].set_visible(True)
    ax.spines['right'].set_visible(True)
    ax.spines['top'].set_visible(True)


    pyplot.xlabel(r'$\tau_{p}$', labelpad=2.5,size=23);
    pyplot.ylabel(r'L$_{2}$ error', labelpad=2.5,size=23);
    pyplot.title(str(r'K$=\pi/')+str(K_den)+'$')

    pyplot.xlim(0.0, 2.0)

    fig.set_size_inches(10, 7.0, forward=True)
    fig.tight_layout(pad=0.0, w_pad=10.0, h_pad=10.0,rect=None)
    #pyplot.tight_layout()
    figname = './figs/' + 'errors_compare_Kpi_'+str(K_den)+'_N'+str(Nelem)\
    +'_sdisc_dt1.0e-4.png'
    #fig.set_size_inches(16.0, 2.0, forward=True)
    pyplot.savefig(figname,format='png')
    figname = './figs/' + 'errors_compare_Kpi_'+str(K_den)+'_N'+str(Nelem)\
    +'_sdisc_dt1.0e-4.eps'
    pyplot.savefig(figname,format='eps')

    pyplot.show()


def compare_diffus_schemes(dir1, mode, DG, RK, CFL, Nelem, N_disc_ppt, tt_, dt_ \
                , Epsilon, gamma_, diffus_scheme, cont_num, disc_num, \
                T, K_den,wave_form_):
    L=2.0;            
    nelem = int(Nelem)
    Np1 = int(int(DG)+1)
    dt = float(dt_)
    if  (mode == 'test') | (mode == 'dt_const'):
        mm_name = str('_dt') + dt_
        m_name = str('dt=') + dt_
    elif(mode == 'CFL_const') | (mode == 'normal'):
        mm_name = str('_CFL')+ str(CFL)     
        m_name = str('CFL=')+ str(CFL)  
    
    tau_p0 = float(gamma_)*float(tt_)*(Np1*nelem/L)**2
    tau_p0 = Decimal(tau_p0)
    tau_p0 =  Decimal(tau_p0.quantize(Decimal('.01')))
    
    if wave_form_=="Trigonometric":   
        spectrum_fname = dir1 + 'tempdata/spectrum_numer_'+diffus_scheme+'p'+str(DG)+'RK'+RK+\
                         '_eta'+str(Epsilon)+'_N'+Nelem+mm_name+str('_')+\
                           str(tau_p0)+str('tau_pi_')+str(K_den)
        spectrum_fname_exact = dir1 + 'tempdata/spectrum_exact_'+diffus_scheme+'p'+str(DG)+'RK'+RK+\
                         '_eta'+str(Epsilon)+'_N'+Nelem+mm_name+str('_')+\
                           str(tau_p0)+str('tau_pi_')+str(K_den)  
                           
        casename0 = diffus_scheme+'p'+str(DG)+'RK'+RK+\
                    '_eta'+str(Epsilon)+'_N'+Nelem+mm_name+str('_')+\
                    str(tau_p0)+str('tau_pi_')+str(K_den) 
        casename1 = diffus_scheme+'p'+str(DG)+'RK'+RK+\
                    '_eta'+str(Epsilon)+'_N'+Nelem+mm_name+str('_')+\
                    str(T)+str('T_pi_')+str(K_den) 
        fig_casetitle_name0 = str(r'K$=\pi/')+str(K_den)+'$'+str(', $\tau_{p}=$') + str(tau_p0)
        fig_casetitle_name1 = str(r'K$=\pi/')+str(K_den)+'$'+str(', T=') + str(T) 
    elif wave_form_=="Gaussian":
        spectrum_fname = dir1 + 'tempdata/spectrum_numer_'+diffus_scheme+'p'+str(DG)+'RK'+RK+\
                         '_eta'+str(Epsilon)+'_N'+Nelem+mm_name+str('_')+\
                           str(tau_p0)+str('tau')
        spectrum_fname_exact = dir1 + 'tempdata/spectrum_exact_'+diffus_scheme+'p'+str(DG)+'RK'+RK+\
                         '_eta'+str(Epsilon)+'_N'+Nelem+mm_name+str('_')+\
                           str(tau_p0)+str('tau')
        casename0 = diffus_scheme+'p'+str(DG)+'RK'+RK+\
                    '_eta'+str(Epsilon)+'_N'+Nelem+mm_name+str('_')+\
                    str(tau_p0)+str('tau')
        casename1 = diffus_scheme+'p'+str(DG)+'RK'+RK+\
                    '_eta'+str(Epsilon)+'_N'+Nelem+mm_name+str('_')+\
                    str(T)+str('T')
        fig_casetitle_name0 = str(r'$\tau_{p}=$') + str(tau_p0) 
        fig_casetitle_name1 = str('T=') + str(T)  
    else:
        print('unrecognized_waveform')  
    
    # Laoding FFT data:
    Epsilon = [0.67, 0.80, 1.00, 1.50, 2.00];
    G_num = [[],[],[],[],[]]
    for i in range(0,size(Epsilon)):
        eps = Decimal(Epsilon[i])
        eps =  Decimal(eps.quantize(Decimal('.01')))
        spectrum_fname = dir1 + 'tempdata/spectrum_'+diffus_scheme+'p'+str(DG)+'RK'+RK+\
                             '_eta'+str(eps)+'_N'+Nelem+mm_name+str('_')+\
                               str(tau_p0)+str('tau')
        k_freq,u_amp,KE,G_ = loadz_spectrum_linearWave(spectrum_fname) 
        G_num[i]=G_
    k_freq_exact, u_amp_exact, KE_exact,G_ex = loadz_spectrum_linearWave(spectrum_fname_exact)
    
    if nelem%2==0:
        normalize_factor = Np1*nelem
    elif nelem%2==1:
        normalize_factor = Np1*nelem-1
    AA = 2*pi/(normalize_factor)
    # Plot(3) , G
    fig, ax = plt.subplots()
    plt.plot(AA*k_freq[:-2], G_num[0]\
        ,'-.ok',markerfacecolor="None",markevery=3,label=r'$\eta='+str(Epsilon[0])+str('$'))
    plt.plot(AA*k_freq[:-2], G_num[1]\
        ,'->b',markerfacecolor="None",markevery=8,label=r'$\eta='+str(Epsilon[1])+str('$'))
    plt.plot(AA*k_freq[:-2], G_num[2]\
        ,'--sc',markerfacecolor="None",markevery=5,label=r'$\eta='+str(Epsilon[2])+str('$'))
    plt.plot(AA*k_freq[:-2], G_num[3]\
        ,'-.*g',markerfacecolor="None",markevery=8,label=r'$\eta='+str(Epsilon[3])+str('$'))
    plt.plot(AA*k_freq[:-2], G_num[4]\
        ,':dm',markerfacecolor="None",markevery=8,label=r'$\eta='+str(Epsilon[4])+str('$'))
    plt.plot(AA*k_freq_exact[:-2], G_ex\
        ,'-k',markevery=1, label=r'Exact, $e^{-K^2 \tau_{p}}$')
    
    xlabels = ['0',r'$\pi$/4',r'$\pi$/2',r'3$\pi$/4',r'$\pi$'];
    xlocs = [0,pi/4,pi/2,3*pi/4,pi];
    plt.xticks(xlocs, xlabels);
    plt.xlim(0,pi)
    ylabels = ['0.85','0.90','0.95','1.00'];
    ylocs = [0.85,0.9,0.95,1.00];
    plt.yticks(ylocs, ylabels);
    if tau_p0 < 0.1:
        plt.ylim(0.85,1.01)
    else:
        plt.ylim(0.0,1.1)
    plt.legend();
    plt.xlabel('K', labelpad=2);
    plt.ylabel(r'G(K;$\tau_{p}$)', labelpad=2);
    plt.grid()
    plt.title(str(r'BR2p2, $\tau_{p}=') + str(tau_p0)+str('$'))

    fig.set_size_inches(8.5, 5.5, forward=True)
    fig.tight_layout(pad=0.0, w_pad=10.0, h_pad=10.0)
    
    figname = dir1 + 'tempfig/G_compareEta_'+casename0+'.png'
    fig.set_size_inches(8.5, 5.5, forward=True)
    plt.savefig(figname,format='png')
    figname = dir1 + 'tempfig/G_compareEta_'+casename0+'.eps'
    plt.savefig(figname,format='eps')
    
    plt.show()
    
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#                                  A D V E C T I O N
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

def plot_advec(dir1, mode, DG, RK, CFL, Nelem, N_disc_ppt, tt_, dt_ \
                     , Beta, Epsilon, cont_num, disc_num, T_period, K_den):

    dt = float(dt_);
    
    if  (mode == 'test') | (mode == 'dt_const'):
        mm_name = str('_dt') + dt_
        m_name = str('dt=') + dt_
    elif(mode == 'CFL_const') | (mode == 'normal'):
        mm_name = str('_CFL')+ str(CFL)     
        m_name = str('CFL=')+ str(CFL)             

    #===================================================
    # Reading continuous data
    #===================================================
    #Numerical:
    fname = dir1 + cont_num + str("_N") + Nelem \
            + mm_name + str("_Beta") + str(Beta)\
            + str('_')+ str(tt_) + str("t.dat")
    data = loadtxt(fname)
    x_cont = data[:, 0]
    u_cont = data[:, 1]
    del fname, data
    k_freq, u_amp, KE = compute_fft(u_cont)
    E_tot = compute_Etotal(k_freq,KE)
    #Exact:
    fname = dir1 + 'time_data/u_cont_exact_'+ str(tt_) + str("t.dat");
    data = loadtxt(fname);  # continuous exact nodal solution
    x_cont_exact = data[:, 0];
    u_cont_exact = data[:, 1];
    del fname, data
    k_freq_exact, u_amp_exact, KE_exact = compute_fft(u_cont_exact)
    E_ex = compute_Etotal(k_freq_exact,KE_exact)
    
    print('\nu_amp_num: ',np.sqrt(E_tot), '\tu_amp_ex: ',np.sqrt(E_ex), '\n')
    #===================================================
    # Reading discontinuous data
    #===================================================
    #Numerical:
    fname = dir1 + disc_num + str("_N") + Nelem \
            + mm_name + str("_Beta") + str(Beta)\
            + str('_')+ str(tt_) + str("t.dat")
    data = loadtxt(fname)
    x_disc = data[:, 0]
    u_disc = data[:, 1]
    nn = size(x_disc)
    Np = N_disc_ppt
    del fname, data
    #Exact:
    fname = dir1 + 'time_data/u_disc_exact_N'+Nelem+"_"+ str(tt_) + str("t.dat");
    data = loadtxt(fname);  # continuous exact nodal solution
    x_disc_exact = data[:, 0];
    u_disc_exact = data[:, 1];
    del fname, data
    
    u_jump0 = -u_disc[0]
    u_jump1 = 0.5*(u_disc[Np] - u_disc[Np-1])
    u_jump = 0.5*(u_jump0+u_jump1) 
    print('u_cont_max: ',max(u_cont))
    print('u_disc_max: ',max(u_disc))
    print('u_jump:     ',u_jump)
    print('est_amp:    ',max(u_disc)-u_jump)
    print('u_cont_exact_max: ',max(u_cont_exact))
    print('u_disc_exact_max: ',max(u_disc_exact))
    print('error: ',abs(1-max(u_disc)))
    print('error_cont: ',abs(max(u_cont_exact)-max(u_cont)))
    print('error_disc: ',abs(max(u_disc_exact)-max(u_disc)))
    #=========================== PLOTTING Solution(1) ============================#
    fig = pyplot.figure();

    ylim_0 = list();
    ylim_1 = list();
    ylim_0.append(min(u_cont_exact));
    ylim_1.append(max(u_cont_exact));
    ylim_0.append(min(u_disc));
    ylim_1.append(max(u_disc));
    
    #plotting continuous data:
    pyplot.plot(x_cont_exact, u_cont_exact, '--k', label='Exact solution');
    ll = str("DGp")+ str(DG) + r'-$\beta$' + str(Beta) \
    +"_RK" + str(RK) +", CFL="+str(CFL)+", t=" + str(tt_);
    pyplot.plot(x_cont, u_cont, '-.b', label=ll);
    
    #plotting discontinuous numerical data:
    nn = size(x_disc);
    Np = N_disc_ppt;
    for i in range(0, nn - 1, Np):
        xx = x_disc[i:i + Np];
        uu = u_disc[i:i + Np];
        pyplot.plot(xx, uu, '-r');
    xx = x_disc[i:i + Np];
    uu = u_disc[i:i + Np];
    pyplot.plot(xx, uu, '-r', label='Discontinuous solution');

    pyplot.legend();
    pyplot.title("Plotting of Continuous solutions");
    pyplot.xlabel('X', labelpad=10);
    pyplot.ylabel('u(x)', labelpad=10);

    pyplot.xlim(min(x_cont_exact), max(x_cont_exact));
    pyplot.ylim(min(ylim_0) * 1.05, max(ylim_1) * 1.05);

    n_divisions = 8;
    xtick_dx = (x_cont_exact[-1] - x_cont_exact[0] ) / n_divisions;
    xlabels = arange(x_cont_exact[0], x_cont_exact[-1]+xtick_dx,xtick_dx);
    #xlabels = [0,10,20,30,40,50,60,70,80];
    xlocs = xlabels;
    pyplot.xticks(xlocs, xlabels);
    pyplot.grid()
    pyplot.xlim(min(x_cont_exact), max(x_cont_exact))
    pyplot.ylim(min(ylim_0), max(ylim_1)*1.3)

    fig.tight_layout()
    figname = dir1 + 'tempfig/' + 'p'+str(DG)+'RK'+RK+\
    '_beta'+str(Beta)+'_N'+Nelem+mm_name+str('_')+\
    str(tt_)+str('t_cont.png')
    fig.set_size_inches(15.0, 9.0, forward=True)
    pyplot.savefig(figname)
    
    #=========================== PLOTTING Solution(2) ============================#
    fig = pyplot.figure();
    
    ylim_0 = list();
    ylim_1 = list();
    ylim_0.append(min(u_disc_exact));
    ylim_1.append(max(u_disc_exact));
    ylim_0.append(min(u_disc));
    ylim_1.append(max(u_disc));

    #plotting discontinuous numerical data:
    label_disc = str("DGp")+ str(DG) + r'-$\beta$' + str(Beta) \
    +"_RK" + str(RK) +", CFL="+str(CFL)+", t=" + str(tt_); # discontinuous label
    for i in range(0, nn - 1, Np):
        xx = x_disc[i:i + Np];
        uu = u_disc[i:i + Np];
        pyplot.plot(xx, uu, '--r');
    xx = x_disc[i:i + Np];
    uu = u_disc[i:i + Np];
    pyplot.plot(xx, uu, '--r', label=label_disc);
    
    #plotting discontinuous exact data:
    nn = size(x_disc_exact)
    Np = N_disc_ppt
    label_disc = "Exact discontinuous" ; # discontinuous exact label
    for i in range(0, nn - 1, Np):
        xx = x_disc_exact[i:i + Np];
        uu = u_disc_exact[i:i + Np];
        pyplot.plot(xx, uu, '-k');
    xx = x_disc_exact[i:i + Np];
    uu = u_disc_exact[i:i + Np];
    pyplot.plot(xx, uu, '-k', label=label_disc);

    pyplot.legend();

    pyplot.title("Plotting of discontinuous solutions");
    pyplot.xlabel('X', labelpad=10);
    pyplot.ylabel('u(x)', labelpad=10);

    pyplot.xlim(min(x_cont_exact), max(x_cont_exact));
    pyplot.ylim(min(ylim_0) * 1.05, max(ylim_1) * 1.05);

    n_divisions = 8;
    xtick_dx = (x_cont_exact[-1] - x_cont_exact[0] ) / n_divisions;
    xlabels = arange(x_cont_exact[0], x_cont_exact[-1]+xtick_dx,xtick_dx);
    #xlabels = [0,10,20,30,40,50,60,70,80];
    xlocs = xlabels;
    pyplot.xticks(xlocs, xlabels);
    pyplot.grid()
    pyplot.xlim(min(x_cont_exact), max(x_cont_exact))
    pyplot.ylim(min(ylim_0), max(ylim_1)*1.3)

    fig.tight_layout()
    figname = dir1 + 'tempfig/' + 'p'+str(DG)+'RK'+RK+\
    '_beta'+str(Beta)+'_N'+Nelem+mm_name+str('_')+\
    str(tt_)+str('t_disc.png')
    fig.set_size_inches(15.0, 9.0, forward=True)
    pyplot.savefig(figname)
    
    #================== Plot FFT ==========================#
    label_dg = 'p'+ str(DG) + r'-$\eta$' + str(Beta) \
    +"_RK" + str(RK) # discontinuous label
    k_max_ = int(k_freq[-1]);
    print('k_max: ',k_max_, '  k_max_ex: ',int(k_freq_exact[-1]))
    print('exact Nelem:',size(x_cont_exact))
    fig, ax = plt.subplots()
    
    plt.plot(2*pi*k_freq_exact/(3*int(Nelem)), u_amp_exact, '-k',markevery=1\
        , label=r'Exact, E$_{tot}$= '+str(np.round(E_ex,4)))
    plt.plot(2*pi*k_freq/(3*int(Nelem)), u_amp, ':om',markevery=1, \
        label=label_dg+r', E$_{tot}$='+str(np.round(E_tot,4)))
    
    xlabels = ['0',r'$\pi$/6',r'$\pi$/4',r'$\pi$/3',r'$\pi$/2',r'3$\pi$/4',r'$\pi$'];
    xlocs = [0,pi/6,pi/4,pi/3,pi/2,3*pi/4,pi];
    plt.xticks(xlocs, xlabels);
    plt.xlim(0,pi)
    plt.legend();
    plt.xlabel('K/(P+1)', labelpad=2);
    plt.ylabel('|u|', labelpad=2);
    plt.grid()
    plt.title(str('t=') + str(tt_))

    fig.set_size_inches(13.0, 10.0, forward=True)
    fig.tight_layout(pad=0.0, w_pad=10.0, h_pad=10.0)
    
    figname = dir1 + 'tempfig/'+ 'fft_'+'p'+str(DG)+'RK'+RK+\
    '_eta'+str(Beta)+'_N'+Nelem+mm_name+str('_')+str(tt_)+str('t_comp_pi_')+str(K_den)+'.png'
    fig.set_size_inches(15.0, 9.0, forward=True)
    plt.savefig(figname,format='png')
    figname = dir1 + 'tempfig/'+ 'fft_'+'p'+str(DG)+'RK'+RK+\
    '_eta'+str(Beta)+'_N'+Nelem+mm_name+str('_')+str(tt_)+str('t_comp_pi_')+str(K_den)+'.eps'
    plt.savefig(figname,format='eps')
    
    #=========================== PLOTTING Errors Evolution ============================#
    # Read error data:
    res_dir = './Results/'
    fname = dir1 + 'errors/errors_N'+str(Nelem)+'_CFL'+str(CFL)+'_Beta1.00_'+str(T_period)+'T.dat'
    data = loadtxt(fname);  # continuous exact nodal solution
    time  = data[:, 0];
    L1err = data[:, 1];
    L2err = data[:, 2];

    #plotting
    fig = pyplot.figure();
    pyplot.plot(time, L2err,'-ok',time,L1err,'-.b')
    pyplot.xlabel('time')
    pyplot.ylabel('L2err')
    #pyplot.xscale('log')
    #pyplot.yscale('log')
    
    figname = dir1 + 'tempfig/' + 'errors_p'+str(DG)+'RK'+RK+\
    '_beta'+str(Beta)+'_N'+Nelem+mm_name+str('_')+\
    str(T_period)+str('T.png')
    fig.set_size_inches(15.0, 9.0, forward=True)
    pyplot.savefig(figname)
    
    #------------------------
    pyplot.show()

    return 'true'

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#                                  A D V E C T I O N--D I F F U S I O N
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
def plot_AdvecDiffus(diffus_scheme, mode, DG, RK, CFL, Nelem, T, dt_\
                , Beta, Epsilon, gamma, dir1, aver, nodal_exact, nodal_comp, discont ):

    Beta = Decimal(Beta.quantize(Decimal('.01')));
    Epsilon = Decimal(Epsilon.quantize(Decimal('.01')));

    #CFL = Decimal(CFL.quantize(Decimal('.001')));
    T = Decimal(T.quantize(Decimal('.001')));
    dt = float(dt_);

    fname = dir1 + nodal_exact + str("_") + str(T) + str("T.dat");
    data = loadtxt(fname);  # continuous exact nodal solution
    xn_exact = data[:, 0];
    un_exact = data[:, 1];

    del fname, data

    if mode == 'CFL_const' or mode == 'normal' or mode == 'test':
        fname = dir1 + discont + str("_N") + Nelem \
                + str("_CFL") + CFL + str("_Beta") + str(Beta) \
                + str("_Eps") + str(Epsilon) \
                + str("_") + str(T) + str("T.dat")
    elif mode == 'dt_const':
        fname = dir1 + discont + str("_N") + Nelem \
                + str("_dt") + dt_ + str("_Eps") + str(Epsilon) \
                + str("_") + str(T) + str("T.dat")

    data = loadtxt(fname)

    x_disc = data[:, 0]
    u_disc = data[:, 1]

    del fname, data

    nn = size(x_disc);

    if int(DG) == 0:
        Np = 2;
    elif int(DG) == 1:
        Np = 2;
    else:
        Np = 10;

    pyplot.figure();

    ylim_0 = list();
    ylim_1 = list();
    ylim_0.append(min(un_exact));
    ylim_1.append(max(un_exact));

    pyplot.plot(xn_exact, un_exact, '-k', label='Exact sol');

    for i in range(0, size(x_disc) - 1, Np):
        xx = x_disc[i:i + Np];
        uu = u_disc[i:i + Np];
        pyplot.plot(xx, uu, '-r');

    ylim_0.append(min(u_disc));
    ylim_1.append(max(u_disc));

    xx = x_disc[i:i + Np];
    uu = u_disc[i:i + Np];
    ll = str("Numerical sol, (") + r'$\varepsilon= $' + str(Epsilon) + " )";
    pyplot.plot(xx, uu, '-r', label=ll);

    pyplot.legend();

    title_a = str("DGp") + DG  + "-RK" + RK + "-" + diffus_scheme \
              + ", and penalty_param (" + r'$\varepsilon= $' \
              + str(Epsilon) + ")\n for CFL=" + str(CFL) \
              + " and at t/T=" + str(T);

    pyplot.title(title_a);
    pyplot.xlabel('X', labelpad=10);
    pyplot.ylabel('u(x)', labelpad=10);

    pyplot.xlim(min(xn_exact), max(xn_exact));
    pyplot.ylim(min(ylim_0) * 1.05, max(ylim_1) * 1.05);

    xlabels = linspace(min(xn_exact), max(xn_exact), 5);
    xlocs = xlabels;
    pyplot.xticks(xlocs, xlabels);
    pyplot.grid();

    pyplot.show()

    return 'true'
    
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#                   B U R G E R S--D E C A Y I N G--T U R B U L E N C E
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
def plot_burgers_decay_turb(dir_input,mode,DG,RK,CFL,Nelem,N_disc_ppt,tt_,dt_\
                         ,Beta,Epsilon,gamma_,diffus_scheme,cont_num,disc_num):
    
    nelem = int(Nelem)
    Np = int(int(DG)+1)
    nDOF = nelem *Np
    
    if  (mode == 'test') | (mode == 'dt_const'):
        mm_name = str('_dt') + dt_
        m_name = str('dt=') + dt_
    elif(mode == 'CFL_const') | (mode == 'normal'):
        mm_name = str('_CFL')+ str(CFL)     
        m_name = str('CFL=')+ str(CFL)             

    dt = float(dt_);

    diffus_scheme_name = diffus_scheme+'eta'+str(Epsilon)
    case_data_string = 'p'+DG+'RK'+RK+'_Ne'+str(Nelem)+'_Beta'+str(Beta)\
                                     + mm_name+ '_t'+ str(tt_)+'_'\
                                     +diffus_scheme_name
    print('case_name: ',case_data_string)
    print('dir_input: ',dir_input)
    # Read continuous data:
    fname = dir_input + cont_num + str("_N") + Nelem \
              + mm_name + str("_Beta") + str(Beta)\
              + str("_Eps") + str(Epsilon) \
              + str('_')+ str(tt_) + str("t.dat")
    data = loadtxt(fname)
    x_cont = data[:, 0]
    u_cont = data[:, 1]
    # compute fft:
    k_freq, u_amp, KEnerg = compute_fft(u_cont);
    del fname, data
    # Read discontinuous data:
    fname = dir_input+disc_num+'_N'+Nelem+mm_name+'_Beta'+str(Beta)+'_Eps'\
               +str(Epsilon)+'_'+str(tt_)+'t.dat'
    data = loadtxt(fname)
    x_disc = data[:, 0]
    u_disc = data[:, 1]
    N_disc_tot = size(x_disc);
    Ndisc_per_elem = N_disc_ppt

    ################ Plotting the solution #####################
    fig, ax = pyplot.subplots(frameon='True')
    ylim_0 = list();
    ylim_1 = list();
    ylim_0.append(min(u_cont));
    ylim_1.append(max(u_cont));
    pyplot.plot(x_cont, u_cont, '-k', label='continuous sol');
    for i in range(0, N_disc_tot - 1, Ndisc_per_elem):
        xx = x_disc[i:i + Np];
        uu = u_disc[i:i + Np];
        #pyplot.plot(xx, uu, '-r');

    ylim_0.append(min(u_disc));
    ylim_1.append(max(u_disc));

    xx = x_disc[i:i + Ndisc_per_elem];
    uu = u_disc[i:i + Ndisc_per_elem];
    ll = str("discontinuous sol");
    #pyplot.plot(xx, uu, '-r', label=ll);

    #pyplot.legend();

    title_a = diffus_scheme_name+("_p") + DG  + "-RK" + RK +" with " + m_name \
              + " and at t=" + str(tt_);
    pyplot.title(title_a);
    pyplot.xlabel('X', labelpad=10);
    pyplot.ylabel('u(x)', labelpad=10);
    #pyplot.xlim(min(x_cont), max(x_cont));
    #pyplot.ylim(min(ylim_0) * 1.05, max(ylim_1) * 1.05);
    xlabels = linspace(min(x_cont), max(x_cont), 5);
    xlocs = xlabels;
    #pyplot.xticks(xlocs, xlabels);
    pyplot.grid();

    fig.set_size_inches(13.0, 9.0, forward=True)
    fig.tight_layout(pad=0, w_pad=10.0, h_pad=10.0,rect=(0.0,0.0,1.0,0.985))
    temp_name = 'sol_vs_x_' +case_data_string
    figname = dir_input + str('/tempfig/eps/') + temp_name +'.eps'
    fig.savefig(figname,format='eps',bbox='tight')
    figname = dir_input + str('/tempfig/png/') + temp_name +'.png'
    plt.savefig(figname,format='png',bbox='tight')
    
    ###################### Plot the fft ##########################
    k_max_ = int(nDOF/2)
    fig1, ax1 = plot_KEnerg(k_freq[0:k_max_+1], KEnerg[0:k_max_+1]) # calling the plot function
    ax1.set_title(r'nDOFs$='+str(nDOF)+'$')
    temp_name = 'KE_'+case_data_string
    figname = dir_input + str('/tempfig/eps/')+temp_name+'.eps'
    fig1.savefig(figname,format='eps',bbox='tight')
    figname = dir_input + str('/tempfig/png/')+temp_name+'.png'
    fig1.savefig(figname,format='png',bbox='tight',dpi=150)
    
    pyplot.show()
    
    return 'true'
                     
def plot_dissipation_rates(dir_input,mode,DG,RK,CFL,Nelem,N_disc_ppt,tt_,dt_\
                    ,Beta,Epsilon,gamma_,diffus_scheme,cont_num,disc_num):

    nelem = int(Nelem)
    Np = int(int(DG)+1)
    nDOF = nelem *Np

    if  (mode == 'test') | (mode == 'dt_const'):
        mm_name = str('_dt') + dt_
        m_name = str('dt=') + dt_
    elif(mode == 'CFL_const') | (mode == 'normal'):
        mm_name = str('_CFL')+ str(CFL)
        m_name = str('CFL=')+ str(CFL)

    dt = float(dt_)
    diffus_scheme_name = diffus_scheme+'eta'+str(Epsilon)
    case_string_name = 'p'+ DG + 'RK' + RK +'_Ne'+ str(Nelem)+ mm_name\
                          +'_Beta'+str(Beta)+'_'+diffus_scheme_name

    # Read KE vs time data:
    fname = dir_input + 'time_data/wave_energy' + str("_N") + Nelem \
             + mm_name + str("_Beta") + str(Beta)\
             + str("_Eps") + str(Epsilon) + '.dat'
    data = loadtxt(fname)
    time_array = data[:, 0]
    KE_array   = data[:, 1]
    G_array    = data[:, 2]
    # Computing the dissipation rates:
    dKE_dt = compute_dEkdt(time_array,KE_array)
    Dt_dissip = zeros(size(time_array))
    for i in range(0,size(time_array)):
        temp_tt = Decimal(time_array[i])
        temp_tt = Decimal(temp_tt.quantize(Decimal('.001')))
        fname = dir_input + cont_num + str("_N") + Nelem \
                  + mm_name + str("_Beta") + str(Beta)\
                  + str("_Eps") + str(Epsilon) \
                  + str('_')+ str(temp_tt) + str("t.dat")
        data = loadtxt(fname)
        u_temp_ = data[:, 1]
        k_freq_temp, u_amp_temp, KEnerg_temp = compute_fft(u_temp_)
        Dt_dissip[i] = compute_Dt_dissipation_rate(k_freq_temp, KEnerg_temp\
                                  , float(gamma_))
    ################ Plotting KE, dKE_dt, & Dt vs time #####################
    fig, ax = plot_E_vs_time(time_array,KE_array)
    temp_name = 'Et_' + case_string_name
    figname = dir_input + str('/tempfig/eps/') + temp_name+'.eps'
    fig.savefig(figname,format='eps',bbox='tight')
    figname = dir_input + str('/tempfig/png/') + temp_name+'.png'
    fig.savefig(figname,format='png',bbox='tight',dpi=150)

    fig, ax = plot_dissip_rates(time_array,dKE_dt,Dt_dissip)
    temp_name = 'dEt_dt_' + case_string_name
    figname = dir_input + str('/tempfig/eps/') + temp_name+'.eps'
    fig.savefig(figname,format='eps',bbox='tight')
    figname = dir_input + str('/tempfig/png/') + temp_name+'.png'
    fig.savefig(figname,format='png',bbox='tight',dpi=150)

    plt.show()
