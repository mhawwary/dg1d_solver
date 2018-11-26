from matplotlib import pyplot as plt   #and the useful plotting library
from matplotlib import pyplot 
from matplotlib import ticker

import numpy as np
from numpy import sin,cos,pi,linspace,ones,zeros,abs,min,max,exp,sqrt\
    , shape, empty_like , size, loadtxt, savetxt, arange , log , hstack, vstack, trapz, append, savez
    
import argparse
from decimal import Decimal
import csv

from subprocess import call, run, PIPE, Popen
from sys_cmd_toolbox import system_process    # locally defined file

from numpy import fft

pyplot.rc('legend',**{'loc':'best'});
pyplot.rcParams[u'legend.fontsize'] = 15
pyplot.rcParams[u'legend.edgecolor']='white'
pyplot.rcParams[u'font.weight']='normal'
#pyplot.rcParams['font.serif']='false'
pyplot.rcParams[u'xtick.labelsize']=14
pyplot.rcParams[u'ytick.labelsize']=14
pyplot.rcParams[u'axes.titlesize']=18
pyplot.rcParams[u'axes.labelsize']=15
pyplot.rcParams[u'axes.spines.right']='false';
pyplot.rcParams[u'axes.spines.top']='false';
pyplot.rcParams[u'lines.linewidth'] = 1.5;
pyplot.rcParams[u'lines.markersize'] = 8;


#==============================================================================================================#
#                                        General FFT Functions
#==============================================================================================================#
def python_inputfile_options(fname,option_name,option_value):
    return;

def load_data(fname):

    data = loadtxt(fname);  # continuous exact nodal solution
    x_data_ = data[:, 0]
    u_data_ = data[:, 1]

    return x_data_, u_data_

def compute_fft(u_data_):

    sample_size = size(u_data_);
    if sample_size%2==0:
        normalize_factor = (sample_size/2.)+1.
    elif sample_size%2==1:
        normalize_factor = (sample_size+1.)/2.   
    
    u_amp = abs(fft.rfft(u_data_))/normalize_factor
    KEnerg = 0.5 * u_amp**2

    sample_rate = 1/sample_size;
    
    #print('normalize_factor,',normalize_factor, '  sample_size,',sample_size)

    k_freq = fft.rfftfreq(sample_size, sample_rate);
    
    #print('kmax, ',k_freq[-1])

    return k_freq, u_amp, KEnerg

def compute_ensemble_fft(x_,u_data_):

    sample_size = size(u_data_[1,:]);
    aver_size = size(u_data_[:,1]);

    if sample_size%2==0:
        normalize_factor = (sample_size/2.)+1.
    elif sample_size%2==1:
        normalize_factor = (sample_size+1.)/2.
    
    sample_rate = 1/sample_size;
    k_freq = fft.rfftfreq(sample_size, sample_rate);

    #u_x = fft.irfft(fft.rfft(u_data_[0,:]),size(x_))
    #print('x: ',size(x_),', u_x: ', size(u_x))
    #pyplot.figure()
    #pyplot.plot(x_,u_x,'--r',x_,u_data_[0,:],'-k')
    #pyplot.show()
    
    u_amp = abs(fft.rfft(u_data_[0,:]))
    for i in range(1,aver_size):
        u_amp += abs(fft.rfft(u_data_[i,:]))

    u_amp = u_amp / (aver_size*normalize_factor)

    KEnerg = 0.5 * u_amp**2.0

    return k_freq, u_amp, KEnerg

def plot_KEnerg(kfreq, KEnerg):

    # Compute Initial and Exact Spectrum:
    kfreq_ex, KEnerg_ex, kfreq_init, KEnerg_init = compute_exact_init_Energy()
    E_tot = compute_Etotal(kfreq,KEnerg)
    Einit_tot = compute_Etotal(kfreq_init,KEnerg_init)

    fig, ax = plt.subplots(frameon='True')
    ax.plot(kfreq_ex, KEnerg_ex*1.4, ':k',lw=0.7)  # theoretical slope of k^-2
    ax.plot(kfreq_init, KEnerg_init, '--r', label=r'KE$_{init}$, E$_{tot}='+str(np.round(Einit_tot,5))+'$',lw=0.7)
    #ax.plot(k_freq_DNS[1:], KE_DNS[1:],'-r', label=DNS_label,lw=0.5)
    markers_on = [35,50,60,70,90,110,125,150,170,200,250,300,400,500]
    ax.plot(kfreq[1:], KEnerg[1:], '-b', label=r'KE$_{num}$, E$_{tot}='+str(np.round(E_tot,5))+'$'\
                   , lw=0.85, markevery=markers_on )

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(fontsize=19,edgecolor='black',loc='upper right')
    ax.set_facecolor('white')
    ax.set_xlabel(r'$k$', labelpad=2,fontsize=24);
    ax.set_ylabel(r'$E$', labelpad=2, fontsize=24);
    ax.set_xlim(10**0, 10**4)
    ax.set_ylim(10**-10, 10 ** -1)
    ax.tick_params(axis='both', which='both', labelsize=24)
    #ax.set_title(r'nDOFs$='+str(nDOF)+'$')
    ax.text(475, 1.5*10**-4, r'k$^{-2}$', fontsize=22)
    #ax.xaxis.label.set_size(20)  
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    ax.spines["bottom"].set_visible(True)

    fig.set_size_inches(13.0, 9.0, forward=True)
    fig.tight_layout(pad=0, w_pad=10.0, h_pad=10.0,rect=(0.0,0.0,1.0,0.965))
    
    return fig, ax

def plot_E_vs_time(time,KE):
    fig, ax = plt.subplots(frameon='True')
    ax.plot(time,KE,'-r')
    ax.set_xlabel('time', labelpad=2,fontsize=18)
    ax.set_ylabel('E', labelpad=2, fontsize=18)
    ax.grid()
    ax.set_xlim(0.0, np.max(time))
    ax.tick_params(axis='both', which='both', labelsize=16)
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    ax.spines["bottom"].set_visible(True)
    fig.set_size_inches(8.0, 6.0, forward=True)
    fig.tight_layout(pad=0, w_pad=10.0, h_pad=10.0\
                  ,rect=(0.0,0.0,1.0,0.985))
    return fig, ax

def plot_dissip_rates(time,dE_dt,Dt):
    fig, ax = plt.subplots(frameon='True')
    ax.plot(time,-dE_dt,'-^r',label='-dE/dt',lw=0.7)
    ax.plot(time,Dt,'--ob',label='D(t)',lw=0.7)
    ax.set_xlabel('time', labelpad=4,fontsize=16)
    ax.set_ylabel('Kinetic energy dissipation rate'\
          ,labelpad=4, fontsize=16)
    ax.set_xlim(0.0, np.max(time))
    #ax.set_ylim(0.4,1.80)
    ax.grid()
    ax.legend(fontsize=16,edgecolor='black',loc='upper right')
    ax.tick_params(axis='both', which='both', labelsize=14)
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    ax.spines["bottom"].set_visible(True)
    fig.set_size_inches(8.0, 6.0, forward=True)
    fig.tight_layout(pad=0, w_pad=10.0, h_pad=10.0\
                ,rect=(0.0,0.0,1.0,0.985))
    return fig, ax
    
def plot_fft(k_freq,u_amp):

    # Plotting the FFT:
    fig, ax = pyplot.subplots(figsize=(10, 8));
    pyplot.plot(k_freq[0:], u_amp[0:], '-b', label='u_amp')
    #pyplot.plot(k_init, E_init, '-r', label='Initial KE')
    #pyplot.plot(k_ex, E_ex, '--k', label=r'E = const * k$^{-2}$')

    #pyplot.xscale('log')
    pyplot.yscale('log')

    pyplot.grid()
    pyplot.legend()

    pyplot.xlabel(r'k $\Delta x$', labelpad=10);
    pyplot.ylabel('|u|', labelpad=10);

    n_divisions = 8
    xtick_dx = (4*pi - 0.0 ) / n_divisions
    xlocs = arange(0.0, 4*pi +xtick_dx,xtick_dx)
    xlabels = ['0.0',r'$\pi$/2',r'$\pi$',r'3$\pi$/2',r'2$\pi$',r'5$\pi$/2',r'3$\pi$',r'7$\pi$/2',r'4$\pi$']
    pyplot.xticks(xlocs, xlabels)
    pyplot.xlim(0,4*pi)
    fig.tight_layout()

    pyplot.show()

    return "true"
    
def compute_exact_init_Energy():
    # preparing the exact spectrum for comparison:
    ko = 10.0
    n_pts = 2000
    kmax = 2000
    k_init = linspace(1, kmax, n_pts)
    A = 2 * ko ** -5 / (3 * sqrt(pi))
    E_init = A * k_init ** 4 * exp(-(k_init / ko) ** 2)
    k_ex = linspace(30, 10**4)
    E_ex = 20 * k_ex ** -2

    return k_ex, E_ex, k_init, E_init

def compute_Etotal(kfreq,KEnerg):

    E = 2.0 * trapz(KEnerg,kfreq)
    
    return E
    
def compute_Dt_dissipation_rate(kfreq,KEnerg, thermal_diffus):

    nu = thermal_diffus
    II = 2.0 * KEnerg*(kfreq**2) 
    Dt = 2.0 * nu  * trapz(II,kfreq)
    
    return Dt

def compute_dEkdt(time,Ek):
    nt = size(time);
    dEk_dt = zeros([nt])

    dEk_dt[0] = (Ek[1]-Ek[0])/(time[1]-time[0])
    for i in range(1,nt):
        dEk_dt[i] = (Ek[i]-Ek[i-1])/(time[i]-time[i-1])

    return dEk_dt
    
def savez_spectrum(k_freq, u_amp, KEnerg, fname):
    savez(fname, k_freq=k_freq, u_amp=u_amp, KEnerg=KEnerg)
    return

def loadz_spectrum(fname):
    data = np.load(fname)
    k_freq = data['k_freq']
    u_amp = data['u_amp']
    KE_fd = data['KEnerg']
    return k_freq, u_amp, KE_fd
    
def savez_spectrum_linearWave(k_freq, u_amp, KEnerg, G, fname):
    # Saving the spectrum along with G for a linear wave, i.e.
    # , a sine or Gaussian
    savez(fname, k_freq=k_freq, u_amp=u_amp, KEnerg=KEnerg, G=G)
    return

def loadz_spectrum_linearWave(fname):
    # Loading the spectrum along with G for a linear wave, i.e.
    # , a sine or Gaussian
    data = np.load(fname+str('.npz'))
    k_freq = data['k_freq']
    u_amp = data['u_amp']
    KEnerg = data['KEnerg']
    G = data['G']
    return k_freq, u_amp, KEnerg, G

#==============================================================================================================#
#                                      Discontinuous Galerkin Functions
#==============================================================================================================#
def DG_fft_input_reader(input_file_):      
    with open(input_file_) as file: 
            reader=csv.reader(file, delimiter=':');

            for row in reader:
                if row[0] =='dir_result':    
                    dir_result = str(row[1]);
                elif row[0] =='dir_input':    
                    dir_input = str(row[1]);
                elif row[0] =='errors':   
                    errors_dir = str(row[1]);
                elif row[0] =='aver':   
                    aver_dir = str(row[1]);
                elif row[0] == 'cont_exact':
                    exact_sol_dir = str(row[1]);
                elif row[0] == 'cont_num':
                    nodal_comp = str(row[1]);
                elif row[0] == 'cont_unsteady_num':
                    dg_cont_sol_dir = str(row[1]);
                elif row[0] == 'disc_unsteady_num':
                    dg_disc_sol_dir = str(row[1]);
                elif row[0] == 'discont': 
                    discont = str(row[1]);
                elif row[0]=='Eqn_set':
                    eqn_set=str(row[1]);
                elif row[0]=='Eqn_type':
                    eqn_type=str(row[1]);
                elif row[0]=='wave_form':
                    wave_form=str(row[1]);
                elif row[0]=='mode':
                    mode=str(row[1]);
                elif row[0]=='Diffusion_scheme':
                    diffus_scheme=str(row[1]);
                elif row[0] == 'thermal_diffus':
                    thermal_diffus = float(row[1])
                elif row[0] == 'Epsilon':
                    Epsilon=Decimal(row[1])
                elif row[0] == 'Beta':
                    Beta=Decimal(row[1])
                elif row[0]=='p':
                    poly_order=str(row[1]);
                elif row[0] == 'RK':    
                    RK=str(row[1]);
                elif row[0] == 'Nelem':    
                    Nelem=str(int(row[1]));
                elif row[0]=='CFL':
                    CFL=Decimal(row[1])
                elif row[0] == 'dt':
                    dt_=str(row[1])
                elif row[0] == 't_end':
                    t_end_=str(row[1]);
                elif row[0] == 'T':     
                    T=Decimal(row[1])
                elif row[0] == 'N_ensembles_fft':
                    N_aver_fft = str(row[1])
                elif row[0] == 'time_for_fft_compute':
                    tt_ = Decimal(row[1])
                
    CFL = Decimal(CFL.quantize(Decimal('0.0001')))
    tt_ = Decimal(tt_.quantize(Decimal('0.001')))
    
    cmd = ['mkdir',dir_result]
    cmd_out,cmd_err=system_process(cmd,5000)
    cmd=['mkdir',dir_result+'data/']
    cmd_out,cmd_err=system_process(cmd,5000)
    cmd=['mkdir',dir_result+'fig/']
    cmd_out,cmd_err=system_process(cmd,5000)
    
    return dir_result, dir_input, mode, poly_order, RK, Beta, Epsilon, Nelem\
    , CFL, dt_, T, tt_, N_aver_fft, dg_cont_sol_dir, dg_disc_sol_dir

def load_ensemble_data_DG(dir_input, mode, poly_order, RK, cont_sol\
        , Nelem, Beta, Epsilon, CFL, dt_, tt_, N_aver_fft_):

    N_aver_fft = int(N_aver_fft_)
   
    if (mode=='test') | (mode =='dt_const'):
        mm_name = '/' + cont_sol + str("_N") + str(Nelem) \
                             + str("_dt") + dt_ + str("_Beta") + str(Beta) \
                             + str("_Eps") + str(Epsilon) \
                             + str("_") + str(tt_) + str("t.dat")
    else:
        mm_name = '/' + cont_sol + str("_N") + str(Nelem) \
                             + str("_CFL") + str(CFL) + str("_Beta") + str(Beta) \
                             + str("_Eps") + str(Epsilon) \
                             + str("_") + str(tt_) + str("t.dat") 
                             
    case_no = str('01')
    mm_dir = dir_input + str('case') + case_no
    fft_inputfile_name = mm_dir + mm_name
    
    x_, u_ = load_data(fft_inputfile_name)

    for i in range(2, N_aver_fft + 1):
        if (i < 10):
            ii = str('0') + str(i)
        else:
            ii = str(i)
            
        case_no = str(ii)
        mm_dir = dir_input + str('case') + case_no
        fft_inputfile_name = mm_dir + mm_name
        x_, u_temp = load_data(fft_inputfile_name)
        u_ = vstack((u_, u_temp))
        
        for jj in range(0,size(u_temp)):
            if (np.isnan(u_temp[jj])):
                print('u_nan at case:',ii, ' jj:',jj,'  x:',x_[jj])
                break;
                

    return x_,u_

def compute_cont_fft_mesh_refine_DG(dir_input, dir_result, mode, poly_order, RK\
, cont_sol, Nelem, Beta, Epsilon, CFL, dt_, tt_, N_aver_fft):

    if (mode=='test') | (mode =='dt_const'):
        sim_name = '_dt'+ dt_
    else:
        sim_name = '_CFL' + str(CFL)
        
    # Compute Initial and Exact Spectrum:
    k_ex, E_ex, k_init, E_init = compute_exact_init_Energy()

    print('ensemble_size: ', N_aver_fft)
    print('Nelem: ',0,',    E_init: ', max(E_init), ',     Et: ',compute_Etotal(k_init,E_init))
    
    # Computing and Plotting FFT:
    fig, ax = pyplot.subplots(figsize=(10, 8))

    for i in range(0,size(Nelem)):
        x_,u_ = load_ensemble_data_DG(dir_input, mode, poly_order, RK, cont_sol\
        , Nelem[i], Beta, Epsilon, CFL, dt_, tt_, N_aver_fft)
        k_freq, u_amp, KEnerg  = compute_ensemble_fft(x_,u_)
        fname = dir_result + 'data/' + str('spectrum_N') + str(Nelem[i]) + sim_name + str('_t') + str(tt_) 
        savez_spectrum(k_freq, u_amp, KEnerg, fname)
        mesh_label = "Nelem " + str(Nelem[i])
        pyplot.plot(k_freq[1:], KEnerg[1:], label=mesh_label)
        print('Nelem: ',Nelem[i],',    E_max:  ',max(KEnerg[1:]),',     Et: ',compute_Etotal(k_freq[1:],KEnerg[1:])\
        , ',    Dt:',compute_Dt_dissipation_rate(k_freq[1:],KEnerg[1:], 2.0e-4))
        
        #Emax_fname = dir_result + 'data/temporal_data' + sim_name + str('_t') + str(tt_) + '.dat'
        #save_text_()
        
        #E_print = numpy.transpose([Nelem,Emax,Etot,Dt]);

#savetxt(Emax_fname, E_print\
#, fmt="%02d"+" %1.10f"+" %1.10f"+" %1.10f"\
#,header="Nelem, , order_L1_aver, order_L2_proj, order_L2_aver"\
#,delimiter=' ');
        
    print('KEnerg[0],',KEnerg[0])
        
    pyplot.plot(k_init, E_init, '-k', label='Initial KE')
    pyplot.plot(k_ex, E_ex, '--k', label=r'E = const * k$^{-2}$')
    pyplot.xscale('log')
    pyplot.yscale('log')
    pyplot.grid()
    pyplot.legend()
    pyplot.xlabel('k', labelpad=10);
    pyplot.ylabel('E(k)', labelpad=10);
    pyplot.xlim(10**0, 10**5)
    pyplot.ylim(10**-12, 10 ** 0)
    fig.tight_layout()
    figname = dir_result + 'fig/' + str('EnergSpect_') + str('p') + poly_order + str('_RK') + RK \
    + sim_name + str('_Ns') + str(N_aver_fft)+ str("_") \
    + str(tt_) + str("t.png")
    pyplot.savefig(figname)
    #pyplot.show()

    return

def perform_dt_convergence_study_DG(dir_input, dir1, DG, RK, cont_comp, Nelem, Beta, Epsilon, dt_, tt_, N_aver_fft):

    # Compute Initial Spectrum:
    k_ex, E_ex, k_init, E_init = compute_exact_init_Energy()
    
    print('ensemble_size: ', N_aver_fft)
    print('             E_init: ', max(E_init), ',     Et: ',compute_Etotal(k_init,E_init))
    
    # Computing and Plotting FFT:
    fig, ax = pyplot.subplots(figsize=(10, 8))

    for i in range(0,size(dt_)):
        x_,u_ = load_ensemble_data_DG(dir_input, DG, RK, cont_comp, Nelem, Beta, Epsilon, dt_[i], tt_, N_aver_fft)
        k_freq, u_amp, KEnerg  = compute_ensemble_fft(x_,u_)
        dt_label = str('dt ') + dt_[i]
        pyplot.plot(k_freq[1:], KEnerg[1:], label=dt_label)
        print('dt: ',dt_[i],',    E_max:  ',max(KEnerg),',     Et: ',compute_Etotal(k_freq,KEnerg)\
        ,',    Dt:',compute_Dt_dissipation_rate(k_freq,KEnerg, 2.0e-4))

    pyplot.plot(k_init, E_init, '-k', label='Initial KE')
    pyplot.plot(k_ex, E_ex, '--k', label=r'E = const * k$^{-2}$')
    pyplot.xscale('log')
    pyplot.yscale('log')
    pyplot.grid()
    pyplot.legend()
    pyplot.xlabel('k', labelpad=10);
    pyplot.ylabel('E(k)', labelpad=10);
    pyplot.ylim(10 ** -12, 10 ** 0)
    fig.tight_layout()
    
    figname = dir1 + str('dt_conv_') + str('p') + DG + str('_RK') + RK \
    + str('_Ne')+str(Nelem) + str('_Ns') + str(N_aver_fft)+ str("_") \
    + str(tt_) + str("t.png")
    pyplot.savefig(figname)
    pyplot.show()
    
def compare_fft_DG(dir_dg, DG1, dg1_sol, DG, dg_sol, Nelem_dg, Beta, Epsilon, RK, dt_, tt_, N_aver_fft):
    # Compute Initial Spectrum:
    k_ex, E_ex, k_init, E_init = compute_exact_init_Energy()
    # Compute DG and FD Spectrums:
    k_freq_dg1, u_amp_dg1, KE_dg1 = compute_fft_DG1(dir_dg, DG1, RK, dg1_sol, Nelem_dg, Beta, Epsilon, dt_, tt_, N_aver_fft)
    k_freq_dg, u_amp_dg, KE_dg = compute_fft_DG(dir_dg, DG, RK, dg_sol, Nelem_dg, Beta, Epsilon, dt_, tt_, N_aver_fft)
    
    #Plotting The Energy Spectrums:
    fig, ax = pyplot.subplots(figsize=(10, 8))
    dg_label = str('DGp') + DG + str('_RK') + RK + str(', Nelem: ') + Nelem_dg 
    dg1_label = str('DG1p') + DG1 + str('_RK') + RK + str(', Nelem: ') + Nelem_dg
    pyplot.plot(k_freq_dg[1:], KE_dg[1:], label=dg_label)
    pyplot.plot(k_freq_dg1[1:], KE_dg1[1:], label=dg1_label)
    
    pyplot.plot(k_init, E_init, '-k', label='Initial KE')
    pyplot.plot(k_ex, E_ex, '--k', label=r'E = const * k$^{-2}$')
    pyplot.xscale('log')
    pyplot.yscale('log')
    pyplot.grid()
    pyplot.legend()
    pyplot.xlabel('k', labelpad=10);
    pyplot.ylabel('E(k)', labelpad=10);
    pyplot.xlim(10**0, 10**5)
    pyplot.ylim(10**-15, 10 ** 0)
    fig.tight_layout()
    figname = str('./Results/') + str('EnergSpect_') + str('DGp') + DG + str('_vs_DG1') +  DG1 \
    + str('_RK') + RK +str('_dt') + dt_+ str('_Ns') + str(N_aver_fft) + str("_") \
    + str(tt_) + str("t.png")
    pyplot.savefig(figname)
    pyplot.show()

#==============================================================================================================#
#                                      Finite-Difference Functions
#==============================================================================================================#
def FD_fft_input_reader(input_file_):
    with open(input_file_) as file: 
        reader=csv.reader(file, delimiter=':');

        for row in reader:
            if row[0] =='dir_result':    
                dir_result = str(row[1])
            elif row[0] =='dir_input':    
                dir_input = str(row[1])
            elif row[0] =='errors':   
                errors_dir = str(row[1])
            elif row[0] == 'exact':
                exact_sol_dir = str(row[1])
            elif row[0] == 'numerical':
                num_sol_dir = str(row[1])
            elif row[0] =='unsteady_num':   
                fd_sol_dir = str(row[1])
            elif row[0]=='Eqn_set':
                eqn_set=str(row[1]);
            elif row[0]=='Eqn_type':
                eqn_type=str(row[1]);
            elif row[0]=='wave_form':
                wave_form=str(row[1]);
            elif row[0]=='mode':
                mode=str(row[1]);
            elif row[0]=='FDOA':
                FD=str(row[1])
                if int(FD)==2:
                    FDOA = FD + str('nd')
                elif int(FD)==3:
                    FDOA = FD + str('rd')
                else:
                    FDOA = FD + str('th')
            elif row[0] == 'RK':    
                RK=str(row[1]);
            elif row[0] == 'Nelem':    
                Nelem=str(int(row[1]));
            elif row[0] == 'Nexact':    
                Nexact=str(int(row[1]));
            elif row[0]=='CFL':
                CFL=Decimal(row[1])
            elif row[0] == 'dt':
                dt_=str(row[1]);
            elif row[0] == 't_end':
                t_end=str(row[1]);
            elif row[0] == 'T':     
                T=Decimal(row[1]);
            elif row[0] == 'N_ensembles_fft':
                N_aver_fft = str(row[1])
            elif row[0] == 'time_for_fft_compute':
                tt_ = Decimal(row[1])

    CFL = Decimal(CFL.quantize(Decimal('0.0001')))
    tt_ = Decimal(tt_.quantize(Decimal('0.001')))
    
    cmd = ['mkdir',dir_result]
    cmd_out,cmd_err=system_process(cmd,5000)
    cmd=['mkdir',dir_result+'data/']
    cmd_out,cmd_err=system_process(cmd,5000)
    cmd=['mkdir',dir_result+'fig/']
    cmd_out,cmd_err=system_process(cmd,5000)
                
    return dir_result, dir_input, mode, FDOA, RK, Nelem, CFL, dt_, tt_, N_aver_fft, fd_sol_dir

def compute_cont_fft_mesh_refine_FD(dir_input, dir_result, mode, FDOA\
            , RK, sol_dir, Nelem, CFL, dt, tt_, N_aver_fft):
    
    if (mode=='test') | (mode =='dt_const'):
        sim_name = '_dt'+ dt
    else:
        sim_name = '_CFL' + str(CFL)

    # Compute Initial Spectrum:
    k_ex, E_ex, k_init, E_init = compute_exact_init_Energy()

    print('ensemble_size: ', N_aver_fft)
    print('Nelem: ',0,',    E_init: ', max(E_init), ',     Et: ',compute_Etotal(k_init,E_init))
    
    # Computing and Plotting FFT:
    fig, ax = pyplot.subplots(figsize=(10, 8))

    for i in range(0,size(Nelem)):
        x_,u_ = load_ensemble_data_FD(dir_input, mode, FDOA, RK, sol_dir, Nelem[i], CFL, dt, tt_, N_aver_fft)
        k_freq, u_amp, KEnerg  = compute_ensemble_fft(x_,u_)
        fname = dir_result + 'data/' + str('spectrum_N') + str(Nelem[i]) + sim_name + str('_t') + str(tt_) 
        savez_spectrum(k_freq, u_amp, KEnerg, fname)
        mesh_label = "Nnodes " + str(Nelem[i]+1)
        pyplot.plot(k_freq[1:], KEnerg[1:], label=mesh_label)
        print('Nelem: ',Nelem[i],',    E_max:  ',max(KEnerg[1:])\
        ,',     Et: ',compute_Etotal(k_freq[1:],KEnerg[1:])\
        ,',    Dt:',compute_Dt_dissipation_rate(k_freq[1:],KEnerg[1:], 2.0e-4))
       
    pyplot.plot(k_init, E_init, '-k', label='Initial KE')
    pyplot.plot(k_ex, E_ex, '--k', label=r'E = const * k$^{-2}$')
    pyplot.xscale('log')
    pyplot.yscale('log')
    pyplot.grid()
    pyplot.legend()
    pyplot.xlabel('k', labelpad=10);
    pyplot.ylabel('E(k)', labelpad=10);
    pyplot.xlim(10**0, 10**5)
    pyplot.ylim(10**-12, 10 ** 0)
    fig.tight_layout()
    figname = dir_result + 'fig/' + str('KE_N') + str(Nelem[0]) + str('_FD') + FDOA + str('_RK') + RK \
    + sim_name + str('_Ns') + str(N_aver_fft)+ str("_") \
    + str(tt_) + str("t.png")
    pyplot.savefig(figname)
    #pyplot.show()

    return

def load_ensemble_data_FD(dir_input, mode, FDOA, RK, sol_dir, Nelem, CFL, dt_, tt_, N_aver_fft_):

    N_aver_fft = int(N_aver_fft_)
   
    if (mode=='test') | (mode =='dt_const'):
        mm_name = '/' + sol_dir + str("_N") + str(Nelem) \
        + str("_dt") + dt_ + str("_") + str(tt_) + str("t.dat")
    else:
        mm_name = '/' + sol_dir + str("_N") + str(Nelem) \
        + str("_CFL") + str(CFL) + str("_") + str(tt_) + str("t.dat")
        
    case_no = str('01')
    mm_dir = dir_input + str('case') + case_no
    fft_inputfile_name = mm_dir + mm_name
    x_, u_ = load_data(fft_inputfile_name)

    for i in range(2, N_aver_fft + 1):
        if (i < 10):
            ii = str('0') + str(i)
        else:
            ii = str(i)
        case_no = str(ii)
        mm_dir = dir_input + str('case') + case_no
        fft_inputfile_name = mm_dir + mm_name
        x_, u_temp = load_data(fft_inputfile_name)
        u_ = vstack((u_, u_temp))

    return x_,u_
      
def perform_dt_convergence_study_FD(dir2, FDOA, RK, sol_comp, Nelem, dt_, tt_, N_aver_fft):

    # Compute Initial Spectrum:
    k_ex, E_ex, k_init, E_init = compute_exact_init_Energy()
    
    # Computing and Plotting FFT:
    fig, ax = pyplot.subplots(figsize=(10, 8))

    for i in range(0,size(dt_)):
        x_,u_ = load_ensemble_data_FD(FDOA, RK, sol_comp, Nelem, dt_[i], tt_, N_aver_fft)
        k_freq, u_amp, KEnerg  = compute_ensemble_fft(x_,u_)
        dt_label = str('dt ') + dt_[i]
        pyplot.plot(k_freq[1:], KEnerg[1:], label=dt_label)
        print('dt: ',dt_[i],',    E_max:  ',max(KEnerg),',     Et: ',compute_Etotal(k_freq,KEnerg)\
        ,',    Dt:',compute_Dt_dissipation_rate(k_freq,KEnerg, 2.0e-4))

    pyplot.plot(k_init, E_init, '-k', label='Initial KE')
    pyplot.plot(k_ex, E_ex, '--k', label=r'E = const * k$^{-2}$')
    pyplot.xscale('log')
    pyplot.yscale('log')
    pyplot.grid()
    pyplot.legend()
    pyplot.xlabel('k', labelpad=10);
    pyplot.ylabel('E(k)', labelpad=10);
    pyplot.ylim(10 ** -12, 10 ** 0)
    fig.tight_layout()
    
    figname = dir2 + str('dt_conv_') + str('FD') + FDOA + str('_RK') + RK \
    + str('_Ne')+str(Nelem)+ str('_Ns') + str(N_aver_fft)+ str("_") \
    + str(tt_) + str("t.png")
    
    pyplot.savefig(figname)
    pyplot.show()
    
    
#==============================================================================================================#
#                                     Compare FD and DG 
#==============================================================================================================#
def compare_FD_DG_fft_input_reader(input_file_):      
    with open(input_file_) as file: 
            reader=csv.reader(file, delimiter=':');

            for row in reader:
                if row[0] =='dir_result_dg':    
                    dir_result_dg = str(row[1]);
                elif row[0] =='dir_result_fd':    
                    dir_result_fd = str(row[1]);
                elif row[0]=='mode':
                    mode=str(row[1]);
                elif row[0] == 'thermal_diffus':
                    thermal_diffus = float(row[1])
                elif row[0]=='p':
                    poly_order=str(row[1]);
                elif row[0] == 'RK':    
                    RK=str(row[1]);
                elif row[0] == 'Nelem_dg':    
                    Nelem_dg=str(int(row[1]));
                elif row[0]=='CFL_dg':
                    CFL_dg=Decimal(row[1])
                elif row[0]=='FDOA':
                    FD=str(row[1])
                    if int(FD)==2:
                        FDOA = FD + str('nd')
                    elif int(FD)==3:
                        FDOA = FD + str('rd')
                    else:
                        FDOA = FD + str('th')
                elif row[0] == 'Nelem_fd':    
                    Nelem_fd=str(int(row[1]));
                elif row[0]=='CFL_fd':
                    CFL_fd=Decimal(row[1])
                elif row[0] == 'dt':
                    dt_=str(row[1])
                elif row[0] == 't_end':
                    t_end_=str(row[1]);
                elif row[0] == 'T':     
                    T=Decimal(row[1])
                elif row[0] == 'N_ensembles_fft':
                    N_aver_fft = str(row[1])
                elif row[0] == 'time_for_fft_compute':
                    tt_ = Decimal(row[1])
                
    CFL_dg = Decimal(CFL_dg.quantize(Decimal('0.0001')))
    CFL_fd = Decimal(CFL_fd.quantize(Decimal('0.0001')))
    tt_ = Decimal(tt_.quantize(Decimal('0.001')))
    
    return dir_result_dg, dir_result_fd, mode, poly_order, RK, Nelem_dg\
    , CFL_dg, FDOA, Nelem_fd, CFL_fd, dt_, T, tt_, N_aver_fft
    
def compute_fft_DG(dir1,DG, RK, cont_comp, Nelem, Beta, Epsilon, dt_, tt_, N_aver_fft, mode):

    if mode=='CFL':
        x_,u_ = load_ensemble_data_DG_cfl(dir1, DG, RK, cont_comp, Nelem, Beta, Epsilon, dt_, tt_, N_aver_fft)
    else:
        x_,u_ = load_ensemble_data_DG(dir1, DG, RK, cont_comp, Nelem, Beta, Epsilon, dt_, tt_, N_aver_fft)
    k_freq, u_amp, KEnerg  = compute_ensemble_fft(x_,u_)

    return k_freq, u_amp, KEnerg

def compute_fft_FD(dir1,FDOA, RK, sol_comp, Nelem, dt_, tt_, N_aver_fft,mode):

    if mode=='CFL':
        x_,u_ = load_ensemble_data_FD_cfl(dir1, FDOA, RK, sol_comp, Nelem, dt_, tt_, N_aver_fft)
    else:
        x_,u_ = load_ensemble_data_FD(dir1, FDOA, RK, sol_comp, Nelem, dt_, tt_, N_aver_fft)
    k_freq, u_amp, KEnerg  = compute_ensemble_fft(x_,u_)

    return k_freq, u_amp, KEnerg
    
def plot_compare_fft_FD_DG(dir_input_fd, dir_result_fd, FDOA, fd_sol, Nelem_fd, dt_fd\
    , dir_input_dg, dir_result_dg, DG, dg_sol, Nelem_dg, Beta, Epsilon, RK, dt_dg\
    , tt_, N_aver_fft):

    # Compute Initial Spectrum:
    k_ex, E_ex, k_init, E_init = compute_exact_init_Energy()
    # Compute DG and FD Spectrums:
    
    k_freq_fd, u_amp_fd, KE_fd = compute_fft_FD(dir_input_fd, FDOA, RK, fd_sol, Nelem_fd, dt_fd, tt_, N_aver_fft)
    k_freq_dg, u_amp_dg, KE_dg = compute_fft_DG(dir_input_dg, DG, RK, dg_sol, Nelem_dg, Beta, Epsilon, dt_dg, tt_, N_aver_fft)
    
    #Plotting The Energy Spectrums:
    fig, ax = pyplot.subplots(figsize=(11, 9))
    dg_label = str('DGp') + DG + str('_RK') + RK + str(', Nelem: ') + Nelem_dg + str(', Nnodes:')+ str(size(KE_dg)*2-1)
    NN = int(Nelem_fd)+1 
    fd_label = str('FD') + FDOA + str('_RK') + RK + str(', Nnodes: ') + str(NN)
    pyplot.plot(k_freq_dg[1:], KE_dg[1:], label=dg_label)
    pyplot.plot(k_freq_fd[1:], KE_fd[1:], label=fd_label)
    
    pyplot.plot(k_init, E_init, '-k', label='Initial KE')
    pyplot.plot(k_ex, E_ex, '--k', label=r'E = const * k$^{-2}$')
    pyplot.xscale('log')
    pyplot.yscale('log')
    pyplot.grid()
    pyplot.legend()
    pyplot.xlabel('k', labelpad=10);
    pyplot.ylabel('E(k)', labelpad=10);
    pyplot.xlim(10**0, 10**3)
    pyplot.ylim(10**-12, 10 ** 0)
    fig.tight_layout()
    figname = str('./Results/') + str('EnergSpect_') + str('DGp') + DG + str('_vs_FD') +  FDOA \
    + str('_RK') + RK +str('_dt') + dt_dg+ str('_Ns') + str(N_aver_fft) + str("_") \
    + str(tt_) + str("t.png")
    pyplot.savefig(figname)
    pyplot.show()
    
def compare_fft_FD_DG2(dir_input_fd, dir_result_fd, FDOA, fd_sol, Nelem_fd, dt_fd\
, dir_input_dg, dir_input_dg1, dir_result_dg, DG, dg_sol, Nelem_dg, Beta, Epsilon, RK, dt_dg\
, tt_, N_aver_fft):

    # Compute Initial Spectrum:
    k_ex, E_ex, k_init, E_init = compute_exact_init_Energy()
    # Compute DG and FD Spectrums:
    k_freq_fd, u_amp_fd, KE_fd = compute_fft_FD(dir_input_fd, FDOA, RK, fd_sol, Nelem_fd, dt_fd, tt_, N_aver_fft)
    k_freq_dg, u_amp_dg, KE_dg = compute_fft_DG(dir_input_dg, DG, RK, dg_sol, Nelem_dg, Beta, Epsilon, dt_dg, tt_, N_aver_fft)
    k_freq_dg1, u_amp_dg, KE_dg1 = compute_fft_DG(dir_input_dg1, DG, RK, dg_sol, Nelem_dg, Beta, Epsilon, dt_dg, tt_, N_aver_fft)
    dir_input_dg2 = './input/DG/p3_RK3_test2/'
    k_freq_dg2, u_amp_dg, KE_dg2 = compute_fft_DG(dir_input_dg2, DG, RK, dg_sol, Nelem_dg, Beta, Epsilon, dt_dg, tt_, N_aver_fft)
    dir_input_dg3 = './input/DG/p3_RK3_test3/'
    k_freq_dg3, u_amp_dg, KE_dg3 = compute_fft_DG(dir_input_dg3, DG, RK, dg_sol, Nelem_dg, Beta, Epsilon, dt_dg, tt_, N_aver_fft)
    dir_input_dg4 = './input/DG/p3_RK3_test4/'
    k_freq_dg4, u_amp_dg, KE_dg4 = compute_fft_DG(dir_input_dg4, DG, RK, dg_sol, Nelem_dg, Beta, Epsilon, dt_dg, tt_, N_aver_fft)
    dir_input_dg5 = './input/DG/p3_RK3_test5/'
    k_freq_dg5, u_amp_dg, KE_dg5 = compute_fft_DG(dir_input_dg5, DG, RK, dg_sol, Nelem_dg, Beta, Epsilon, dt_dg, tt_, N_aver_fft)
    dir_input_dg = './input/DG/p3_RK3_test6/'
    k_freq_dg6, u_amp_dg, KE_dg6 = compute_fft_DG(dir_input_dg, DG, RK, dg_sol, Nelem_dg, Beta, Epsilon, dt_dg, tt_, N_aver_fft)
    dir_input_dg = './input/DG/p3_RK3_test7/'
    k_freq_dg7, u_amp_dg, KE_dg7 = compute_fft_DG(dir_input_dg, DG, RK, dg_sol, Nelem_dg, Beta, Epsilon, dt_dg, tt_, N_aver_fft)
    dir_input_dg = './input/DG/p3_RK3_test8/'
    k_freq_dg8, u_amp_dg, KE_dg8 = compute_fft_DG(dir_input_dg, DG, RK, dg_sol, Nelem_dg, Beta, Epsilon, dt_dg, tt_, N_aver_fft)
    dir_input_dg = './input/DG/p3_RK3_test9/'
    k_freq_dg9, u_amp_dg, KE_dg9 = compute_fft_DG(dir_input_dg, DG, RK, dg_sol, Nelem_dg, Beta, Epsilon, dt_dg, tt_, N_aver_fft)
    
    #Plotting The Energy Spectrums:
    fig, ax = pyplot.subplots(figsize=(11, 9))
    dg_label = str('DGp') + DG + str('_RK') + RK + str(', Nelem: ') + Nelem_dg + str(', Nnodes:')+ str(size(KE_dg)*2-1)   + str(', (P+') + str(1) + str(')pts')
    dg1_label = str('DGp') + DG + str('_RK') + RK + str(', Nelem: ') + Nelem_dg + str(', Nnodes:')+ str(size(KE_dg1)*2-1) + str(', (P+') + str(2) + str(')pts')
    dg2_label = str('DGp') + DG + str('_RK') + RK + str(', Nelem: ') + Nelem_dg + str(', Nnodes:')+ str(size(KE_dg2)*2-1) + str(', (P+') + str(3) + str(')pts')
    dg3_label = str('DGp') + DG + str('_RK') + RK + str(', Nelem: ') + Nelem_dg + str(', Nnodes:')+ str(size(KE_dg3)*2-1) + str(', (P+') + str(4) + str(')pts')
    dg4_label = str('DGp') + DG + str('_RK') + RK + str(', Nelem: ') + Nelem_dg + str(', Nnodes:')+ str(size(KE_dg4)*2-1) + str(', (P+') + str(5) + str(')pts')
    dg5_label = str('DGp') + DG + str('_RK') + RK + str(', Nelem: ') + Nelem_dg + str(', Nnodes:')+ str(size(KE_dg5)*2-1) + str(', (P+') + str(6) + str(')pts')
    dg6_label = str('DGp') + DG + str('_RK') + RK + str(', Nelem: ') + Nelem_dg + str(', Nnodes:')+ str(size(KE_dg6)*2-1) + str(', (P+') + str(9) + str(')pts')
    dg7_label = str('DGp') + DG + str('_RK') + RK + str(', Nelem: ') + Nelem_dg + str(', Nnodes:')+ str(size(KE_dg7)*2-1) + str(', (P+') + str(11) + str(')pts')
    dg8_label = str('DGp') + DG + str('_RK') + RK + str(', Nelem: ') + Nelem_dg + str(', Nnodes:')+ str(size(KE_dg8)*2-1) + str(', (P+') + str(12) + str(')pts')
    dg9_label = str('DGp') + DG + str('_RK') + RK + str(', Nelem: ') + Nelem_dg + str(', Nnodes:')+ str(size(KE_dg9)*2-1) + str(', (P+') + str(13) + str(')pts')
    NN = int(Nelem_fd)+1 
    fd_label = str('FD') + FDOA + str('_RK') + RK + str(', Nnodes: ') + str(NN)
    
    ax.plot(k_freq_dg [1:], KE_dg [1:], '-b', label=dg_label )
    ax.plot(k_freq_dg1[1:], KE_dg1[1:], label=dg1_label)
    ax.plot(k_freq_dg2[1:], KE_dg2[1:], label=dg2_label)
    ax.plot(k_freq_dg3[1:], KE_dg3[1:], label=dg3_label)
    ax.plot(k_freq_dg4[1:], KE_dg4[1:], label=dg4_label)
    ax.plot(k_freq_dg5[1:], KE_dg5[1:], label=dg5_label)
    ax.plot(k_freq_dg6[1:], KE_dg6[1:], 'c', label=dg6_label)
    ax.plot(k_freq_dg7[1:], KE_dg7[1:], label=dg7_label)
    ax.plot(k_freq_dg8[1:], KE_dg8[1:],  label=dg8_label)
    ax.plot(k_freq_dg9[1:], KE_dg9[1:], '-g', label=dg9_label)
    ax.plot(k_freq_fd [1:], KE_fd [1:], '-r', label=fd_label )
    
    ax.plot(k_init, E_init, '-k', label='Initial KE')
    ax.plot(k_ex, E_ex, '--k', label=r'E = const * k$^{-2}$')
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid()
    ax.legend()
    ax.set_xlabel('k', labelpad=10);
    ax.set_ylabel('E(k)', labelpad=10);
    ax.set_xlim(10**0, 10**4)
    ax.set_ylim(10**-12, 10 ** 0)
    fig.tight_layout()
    
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes
    inset_ax = zoomed_inset_axes(ax,
                               1.35, # zoom = 2.5
                               loc=10,
                               bbox_to_anchor=(900,750,10,10))
    inset_ax.plot(k_freq_dg[1:], KE_dg[1:], '-b', label=dg_label)
    inset_ax.plot(k_freq_dg1[1:], KE_dg1[1:], label=dg1_label)
    inset_ax.plot(k_freq_dg5[1:], KE_dg5[1:], label=dg5_label)
    inset_ax.plot(k_freq_dg6[1:], KE_dg6[1:], 'c', label=dg6_label)
    inset_ax.plot(k_freq_dg9[1:], KE_dg9[1:], '-g', label=dg9_label)
    inset_ax.plot(k_freq_fd[1:], KE_fd[1:], '-r', label=fd_label)
    
    #inset_ax.grid()
    #inset_ax.legend()
    #inset_ax.set_xlabel('k', labelpad=10);
    #inset_ax.set_ylabel('E(k)', labelpad=10);
    inset_ax.set_xscale('log')
    inset_ax.set_yscale('log')
    inset_ax.set_xlim(1.5*10**2, 6*10**2)
    inset_ax.set_ylim(10**-4, 0.5*10**-6)
    pyplot.setp( inset_ax.get_xticklabels(which='minor'), visible=False)
    pyplot.setp( inset_ax.get_yticklabels(), visible=False)
    inset_ax.set_xticklabels([])
    inset_ax.set_xticks([])
    inset_ax.set_yticks([])
    pyplot.xticks(visible=False)
    pyplot.yticks(visible=False)
    #inset_ax.minorticks_off()
    #inset_ax.tick_params(axis='both', which='both', length=0)

    
    figname = str('./Results/') + str('EnergSpect_') + str('DGp') + DG + str('_vs_FD') +  FDOA \
    + str('_RK') + RK +str('_dt') + dt_dg+ str('_Ns') + str(N_aver_fft) + str("_") \
    + str(tt_) + str("t.png")
    pyplot.savefig(figname)
    pyplot.show()
    
    
def compare_fft_FD_DG3(dir_input_fd, dir_result_fd, FDOA, fd_sol, Nelem_fd, CFL_fd\
, dir_input_dg, dir_result_dg, DG, dg_sol, Nelem_dg, Beta, Epsilon, RK, CFL_dg\
, tt_, N_aver_fft):

    mode = 'CFL';

    # Compute Initial Spectrum:
    k_ex, E_ex, k_init, E_init = compute_exact_init_Energy()
    # Compute DG and FD Spectrums:
    k_freq_fd, u_amp_fd, KE_fd = compute_fft_FD(dir_input_fd, FDOA, RK, fd_sol, Nelem_fd, CFL_fd, tt_, N_aver_fft,mode)
    k_freq_dg1, u_amp_dg, KE_dg1 = compute_fft_DG(dir_input_dg, DG, RK, dg_sol, Nelem_dg, Beta, Epsilon, CFL_dg, tt_, N_aver_fft,mode)
    dir_input_dg = './input/DG/p3_RK3_5pts/'
    k_freq_dg2, u_amp_dg, KE_dg2 = compute_fft_DG(dir_input_dg, DG, RK, dg_sol, Nelem_dg, Beta, Epsilon, CFL_dg, tt_, N_aver_fft,mode)
    dir_input_dg = './input/DG/p3_RK3_8pts/'
    k_freq_dg3, u_amp_dg, KE_dg3 = compute_fft_DG(dir_input_dg, DG, RK, dg_sol, Nelem_dg, Beta, Epsilon, CFL_dg, tt_, N_aver_fft,mode)
    dir_input_dg = './input/DG/p3_RK3_15pts/'
    k_freq_dg4, u_amp_dg, KE_dg4 = compute_fft_DG(dir_input_dg, DG, RK, dg_sol, Nelem_dg, Beta, Epsilon, CFL_dg, tt_, N_aver_fft,mode)
   
    
    #Plotting The Energy Spectrums:
    fig, ax = pyplot.subplots(figsize=(11, 9))
    dg1_label = str('DGp') + DG + str('_RK') + RK + str(', Nelem: ') + Nelem_dg + str(', Nnodes:')+ str(size(KE_dg1)*2-1)   + str(', (P+') + str(1) + str(')pts')
    dg2_label = str('DGp') + DG + str('_RK') + RK + str(', Nelem: ') + Nelem_dg + str(', Nnodes:')+ str(size(KE_dg2)*2-1) + str(', (P+') + str(2) + str(')pts')
    dg3_label = str('DGp') + DG + str('_RK') + RK + str(', Nelem: ') + Nelem_dg + str(', Nnodes:')+ str(size(KE_dg3)*2-1) + str(', (P+') + str(5) + str(')pts')
    dg4_label = str('DGp') + DG + str('_RK') + RK + str(', Nelem: ') + Nelem_dg + str(', Nnodes:')+ str(size(KE_dg4)*2-1) + str(', (P+') + str(13) + str(')pts')
    NN = int(Nelem_fd)+1 
    fd_label = str('FD') + FDOA + str('_RK') + RK + str(', Nnodes: ') + str(NN)
    
    ax.plot(k_freq_dg1[1:], KE_dg1[1:], label=dg1_label)
    ax.plot(k_freq_dg2[1:], KE_dg2[1:], label=dg2_label)
    ax.plot(k_freq_dg3[1:], KE_dg3[1:], label=dg3_label)
    ax.plot(k_freq_dg4[1:], KE_dg4[1:], label=dg4_label)
    ax.plot(k_freq_fd [1:], KE_fd [1:], label=fd_label )
    
    ax.plot(k_init, E_init, '-k', label='Initial KE')
    ax.plot(k_ex, E_ex, '--k', label=r'E = const * k$^{-2}$')
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid()
    ax.legend()
    ax.set_xlabel('k', labelpad=10);
    ax.set_ylabel('E(k)', labelpad=10);
    ax.set_xlim(10**0, 10**4)
    ax.set_ylim(10**-12, 10 ** 0)
    fig.tight_layout()
    
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes
    inset_ax = zoomed_inset_axes(ax,
                               1.35, # zoom = 2.5
                               loc=10,
                               bbox_to_anchor=(900,750,10,10))

    inset_ax.plot(k_freq_dg1[1:], KE_dg1[1:], label=dg1_label)
    inset_ax.plot(k_freq_dg2[1:], KE_dg2[1:], label=dg2_label)
    inset_ax.plot(k_freq_dg3[1:], KE_dg3[1:], label=dg3_label)
    inset_ax.plot(k_freq_dg4[1:], KE_dg4[1:], label=dg4_label)
    inset_ax.plot(k_freq_fd[1:], KE_fd[1:], label=fd_label)
    
    #inset_ax.grid()
    #inset_ax.legend()
    #inset_ax.set_xlabel('k', labelpad=10);
    #inset_ax.set_ylabel('E(k)', labelpad=10);
    inset_ax.set_xscale('log')
    inset_ax.set_yscale('log')
    inset_ax.set_xlim(1.5*10**2, 6*10**2)
    inset_ax.set_ylim(10**-4, 0.5*10**-6)
    pyplot.setp( inset_ax.get_xticklabels(which='minor'), visible=False)
    pyplot.setp( inset_ax.get_yticklabels(), visible=False)
    inset_ax.set_xticklabels([])
    inset_ax.set_xticks([])
    inset_ax.set_yticks([])
    pyplot.xticks(visible=False)
    pyplot.yticks(visible=False)
    #inset_ax.minorticks_off()
    #inset_ax.tick_params(axis='both', which='both', length=0)

    
    figname = str('./Results/') + str('EnergSpect_') + str('DGp') + DG + str('_vs_FD') +  FDOA \
    + str('_RK') + RK +str('_CFL') + str(CFL_dg)+ str('_Ns') + str(N_aver_fft) + str("_") \
    + str(tt_) + str("t.png")
    pyplot.savefig(figname)
    pyplot.show()
    
def compare_fft_FD_DG4(dir_input_fd, dir_result_fd, FDOA, fd_sol, Nelem_fd, dt_fd\
, dir_input_dg, dir_result_dg, DG, dg_sol, Nelem_dg, Beta, Epsilon, RK, dt_dg\
, tt_, N_aver_fft):

    mode = 'dt';

    # Compute Initial Spectrum:
    k_ex, E_ex, k_init, E_init = compute_exact_init_Energy()
    # Compute DG and FD Spectrums:
    k_freq_fd, u_amp_fd, KE_fd = compute_fft_FD(dir_input_fd, FDOA, RK, fd_sol, Nelem_fd, dt_fd, tt_, N_aver_fft,mode)
    k_freq_dg1, u_amp_dg, KE_dg1 = compute_fft_DG(dir_input_dg, DG, RK, dg_sol, Nelem_dg, Beta, Epsilon, dt_dg, tt_, N_aver_fft,mode)
    dir_input_dg = './input/DG/p3_RK3_5pts/'
    k_freq_dg2, u_amp_dg, KE_dg2 = compute_fft_DG(dir_input_dg, DG, RK, dg_sol, Nelem_dg, Beta, Epsilon, dt_dg, tt_, N_aver_fft,mode)
    dir_input_dg = './input/DG/p3_RK3_8pts/'
    k_freq_dg3, u_amp_dg, KE_dg3 = compute_fft_DG(dir_input_dg, DG, RK, dg_sol, Nelem_dg, Beta, Epsilon, dt_dg, tt_, N_aver_fft,mode)
    dir_input_dg = './input/DG/p3_RK3_15pts/'
    k_freq_dg4, u_amp_dg, KE_dg4 = compute_fft_DG(dir_input_dg, DG, RK, dg_sol, Nelem_dg, Beta, Epsilon, dt_dg, tt_, N_aver_fft,mode)
   
    
    #Plotting The Energy Spectrums:
    fig, ax = pyplot.subplots(figsize=(11, 9))
    dg1_label = str('DGp') + DG + str('_RK') + RK + str(', Nelem: ') + Nelem_dg + str(', Nnodes:')+ str(size(KE_dg1)*2-1)   + str(', (P+') + str(1) + str(')pts')
    dg2_label = str('DGp') + DG + str('_RK') + RK + str(', Nelem: ') + Nelem_dg + str(', Nnodes:')+ str(size(KE_dg2)*2-1) + str(', (P+') + str(2) + str(')pts')
    dg3_label = str('DGp') + DG + str('_RK') + RK + str(', Nelem: ') + Nelem_dg + str(', Nnodes:')+ str(size(KE_dg3)*2-1) + str(', (P+') + str(5) + str(')pts')
    dg4_label = str('DGp') + DG + str('_RK') + RK + str(', Nelem: ') + Nelem_dg + str(', Nnodes:')+ str(size(KE_dg4)*2-1) + str(', (P+') + str(13) + str(')pts')
    NN = int(Nelem_fd)+1 
    fd_label = str('FD') + FDOA + str('_RK') + RK + str(', Nnodes: ') + str(NN)
    
    ax.plot(k_freq_dg1[1:], KE_dg1[1:], label=dg1_label)
    ax.plot(k_freq_dg2[1:], KE_dg2[1:], label=dg2_label)
    ax.plot(k_freq_dg3[1:], KE_dg3[1:], label=dg3_label)
    ax.plot(k_freq_dg4[1:], KE_dg4[1:], label=dg4_label)
    ax.plot(k_freq_fd [1:], KE_fd [1:], label=fd_label )
    
    ax.plot(k_init, E_init, '-k', label='Initial KE')
    ax.plot(k_ex, E_ex, '--k', label=r'E = const * k$^{-2}$')
    
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.grid()
    ax.legend()
    ax.set_xlabel('k', labelpad=10);
    ax.set_ylabel('E(k)', labelpad=10);
    ax.set_xlim(10**0, 10**4)
    ax.set_ylim(10**-12, 10 ** 0)
    fig.tight_layout()
    
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes
    inset_ax = zoomed_inset_axes(ax,
                               1.35, # zoom = 2.5
                               loc=10,
                               bbox_to_anchor=(900,750,10,10))

    inset_ax.plot(k_freq_dg1[1:], KE_dg1[1:], label=dg1_label)
    inset_ax.plot(k_freq_dg2[1:], KE_dg2[1:], label=dg2_label)
    inset_ax.plot(k_freq_dg3[1:], KE_dg3[1:], label=dg3_label)
    inset_ax.plot(k_freq_dg4[1:], KE_dg4[1:], label=dg4_label)
    inset_ax.plot(k_freq_fd[1:], KE_fd[1:], label=fd_label)
    
    #inset_ax.grid()
    #inset_ax.legend()
    #inset_ax.set_xlabel('k', labelpad=10);
    #inset_ax.set_ylabel('E(k)', labelpad=10);
    inset_ax.set_xscale('log')
    inset_ax.set_yscale('log')
    inset_ax.set_xlim(1.5*10**2, 6*10**2)
    inset_ax.set_ylim(10**-4, 0.5*10**-6)
    pyplot.setp( inset_ax.get_xticklabels(which='minor'), visible=False)
    pyplot.setp( inset_ax.get_yticklabels(), visible=False)
    inset_ax.set_xticklabels([])
    inset_ax.set_xticks([])
    inset_ax.set_yticks([])
    pyplot.xticks(visible=False)
    pyplot.yticks(visible=False)
    #inset_ax.minorticks_off()
    #inset_ax.tick_params(axis='both', which='both', length=0)

    figname = str('./Results/') + str('EnergSpect_') + str('DGp') + DG + str('_vs_FD') +  FDOA \
    + str('_RK') + RK +str('_dt') + dt_dg + str('_Ns') + str(N_aver_fft) + str("_") \
    + str(tt_) + str("t.png")
    pyplot.savefig(figname)
    pyplot.show()
    
    
#===============================================================================================================#
 
 
















