from matplotlib import pyplot ,ticker   #and the useful plotting library
from numpy import sin,cos,pi,linspace,ones,zeros,abs,min,max,exp, shape, empty_like , size, loadtxt, arange , log
from decimal import Decimal
import csv
from matplotlib import pyplot as plt

import numpy as np

from fft_toolbox_python import load_data, compute_fft, compute_Etotal

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


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#                                  D I F F U S I ON
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
def plot_diffus(dir1, mode, DG, RK, CFL, Nelem, N_disc_ppt, tt_, dt_ \
                     , Epsilon, gamma_, diffus_scheme, cont_num, disc_num, T):

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
            + mm_name + str("_Eps") + str(Epsilon)\
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
            + mm_name + str("_Eps") + str(Epsilon)\
            + str('_')+ str(tt_) + str("t.dat")
    data = loadtxt(fname)
    x_disc = data[:, 0]
    u_disc = data[:, 1]
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
    print('error: ',abs(1-max(u_disc)))
    print('error_cont: ',abs(max(u_cont_exact)-max(u_cont)))
    print('error_disc: ',abs(max(u_disc_exact)-max(u_disc)))

    #=========================== PLOTTING Solution(1) ============================#
    fig = pyplot.figure();
    #plotting continuous data:
    pyplot.plot(x_cont_exact, u_cont_exact, '-k', label='Exact solution');
    label_cont = str("DGp")+ str(DG) + r'-$\eta$' + str(Epsilon) \
    +"_RK" + str(RK) +", CFL="+str(CFL)+", t=" + str(tt_); # discontinuous label
    pyplot.plot(x_cont, u_cont, '--b', label=label_cont,lw=1.2);
    
    ylim_0 = list();
    ylim_1 = list();
    ylim_0.append(min(u_cont_exact));
    ylim_1.append(max(u_cont_exact));

    #plotting discontinuous numerical data:
    nn = size(x_disc)
    Np = N_disc_ppt
    #markers_on=[0,5,10,15,19]
    markers_on=1;
    for i in range(0, size(x_disc) - 1, Np):
        xx = x_disc[i:i + Np];
        uu = u_disc[i:i + Np];
        pyplot.plot(xx, uu, ':om',lw=0.7,markevery=markers_on);

    ylim_0.append(min(u_disc));
    ylim_1.append(max(u_disc));

    xx = x_disc[i:i + Np];
    uu = u_disc[i:i + Np];
    pyplot.plot(xx, uu, ':om', markevery=markers_on, label='Discontinuous DG solution',lw=0.7);

    pyplot.legend();

    #pyplot.title(title_a);
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
    '_eta'+str(Epsilon)+'_N'+Nelem+mm_name+str('_')+\
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
    label_disc = str("DGp")+ str(DG) + r'-$\eta$' + str(Epsilon) \
    +"_RK" + str(RK) +", CFL="+str(CFL)+", t=" + str(tt_); # discontinuous label
    for i in range(0, size(x_disc) - 1, Np):
        xx = x_disc[i:i + Np];
        uu = u_disc[i:i + Np];
        pyplot.plot(xx, uu, ':m');
    xx = x_disc[i:i + Np];
    uu = u_disc[i:i + Np];
    pyplot.plot(xx, uu, ':m', label=label_disc);
    
    #plotting discontinuous exact data:
    nn = size(x_disc_exact)
    Np = N_disc_ppt
    label_disc = "exact discontinuous" ; # discontinuous exact label
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
    '_eta'+str(Epsilon)+'_N'+Nelem+mm_name+str('_')+\
    str(tt_)+str('t_disc.png')
    fig.set_size_inches(15.0, 9.0, forward=True)
    pyplot.savefig(figname)
    
    #================== Plot FFT ==========================#
    label_dg = str("DGp")+ str(DG) + r'-$\eta$' + str(Epsilon) \
    +"_RK" + str(RK) +", CFL="+str(CFL)+", t=" + str(tt_); # discontinuous label
    k_max_ = int(k_freq[-1]);
    print('k_max: ',k_max_, '  k_max_ex: ',int(k_freq_exact[-1]))
    print('exact Nelem:',size(x_cont_exact))
    fig, ax = plt.subplots()
    
    plt.plot(2*pi*k_freq_exact/(3*int(Nelem)), u_amp_exact, '-k',markevery=1\
        , label=r'Exact, E$_{tot}$= '+str(np.round(E_ex,4)))
    plt.plot(2*pi*k_freq/(3*int(Nelem)), u_amp, ':om',markevery=1, \
        label=r'Numerical DG, E$_{tot}$='+str(np.round(E_tot,4)))
    
    xlabels = ['0',r'$\pi$/8',r'$\pi$/4',r'$\pi$/3',r'$\pi$/2',r'3$\pi$/4',r'$\pi$'];
    xlocs = [0,pi/8,pi/4,pi/3,pi/2,3*pi/4,pi];
    plt.xticks(xlocs, xlabels);
    plt.xlim(0,pi)
    plt.legend();
    plt.xlabel('K/(P+1)', labelpad=2);
    plt.ylabel('|u|', labelpad=2);
    plt.grid()

    fig.set_size_inches(13.0, 10.0, forward=True)
    fig.tight_layout(pad=0.0, w_pad=10.0, h_pad=10.0)
    
    figname = dir1 + 'tempfig/'+ 'fft_N'+Nelem\
    +mm_name+str('_')+str(tt_)+str('t_comp.png')
    fig.set_size_inches(15.0, 9.0, forward=True)
    plt.savefig(figname,format='png')
    figname = dir1 + 'tempfig/'+ 'fft_N'+Nelem\
    +mm_name+str('_')+str(tt_)+str('t_comp.eps')
    plt.savefig(figname,format='eps')
    
    #=========================== PLOTTING Errors Evolution ============================#
    # Read error data:
    res_dir = './Results/'
    fname = dir1 + 'errors/errors_N'+str(Nelem)+'_CFL'+str(CFL)+'_Eps'+str(Epsilon)+'_'+str(T)+'T.dat'
    data = loadtxt(fname);  # continuous exact nodal solution
    time  = data[:, 0];
    L1err = data[:, 1];
    L2err = data[:, 2];
    L1err_nodal = data[:, 3];
    L2err_nodal = data[:, 4];

    #plotting
    fig = pyplot.figure();
    
    tau_p = float(gamma_)*time*((int(DG)+1)*int(Nelem))**2
    pyplot.plot(tau_p,L1err,'-.sb',label=r'$L_{1}(u(\xi))$')
    pyplot.plot(tau_p, L2err,'-ok',label=r'$L_{2}(u(\xi))$')
    pyplot.plot(tau_p,L1err_nodal,':vm',label=r'$L_{1}$, nodal')
    pyplot.plot(tau_p, L2err_nodal,'--c',label=r'$L_{2}$, nodal')
    
    pyplot.xlabel('time')
    pyplot.ylabel('errors')
    pyplot.legend();
    #pyplot.xscale('log')
    #pyplot.yscale('log')
    
    figname = dir1 + 'tempfig/' + 'errors_p'+str(DG)+'RK'+RK+\
    '_eta'+str(Epsilon)+'_N'+Nelem+mm_name+str('_')+\
    str(T)+str('T.png')
    fig.set_size_inches(15.0, 9.0, forward=True)
    pyplot.savefig(figname)
    
    #------------------------
    pyplot.show()

    return 'true'

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#                                  A D V E C T I O N
#
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

def plot_advec(dir1, mode, DG, RK, CFL, Nelem, N_disc_ppt, tt_, dt_ \
                     , Beta, Epsilon, cont_num, disc_num, T_period):

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
    #Exact:
    fname = dir1 + 'time_data/u_cont_exact_'+ str(tt_) + str("t.dat");
    data = loadtxt(fname);  # continuous exact nodal solution
    x_cont_exact = data[:, 0];
    u_cont_exact = data[:, 1];
    del fname, data
    
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
 
def plot_burgers_decay_turb(dir_input, mode, DG, RK, CFL, Nelem, tt_, dt_, Beta, Epsilon, gamma_, cont_num, disc_num):
    
    if  (mode == 'test') | (mode == 'dt_const'):
        mm_name = str('_dt') + dt_
        m_name = str('dt=') + dt_
    elif(mode == 'CFL_const') | (mode == 'normal'):
        mm_name = str('_CFL')+ str(CFL)     
        m_name = str('CFL=')+ str(CFL)             

    dt = float(dt_);
    
    print('dir_input: ',dir_input)

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
    
    fname = dir_input + disc_num + str("_N") + Nelem \
    + mm_name + str("_Beta") + str(Beta)\
    + str("_Eps") + str(Epsilon) \
    + str('_')+ str(tt_) + str("t.dat")

    data = loadtxt(fname)

    x_disc = data[:, 0]
    u_disc = data[:, 1]

    nn = size(x_disc);

    if int(DG) == 0:
        Np = 2;
    elif int(DG) == 1:
        Np = 2;
    else:
        Np = 10;

    ################ Plotting the solution #####################
    fig, ax = pyplot.subplots(frameon='True')

    ylim_0 = list();
    ylim_1 = list();
    ylim_0.append(min(u_cont));
    ylim_1.append(max(u_cont));

    pyplot.plot(x_cont, u_cont, '-k', label='continuous sol');

    for i in range(0, size(x_disc) - 1, Np):
        xx = x_disc[i:i + Np];
        uu = u_disc[i:i + Np];
        #pyplot.plot(xx, uu, '-r');

    ylim_0.append(min(u_disc));
    ylim_1.append(max(u_disc));

    xx = x_disc[i:i + Np];
    uu = u_disc[i:i + Np];
    ll = str("discontinuous sol");
    #pyplot.plot(xx, uu, '-r', label=ll);

    #pyplot.legend();

    title_a = str("DGp") + DG  + "-RK" + RK +" with " + m_name \
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
    
    temp_name = 'sol_vs_x_p' + DG + 'RK' + RK +'_Ne'+ str(Nelem) +'_'+ m_name \
              + '_t'+ str(tt_);
           
    figname = dir_input + str('/tempfig/') + temp_name +'.eps'
    fig.savefig(figname,format='eps',bbox='tight')
    figname = dir_input + str('/tempfig/') + temp_name +'.png'
    plt.savefig(figname,format='png',bbox='tight')
    
    ###################### Plot the fft ##########################
    fig, ax = pyplot.subplots(frameon='True')
    ax.plot(k_freq, KEnerg)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.legend(fontsize=19,edgecolor='black')
    ax.set_facecolor('white')
    ax.set_xlabel(r'$k$', labelpad=2,fontsize=24);
    ax.set_ylabel(r'$E$', labelpad=2, fontsize=24);
    ax.set_xlim(10**0, 10**4)
    ax.set_ylim(10**-10, 10 ** -1)
    ax.tick_params(axis='both', which='both', labelsize=24)
    
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    ax.spines["bottom"].set_visible(True)
    
    fig.set_size_inches(13.0, 9.0, forward=True)
    fig.tight_layout(pad=0, w_pad=10.0, h_pad=10.0,rect=(0.0,0.0,1.0,0.985))
    
    temp_name = 'KE_p' + DG + 'RK' + RK +'_Ne'+ str(Nelem)+ '_' + m_name \
              + '_t'+ str(tt_);
    figname = dir_input + str('/tempfig/') + temp_name+'.eps'
    fig.savefig(figname,format='eps',bbox='tight')
    figname = dir_input + str('/tempfig/') + temp_name+'.png'
    plt.savefig(figname,format='png',bbox='tight')

    pyplot.show()

    return 'true'
                     
                     
