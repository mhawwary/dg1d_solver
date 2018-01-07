from matplotlib import pyplot ,ticker   #and the useful plotting library
from numpy import sin,cos,pi,linspace,ones,zeros,abs,min,max,exp, shape, empty_like , size, loadtxt, arange , log
from decimal import Decimal
import csv

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

def plot_diffus(diffus_scheme, mode, DG, RK, CFL, Nelem, T, dt_\
                , Epsilon, dir1, aver, nodal_exact, nodal_comp, discont ):

    Epsilon = Decimal(Epsilon.quantize(Decimal('.01')));

    CFL = Decimal(CFL.quantize(Decimal('.001')));
    T = Decimal(T.quantize(Decimal('.001')));
    dt = float(dt_);

    fname = dir1 + nodal_exact + str("_") + str(T) + str("T.dat");
    data = loadtxt(fname);  # continuous exact nodal solution
    xn_exact = data[:, 0];
    un_exact = data[:, 1];

    del fname, data

    # fname = dir1+nodal_comp+str("_N")+Nelem\
    # +str("_CFL")+str(CFL)+str("_Beta")+str(Beta)\
    # +str("_")+str(T)+str("T.dat");
    # data = loadtxt(fname);   # continuous numerical nodal solution
    # xn_comp = data[:,0];
    # un_comp = data[:,1];
    # del fname,data

    if mode == 'CFL_const' or mode == 'normal' or mode == 'test':
        fname = dir1 + discont + str("_N") + Nelem \
                + str("_CFL") + str(CFL) + str("_Eps") + str(Epsilon) \
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

    pyplot.show()

    return 'true'

def plot_advec(dir1, mode, DG, RK, CFL, Nelem, N_disc_ppt, tt_, dt_ \
                     , Beta, Epsilon, cont_num, disc_num):

    dt = float(dt_);
    
    if  (mode == 'test') | (mode == 'dt_const'):
        mm_name = str('_dt') + dt_
        m_name = str('dt=') + dt_
    elif(mode == 'CFL_const') | (mode == 'normal'):
        mm_name = str('_CFL')+ str(CFL)     
        m_name = str('CFL=')+ str(CFL)             


    fname = dir1 + 'time_data/u_cont_exact_'+ str(tt_) + str("t.dat");
    data = loadtxt(fname);  # continuous exact nodal solution
    xn_exact = data[:, 0];
    un_exact = data[:, 1];
    del fname, data
    
    fname = dir1 + cont_num + str("_N") + Nelem \
            + mm_name + str("_Beta") + str(Beta)\
            + str('_')+ str(tt_) + str("t.dat")

    data = loadtxt(fname)

    x_cont = data[:, 0]
    u_cont = data[:, 1]
    del fname, data

    fname = dir1 + disc_num + str("_N") + Nelem \
            + mm_name + str("_Beta") + str(Beta)\
            + str('_')+ str(tt_) + str("t.dat")

    data = loadtxt(fname)

    x_disc = data[:, 0]
    u_disc = data[:, 1]
    
    print('u_cont_max: ',max(u_cont))
    print('u_disc_max: ',max(u_disc))
    print('error: ',abs(1-max(u_disc)))

    del fname, data

    nn = size(x_disc);
        
    Np = N_disc_ppt;

    fig = pyplot.figure();

    ylim_0 = list();
    ylim_1 = list();
    ylim_0.append(min(un_exact));
    ylim_1.append(max(un_exact));

    pyplot.plot(xn_exact, un_exact, '--k', label='exact sol');
    pyplot.plot(x_cont, u_cont, '-.b', label='cont sol');
    

    for i in range(0, size(x_disc) - 1, Np):
        xx = x_disc[i:i + Np];
        uu = u_disc[i:i + Np];
        pyplot.plot(xx, uu, '-r');

    ylim_0.append(min(u_disc));
    ylim_1.append(max(u_disc));

    xx = x_disc[i:i + Np];
    uu = u_disc[i:i + Np];
    ll = str("DGp")+ str(DG) + r'-$\beta$' + str(Beta) \
    +"_RK" + str(RK) +", CFL="+str(CFL)+", t=" + str(tt_);
    pyplot.plot(xx, uu, '-r', label=ll);

    pyplot.legend();

    #title_a = str("DGp") + DG + " RK" + RK \
     #         + ", and upwind_param (" + r'$\beta= $' \
      #        + str(Beta) + ")\n for CFL=" + str(CFL) \
       #       + " and at t/T=" + str(T);

    #pyplot.title(title_a);
    pyplot.xlabel('X', labelpad=10);
    pyplot.ylabel('u(x)', labelpad=10);

    pyplot.xlim(min(xn_exact), max(xn_exact));
    pyplot.ylim(min(ylim_0) * 1.05, max(ylim_1) * 1.05);

    n_divisions = 8;
    xtick_dx = (xn_exact[-1] - xn_exact[0] ) / n_divisions;
    xlabels = arange(xn_exact[0], xn_exact[-1]+xtick_dx,xtick_dx);
    #xlabels = [0,10,20,30,40,50,60,70,80];
    xlocs = xlabels;
    pyplot.xticks(xlocs, xlabels);
    pyplot.grid()
    pyplot.xlim(min(xn_exact), max(xn_exact))
    pyplot.ylim(min(ylim_0), max(ylim_1)*1.3)

    fig.tight_layout()
    figname = dir1 + 'tempfig/' + 'p'+str(DG)+'RK'+RK+\
    '_beta'+str(Beta)+'_N'+Nelem+mm_name+str('_')+\
    str(tt_)+str('t.png')
    fig.set_size_inches(15.0, 9.0, forward=True)
    pyplot.savefig(figname)
    pyplot.show()

    return 'true'

def plot_AdvecDiffus(diffus_scheme, mode, DG, RK, CFL, Nelem, T, dt_\
                , Beta, Epsilon, dir1, aver, nodal_exact, nodal_comp, discont ):

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
    
 #==================================================================================================================#
 #                                Decay Burger's Turbulence 
 #==================================================================================================================#
 
def plot_burgers_decay_turb(dir_input, mode, DG, RK, CFL, Nelem, tt_, dt_, Beta, Epsilon, cont_num, disc_num):
    
    if  (mode == 'test') | (mode == 'dt_const'):
        mm_name = str('_dt') + dt_
        m_name = str('dt=') + dt_
    elif(mode == 'CFL_const') | (mode == 'normal'):
        mm_name = str('_CFL')+ str(CFL)     
        m_name = str('CFL=')+ str(CFL)             

    dt = float(dt_);
    
    #print('Epsilon: ',Epsilon)

    fname = dir_input + cont_num + str("_N") + Nelem \
    + mm_name + str("_Beta") + str(Beta)\
    + str("_Eps") + str(Epsilon) \
    + str('_')+ str(tt_) + str("t.dat")

    data = loadtxt(fname)

    x_cont = data[:, 0]
    u_cont = data[:, 1]

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

    pyplot.figure();

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

    pyplot.show()

    return 'true'
                     
                     
