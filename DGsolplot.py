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

def plot_advect(mode, DG, RK, CFL, Nelem, T, dt_\
                , Beta, dir1, aver, nodal_exact, nodal_comp, discont ):

    Beta = Decimal(Beta.quantize(Decimal('.01')));

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
                + str("_CFL") + str(CFL) + str("_Beta") + str(Beta) \
                + str("_") + str(T) + str("T.dat")
    elif mode == 'dt_const':
        fname = dir1 + discont + str("_N") + Nelem \
                + str("_dt") + dt_ + str("_Beta") + str(Beta) \
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

    pyplot.plot(xn_exact, un_exact, '--k', label='Exact sol');

    for i in range(0, size(x_disc) - 1, Np):
        xx = x_disc[i:i + Np];
        uu = u_disc[i:i + Np];
        pyplot.plot(xx, uu, '-r');

    ylim_0.append(min(u_disc));
    ylim_1.append(max(u_disc));

    xx = x_disc[i:i + Np];
    uu = u_disc[i:i + Np];
    ll = str("Numerical sol, (") + r'$\beta= $' + str(Beta) + " )";
    pyplot.plot(xx, uu, '-r', label=ll);

    pyplot.legend();

    title_a = str("DGp") + DG + " RK" + RK \
              + ", and upwind_param (" + r'$\beta= $' \
              + str(Beta) + ")\n for CFL=" + str(CFL) \
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