import numpy                 #loading our favorite library
from matplotlib import pyplot    #and the useful plotting library
from numpy import sin,cos,pi,linspace,ones,zeros,abs,min,max,exp, shape, empty_like , size , arange, loadtxt
import argparse
from decimal import Decimal
import csv


pyplot.rc('legend',**{'loc':'upper left'});
pyplot.rcParams[u'legend.fontsize'] = 16
pyplot.rcParams[u'legend.edgecolor']='white'
pyplot.rcParams[u'legend.facecolor']='0.8'
pyplot.rcParams[u'font.weight']='normal'
#pyplot.rcParams['font.serif']='false'
pyplot.rcParams[u'xtick.labelsize']=15
pyplot.rcParams[u'ytick.labelsize']=15
pyplot.rcParams[u'axes.titlesize']=16
pyplot.rcParams[u'axes.labelsize']=16
pyplot.rcParams[u'axes.spines.right']='false';
pyplot.rcParams[u'axes.spines.top']='false';
pyplot.rcParams[u'lines.linewidth'] = 1.5;
pyplot.rcParams[u'lines.markersize'] = 8;

parser = argparse.ArgumentParser(description='python_DG_argument_parsing');

parser.add_argument('-f', type=str, dest='python_input');

args = parser.parse_args();

Beta = list();
CFL = list();

with open(args.python_input) as file: 
    
    reader=csv.reader(file, delimiter=':');

    for row in reader:
        if row[0]=='mode':
            mode=str(row[1]);
        elif row[0]=='CFL':
            i=0;
            for i,j in enumerate(row[1:]):
                CFL.append(float(j));
        elif row[0]=='DGp':
            DG=str(row[1]);
        elif row[0] == 'RK':    
            RK=str(row[1]);
        elif row[0] == 'Nelem':    
            Nelem=str(int(row[1]));
        elif row[0] == 'T':     
            T=Decimal(row[1]);
        elif row[0] =='dir':    
            dir1 = str(row[1]);
        elif row[0] =='aver':   
            aver = str(row[1]);
        elif row[0] == 'cont_exact':
            nodal_exact = str(row[1]);
        elif row[0] == 'cont_num':
            nodal_comp = str(row[1]);
        elif row[0] == 'discont': 
            discont = str(row[1]);
        elif row[0]=='Beta':
            i=0;
            for i,j in enumerate(row[1:]):
                Beta.append(float(j));

for ii in range(0,size(Beta)):
    Beta[ii] = Decimal(Beta[ii]);
    Beta[ii] = Decimal(Beta[ii].quantize(Decimal('.01')));
    CFL[ii] = Decimal(CFL[ii]);
    CFL[ii] = Decimal(CFL[ii].quantize(Decimal('.001')));

T=Decimal(T.quantize(Decimal('.001')));

fname_un_ex = dir1+nodal_exact+str("_")+str(T)+str("T.dat");
data_exact_nodal= loadtxt(fname_un_ex);  # continuous exact nodal solution
xn_exact = data_exact_nodal[:,0]; 
un_exact= data_exact_nodal[:,1];

fname_un_disc = dir1+discont+str("_N")+Nelem\
+str("_CFL")+str(CFL[0])+str("_Beta")+str(Beta[0])\
+str("_")+str(T)+str("T.dat");
data_disc= loadtxt(fname_un_disc);

x_disc = data_disc[:,0];
u_disc = data_disc[:,1];

fname_un_disc_upw = dir1+discont+str("_N")+Nelem\
+str("_CFL")+str(CFL[1])+str("_Beta")+str(Beta[1])\
+str("_")+str(T)+str("T.dat");
data_disc_upw= loadtxt(fname_un_disc_upw);

x_disc_upw = data_disc_upw[:,0];
u_disc_upw = data_disc_upw[:,1];

#! Plotting The figure:
#---------------------------
if int(DG)<=1:
    Np = 2;
else:
    Np=10;

fig, ax = pyplot.subplots(figsize=(10.5, 7.5)) ;
pyplot.plot(xn_exact,un_exact,'--k',label='exact sol');

ylim_0 = list();
ylim_1 = list();

ylim_0.append(min(un_exact));
ylim_1.append(max(un_exact));

ylim_0.append(min(u_disc));
ylim_1.append(max(u_disc));

ylim_0.append(min(u_disc_upw));
ylim_1.append(max(u_disc_upw));

nn = size(x_disc);

for i in range(0,size(x_disc)-1,Np):
    xx = x_disc[i:i+Np];
    uu = u_disc[i:i+Np];
    pyplot.plot(xx,uu,'-r');
    xx = x_disc_upw[i:i + Np];
    uu = u_disc_upw[i:i + Np];
    pyplot.plot(xx, uu, '-b');

xx = x_disc[i:i+Np];
uu = u_disc[i:i+Np];
ll = r'$\beta= $'+str(Beta[0])+', CFL='+str(CFL[0]);
pyplot.plot(xx,uu,'-r',label=ll);

xx = x_disc_upw[i:i+Np];
uu = u_disc_upw[i:i+Np];
ll = r'$\beta= $'+str(Beta[1])+', CFL='+str(CFL[1]);
pyplot.plot(xx,uu,'-b',label=ll);

pyplot.legend(loc='upper left');
#CFL=Decimal(CFL.quantize(Decimal('.1')));
title_a = str("DGp")+ DG + " RK"+ RK +" and at t/T="+str(T);
#pyplot.title(title_a);

pyplot.xlabel('X',labelpad=10);
pyplot.ylabel('u',labelpad=10);

pyplot.xlim(-1.0,1.0);
pyplot.ylim(min(ylim_0)*1.05,max(ylim_1)*1.05);

from matplotlib.ticker import FormatStrFormatter
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

xlabels=[-1.0,-0.5,0,0.5,1.0];
xlocs=[-1.0,-0.5,0,0.5,1.0];
pyplot.xticks(xlocs, xlabels);

ylabels=arange(-0.2,1.2,0.2);
ylocs=arange(-0.2,1.2,0.2);
pyplot.yticks(ylocs,ylabels);

#pyplot.savefig('fig1.png',bbox='tight')
fig.tight_layout()
pyplot.show()












