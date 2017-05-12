import numpy                 #loading our favorite library
from matplotlib import pyplot    #and the useful plotting library
from numpy import sin,cos,pi,linspace,ones,zeros,abs,min,max,exp, shape, empty_like , size , arange
import argparse
from decimal import Decimal
import csv

pyplot.rc('legend',**{'loc':'upper left'});
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

parser = argparse.ArgumentParser(description='python_DG_argument_parsing');

parser.add_argument('-f', type=str, dest='python_input');

args = parser.parse_args();

Beta = list();

with open(args.python_input) as file: 
    
    reader=csv.reader(file, delimiter=':');

    for row in reader:
        if row[0]=='CFL':    
            CFL=Decimal(row[1]);
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
        elif row[0] == 'nodal_exact': 
            nodal_exact = str(row[1]);
        elif row[0] == 'nodal_num': 
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

CFL=Decimal(CFL.quantize(Decimal('.01')));
T=Decimal(T.quantize(Decimal('.1')));

fname_un_ex = dir1+nodal_exact+str("_N")+Nelem\
+str("_CFL")+str(CFL)+str("_Beta")\
+str(Beta[0])+str("_")+str(T)+str("T.dat");

data_exact_nodal= numpy.loadtxt(fname_un_ex);
xn_exact = data_exact_nodal[:,0]; 
un_exact= data_exact_nodal[:,1]; 

import matplotlib.colors as colors
import matplotlib.cm as cmx

jet = cm = pyplot.get_cmap('jet') 
cNorm  = colors.Normalize(vmin=float(Beta[0]), vmax=float(Beta[-1]))
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet);

#cluster = ['o','^','v','p','s'];
#lines_ = ['-','--','-','--','-'];

#! Plotting The figure:
#---------------------------
if int(DG)<=1:
    Np = 2;
else:
    Np=10;


fig, ax = pyplot.subplots(figsize=(15, 10.0)) ;
#pyplot.figure(figsize=(15, 10.0));
pyplot.plot(xn_exact,un_exact,'-k',label='exact sol'); 

ylim_0 = list();
ylim_1 = list();

ylim_0.append(min(un_exact));
ylim_1.append(max(un_exact));

for ii in range(0,size(Beta)):

    fname_un_disc = dir1+discont+str("_N")+Nelem\
    +str("_CFL")+str(CFL)+str("_Beta")+str(Beta[ii])\
    +str("_")+str(T)+str("T.dat");

    data_disc= numpy.loadtxt(fname_un_disc);
    x_disc = data_disc[:,0];
    u_disc = data_disc[:,1];

    ylim_0.append(min(u_disc));
    ylim_1.append(max(u_disc));

    nn = size(x_disc);

    for i in range(0,size(x_disc)-1,Np):
        xx = x_disc[i:i+Np];
        uu = u_disc[i:i+Np];
        colorVal = scalarMap.to_rgba(float(Beta[ii]));
        pyplot.plot(xx,uu,color=colorVal); 

    xx = x_disc[i:i+Np];
    uu = u_disc[i:i+Np];
    #ll = str("Numerical sol, (")+r'$\beta= $'+str(Beta[ii])+" )";
    ll = r'$\beta= $'+str(Beta[ii]);  
    pyplot.plot(xx,uu,color=colorVal,label=ll); 
    #pyplot.plot(xx,uu,marker=cluster[ii],color=colorVal,label=ll,markevery=5);
    #pyplot.plot(xx,uu,ls=lines_[ii],color=colorVal,label=ll); 


pyplot.legend(loc='upper left');
CFL=Decimal(CFL.quantize(Decimal('.1')));
title_a = str("DGp")+ DG + " RK"+ RK \
+" for CFL="+ str(CFL)+ " and at t/T="+str(T);
#pyplot.title(title_a);

pyplot.xlabel('X');
pyplot.ylabel('u');

pyplot.xlim(-1.0,1.0);
pyplot.ylim(min(ylim_0)*1.05,max(ylim_1)*1.05);

from matplotlib.ticker import FormatStrFormatter
ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

xlabels=[-1.0,-0.5,0,0.5,1.0];
xlocs=[-1.0,-0.5,0,0.5,1.0];
pyplot.xticks(xlocs, xlabels);

ylabels=arange(0,1.2,0.2);
ylocs=arange(0,1.2,0.2);
pyplot.yticks(ylocs,ylabels);

pyplot.savefig('fig1.png',bbox='tight')

pyplot.show()












