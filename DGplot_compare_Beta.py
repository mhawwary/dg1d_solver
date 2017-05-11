import numpy                 #loading our favorite library
from matplotlib import pyplot    #and the useful plotting library
# %matplotlib inline
#from numpy import *
from numpy import sin,cos,pi,linspace,ones,zeros,abs,min,max,exp, shape, empty_like , size 
import argparse
from decimal import Decimal
import csv

pyplot.rc('legend',**{'loc':'upper left'});
pyplot.rcParams[u'legend.fontsize'] = 18
pyplot.rcParams[u'font.weight']='normal'
pyplot.rcParams[u'xtick.labelsize']=15
pyplot.rcParams[u'ytick.labelsize']=15
pyplot.rcParams[u'axes.titlesize']=22
pyplot.rcParams[u'axes.labelsize']=22
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
            T=int(row[1]);
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

fname_un_ex = dir1+nodal_exact+str("_N")+Nelem\
+str("_CFL")+str(CFL)+str("_Beta")\
+str(Beta[0])+str("_")+str(T)+str("T.dat");

data_exact_nodal= numpy.loadtxt(fname_un_ex);
xn_exact = data_exact_nodal[:,0]; 
un_exact= data_exact_nodal[:,1]; 

import matplotlib.colors as colors
import matplotlib.cm as cmx

jet = cm = pyplot.get_cmap('jet') 
cNorm  = colors.Normalize(vmin=0, vmax=float(Beta[-1]))
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet);

#! Plotting The figure:
#---------------------------
fig = pyplot.figure(figsize=(15, 10.0));

pyplot.plot(xn_exact,un_exact,'-k',label='Exact sol'); 

for ii in range(0,size(Beta)):
    fname_un_disc = dir1+discont+str("_N")+Nelem+str("_CFL")+str(CFL)+str("_Beta")+str(Beta[ii])+str("_")+str(T)+str("T.dat");
    data_disc= numpy.loadtxt(fname_un_disc);
    x_disc = data_disc[:,0];
    u_disc = data_disc[:,1];

    nn = size(x_disc);
    Np = 2;

    for i in range(0,size(x_disc)-1,Np):
        xx = x_disc[i:i+Np];
        uu = u_disc[i:i+Np];
        colorVal = scalarMap.to_rgba(float(Beta[ii]));
        pyplot.plot(xx,uu,color=colorVal); 

    xx = x_disc[i:i+Np];
    uu = u_disc[i:i+Np];
    ll = str("Numerical sol, (")+r'$\beta= $'+str(Beta[ii])+" )"; 
    pyplot.plot(xx,uu,color=colorVal,label=ll); 


pyplot.legend(loc='upper left');
title_a = str("DGp")+ DG + " RK"+ RK +" for CFL="+ str(CFL)+ " and at t/T="+str(T);
pyplot.title(title_a);
pyplot.xlabel('X');
pyplot.ylabel('u(x)');

pyplot.savefig('fig1.png',bbox='tight')

pyplot.show()












