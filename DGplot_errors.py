import numpy                 #loading our favorite library
from matplotlib import pyplot    #and the useful plotting library
from numpy import sin,cos,pi,linspace,ones,zeros,abs,min,max,exp, shape, empty_like , size , arange, log
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

#Beta = list();

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
        elif row[0] == 'errors': 
            errors = str(row[1]);

CFL=Decimal(CFL.quantize(Decimal('.01')));
T=Decimal(T.quantize(Decimal('.1')));

fname_errors = dir1+errors+str("_N")+Nelem\
+str("_CFL")+str(CFL)+"_allBeta_"\
+str(T)+str("T.dat");

data_errors= numpy.loadtxt(fname_errors);
Beta = data_errors[:,0]; 
L2_proj= data_errors[:,1]; 
L2_aver= data_errors[:,2];

#! Plotting The figure:
#---------------------------

fig, ax = pyplot.subplots(figsize=(15, 10.0)) ;

pyplot.plot(Beta,L2_proj,'-or',label='L2 of projected polynomials'); 
pyplot.yscale('log');

CFL=Decimal(CFL.quantize(Decimal('.01')));
title_a = str("DGp")+ DG + " RK"+ RK \
+" for CFL="+ str(CFL)+ " and at t/T="+str(T);
pyplot.title(title_a);
pyplot.xlabel(r'$\beta$');
pyplot.ylabel('L2 error');

pyplot.ylim(1.0e-2,6.5e-2)
from matplotlib import ticker 
#ax.yaxis.set_major_formatter(ticker.LogFormatter())
ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=15))



pyplot.show()


