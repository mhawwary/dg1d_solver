import numpy                 #loading our favorite library
from matplotlib import pyplot    #and the useful plotting library
# %matplotlib inline
#from numpy import *
from numpy import sin,cos,pi,linspace,ones,zeros,abs,min,max,exp, shape, empty_like , size 
import argparse
from decimal import Decimal
import csv

parser = argparse.ArgumentParser(description='python_DG_argument_parsing');

parser.add_argument('-f', type=str, dest='python_input');

args = parser.parse_args();

with open(args.python_input) as file: 
	reader=csv.reader(file, delimiter=':');
    	for row in reader:
		if row[0]=='CFL':       CFL=Decimal(row[1]);
		elif row[0] == 'DGp':   DG=str(row[1]);
		elif row[0] == 'RK':    RK=str(row[1]);
		elif row[0] == 'Nelem':  Nelem=str(int(row[1]));
   		elif row[0] == 'Beta':  Beta=Decimal(row[1]);
		elif row[0] == 'T':     T=int(row[1]);
		elif row[0] =='dir':    dir1 = str(row[1]);
		elif row[0] =='aver':   aver = str(row[1]);
		elif row[0] == 'nodal_exact': nodal_exact = str(row[1]);
		elif row[0] == 'nodal_num': nodal_comp = str(row[1]);
		elif row[0] == 'discont': discont = str(row[1]);

Beta= Decimal(Beta.quantize(Decimal('.01')));
CFL=Decimal(CFL.quantize(Decimal('.01')));

fname_u_aver = dir1+aver+str("_N")+Nelem+str("_CFL")+str(CFL)+str("_Beta")+str(Beta)+str("_")+str(T)+str("T.dat");
data_aver= numpy.loadtxt(fname_u_aver);

fname_un_ex = dir1+nodal_exact+str("_N")+Nelem+str("_CFL")+str(CFL)+str("_Beta")+str(Beta)+str("_")+str(T)+str("T.dat");
data_exact_nodal= numpy.loadtxt(fname_un_ex);

fname_un_comp = dir1+nodal_comp+str("_N")+Nelem+str("_CFL")+str(CFL)+str("_Beta")+str(Beta)+str("_")+str(T)+str("T.dat");
data_num_nodal= numpy.loadtxt(fname_un_comp);

fname_un_disc = dir1+discont+str("_N")+Nelem+str("_CFL")+str(CFL)+str("_Beta")+str(Beta)+str("_")+str(T)+str("T.dat");
data_disc= numpy.loadtxt(fname_un_disc);

fname_un_disc_upw = dir1+discont+str("_N")+Nelem+str("_CFL")+str(CFL)+str("_Beta")+"1.00"+str("_")+str(T)+str("T.dat");
data_disc_upw= numpy.loadtxt(fname_un_disc_upw);


xc = data_aver[:,0];
u_aver_comp = data_aver[:,1];
u_aver_exact = data_aver[:,2];

xn_exact = data_exact_nodal[:,0]; 
un_exact= data_exact_nodal[:,1]; 

xn_comp = data_num_nodal[:,0];
un_comp = data_num_nodal[:,1];

x_disc = data_disc[:,0];
u_disc = data_disc[:,1];

x_disc_upw = data_disc_upw[:,0];
u_disc_upw = data_disc_upw[:,1];

#print pyplot.rcParams.keys() 
pyplot.rc('legend',**{'loc':'upper left'});
#pyplot.rcParams['font.size'] = 15
pyplot.rcParams[u'legend.fontsize'] = 18
#pyplot.rcParams['figure.titlesize'] = 20
pyplot.rcParams[u'font.weight']='normal'
pyplot.rcParams[u'xtick.labelsize']=15
pyplot.rcParams[u'ytick.labelsize']=15
pyplot.rcParams[u'axes.titlesize']=22
pyplot.rcParams[u'axes.labelsize']=22
pyplot.rcParams[u'lines.linewidth'] = 1.5;
pyplot.rcParams[u'lines.markersize'] = 8;
#pyplot.rcParams[u'font.family'] = 'SanSerif';

nn = size(x_disc);

Np = 5;

pyplot.figure();

pyplot.plot(xn_exact,un_exact,'--k',label='Exact sol'); 
pyplot.hold(True);
for i in range(0,size(x_disc)-2,Np):
	xx = x_disc[i:i+Np];
 	uu = u_disc[i:i+Np];
	pyplot.plot(xx,uu,'-r'); 
	pyplot.hold(True);

xx = x_disc[i:i+Np];
uu = u_disc[i:i+Np];
ll = str("Numerical sol, (")+r'$\beta= $'+str(Beta)+" )"; 
pyplot.plot(xx,uu,'-r',label=ll); 

for i in range(0,size(x_disc_upw)-2,Np):
	xx = x_disc_upw[i:i+Np];
 	uu = u_disc_upw[i:i+Np];
	pyplot.plot(xx,uu,'--b'); 
	pyplot.hold(True);

xx = x_disc_upw[i:i+Np];
uu = u_disc_upw[i:i+Np];
pyplot.plot(xx,uu,'--b',label='Numerical sol, upwind'); 

pyplot.legend();
title_a = str("DGp")+ DG + " RK"+ RK + ", and upwind_param ("+r'$\beta= $'+str(Beta)+")\n for CFL="+ str(CFL)+ " and at t/T="+str(T);
pyplot.title(title_a);
pyplot.xlabel('X');
pyplot.ylabel('u(x)');

#pyplot.figure();

#pyplot.plot(xc,u_aver_comp,'or',label='Numerical average sol');
#pyplot.plot(xc,u_aver_exact,'sk',label='Exact average sol');
#pyplot.grid();
#pyplot.legend();
#pyplot.xlabel(r'$X_{center}$');
#pyplot.ylabel(r'$\bar{u}$',fontsize=30);
#title_a = str("DGp")+ DG + " RK"+ RK + ", and upwind_param ("+r'$\beta= $'+str(Beta)+")\n for CFL="+ str(CFL)+ " and at t/T="+str(T);
#pyplot.title(title_a);
#pyplot.xticks();
#pyplot.yticks();
#pyplot.style.use('presentation');
#print(pyplot.style.available)

#import matplotlib as matlib
#print matlib.get_configdir()


pyplot.show()










