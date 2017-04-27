import numpy                 #loading our favorite library
from matplotlib import pyplot    #and the useful plotting library
# %matplotlib inline
#from numpy import *
from numpy import sin,cos,pi,linspace,ones,zeros,abs,min,max,exp, shape, empty_like 
import argparse
from decimal import Decimal
import os.path
import csv

parser = argparse.ArgumentParser(description='Short sample app');

parser.add_argument('-i', metavar='N', type=float, nargs=3, dest='data_args',
                    help='an integer for the accumulator');

parser.add_argument('-f', type=str, dest='python_input');

#parser.add_argument('-f', type=argparse.FileType('r'),dest='inputfile');
#if args.inputfile is not None:	
#	print args.inputfile.readlines();

args = parser.parse_args();

with open(args.python_input) as file: 
	reader=csv.reader(file, delimiter=':');
    	for row in reader:
		if row[0]=='CFL':
			CFL=row[1];
		elif row[0] == 'DGp': DG=row[1];
		elif row[0] == 'RK': RK=row[1];
   		elif row[0] == 'Beta': Beta=row[1];
		elif row[0] == 'T': T=row[1];

if args.data_args is not None:
	CFL = Decimal(args.data_args[0]); 
	CFL=Decimal(CFL.quantize(Decimal('.01')));
	Beta = Decimal(args.data_args[1]);
	Beta=Decimal(Beta.quantize(Decimal('.01')));
	T = args.data_args[2];
	T = int(T);

DG=int(DG);
RK=int(RK);
Beta=Decimal(Beta);
Beta=Decimal(Beta.quantize(Decimal('.01')));
CFL = Decimal(CFL);
CFL=Decimal(CFL.quantize(Decimal('.01')));
T=int(T);

fname_u_aver = str("u_aver_CFL")+str(CFL)+str("_Beta")+str(Beta)+str("_")+str(T)+str("T.dat");
data_aver= numpy.loadtxt(fname_u_aver);

fname_un_ex = str("u_nodal_exact_CFL")+str(CFL)+str("_Beta")+str(Beta)+str("_")+str(T)+str("T.dat");
data_exact_nodal= numpy.loadtxt(fname_un_ex);

fname_un_comp = str("u_nodal_CFL")+str(CFL)+str("_Beta")+str(Beta)+str("_")+str(T)+str("T.dat");
data_num_nodal= numpy.loadtxt(fname_un_comp);


xc = data_aver[:,0];
u_aver_comp = data_aver[:,1];
u_aver_exact = data_aver[:,2];

xn_exact = data_exact_nodal[:,0]; 
un_exact= data_exact_nodal[:,1]; 

xn_comp = data_num_nodal[:,0];
un_comp = data_num_nodal[:,1];


pyplot.figure();

pyplot.plot(xn_comp,un_comp,'-or',label='Numerical sol');
pyplot.plot(xn_exact,un_exact,'--k',label='Exact sol');
pyplot.grid();
pyplot.legend();
pyplot.xlabel('x');
pyplot.ylabel('u');
pyplot.title('Nodal Solution Comparison');

pyplot.figure();

pyplot.plot(xc,u_aver_comp,'or',label='Numerical average sol');
pyplot.plot(xc,u_aver_exact,'sk',label='Exact average sol');
pyplot.grid();
pyplot.legend();
pyplot.xlabel('xc');
pyplot.ylabel('u_average');
pyplot.title('Average Solution Comparison');

pyplot.figure();

pyplot.plot(xc,u_aver_comp,'-or',label='Numerical average sol');
pyplot.grid();
pyplot.legend();
pyplot.xlabel('xc');
pyplot.ylabel('u_average');


pyplot.show()










