import numpy                 #loading our favorite library
from matplotlib import pyplot    #and the useful plotting library
# %matplotlib inline
#from numpy import *
from numpy import sin,cos,pi,linspace,ones,zeros,abs,min,max,exp, shape, empty_like 
import argparse

parser = argparse.ArgumentParser(description='Short sample app');

parser.add_argument('-i', metavar='N', type=float, nargs=3, dest='b',
                    help='an integer for the accumulator')

results = parser.parse_args();

from decimal import Decimal
CFL = Decimal(results.b[0]); 
CFL=Decimal(CFL.quantize(Decimal('.01')));
Beta = Decimal(results.b[1]);
Beta=Decimal(Beta.quantize(Decimal('.01')));
T = results.b[2];
T = int(T);

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










