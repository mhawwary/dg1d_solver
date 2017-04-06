import numpy                 #loading our favorite library
from matplotlib import pyplot    #and the useful plotting library
# %matplotlib inline
#from numpy import *
from numpy import sin,cos,pi,linspace,ones,zeros,abs,min,max,exp, shape, empty_like 

x = numpy.loadtxt('x.dat');
xc = numpy.loadtxt('Xc.dat');

u_comp = numpy.loadtxt('u_final.dat');
u_exact = numpy.loadtxt('u_exact.dat');

u_aver_exact = numpy.loadtxt('u_aver_exact.dat');
u_aver_comp = numpy.loadtxt('u_aver_final.dat');


pyplot.figure();

pyplot.plot(x,u_comp,'-or',label='Numerical sol');
pyplot.plot(x,u_exact,'--k',label='Exact sol');
pyplot.grid();
pyplot.legend();
pyplot.xlabel('x');
pyplot.ylabel('u');
pyplot.title('Nodal Solution Comparison');

pyplot.figure();

pyplot.plot(xc,u_aver_comp,'-or',label='Numerical average sol');
pyplot.plot(xc,u_aver_exact,'--k',label='Exact average sol');
pyplot.grid();
pyplot.legend();
pyplot.xlabel('xc');
pyplot.ylabel('u_average');
pyplot.title('Average Solution Comparison');

pyplot.show()






