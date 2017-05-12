from matplotlib import pyplot ,ticker   #and the useful plotting library
# %matplotlib inline
#from numpy import *
from numpy import sin,cos,pi,linspace,ones,zeros,abs,min,max,exp, shape, empty_like , size, loadtxt, arange , log
import argparse
from decimal import Decimal
import csv

parser = argparse.ArgumentParser(description='python_DG_argument_parsing');

parser.add_argument('-f', type=str, dest='python_input');

args = parser.parse_args();

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
        elif row[0] == 'Beta':
            Beta=Decimal(row[1]); 

Beta= Decimal(Beta.quantize(Decimal('.01')));
CFL=Decimal(CFL.quantize(Decimal('.01')));
T=Decimal(T.quantize(Decimal('.1')));

fname_u_aver = dir1+aver+str("_N")+Nelem+str("_CFL")+str(CFL)+str("_Beta")+str(Beta)+str("_")+str(T)+str("T.dat");
data_aver= loadtxt(fname_u_aver);

fname_un_ex = dir1+nodal_exact+str("_N")+Nelem+str("_CFL")+str(CFL)+str("_Beta")+str(Beta)+str("_")+str(T)+str("T.dat");
data_exact_nodal= loadtxt(fname_un_ex);

fname_un_comp = dir1+nodal_comp+str("_N")+Nelem+str("_CFL")+str(CFL)+str("_Beta")+str(Beta)+str("_")+str(T)+str("T.dat");
data_num_nodal= loadtxt(fname_un_comp);

fname_un_disc = dir1+discont+str("_N")+Nelem+str("_CFL")+str(CFL)+str("_Beta")+str(Beta)+str("_")+str(T)+str("T.dat");
data_disc= loadtxt(fname_un_disc);

xc = data_aver[:,0];
u_aver_comp = data_aver[:,1];
u_aver_exact = data_aver[:,2];

xn_exact = data_exact_nodal[:,0]; 
un_exact= data_exact_nodal[:,1]; 

xn_comp = data_num_nodal[:,0];
un_comp = data_num_nodal[:,1];

x_disc = data_disc[:,0];
u_disc = data_disc[:,1];

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

nn = size(x_disc);

if int(DG)<=1:
    Np = 2;
else:
    Np=10;

pyplot.figure();

ylim_0 = list();
ylim_1 = list();
ylim_0.append(min(un_exact));
ylim_1.append(max(un_exact));

pyplot.plot(xn_exact,un_exact,'--k',label='Exact sol'); 

for i in range(0,size(x_disc)-1,Np):
    xx = x_disc[i:i+Np];
    uu = u_disc[i:i+Np];
    pyplot.plot(xx,uu,'-r'); 

ylim_0.append(min(u_disc));
ylim_1.append(max(u_disc));

xx = x_disc[i:i+Np];
uu = u_disc[i:i+Np];
ll = str("Numerical sol, (")+r'$\beta= $'+str(Beta)+" )"; 
pyplot.plot(xx,uu,'-r',label=ll); 

pyplot.legend();

title_a = str("DGp")+ DG + " RK"+ RK \
+ ", and upwind_param ("+r'$\beta= $'\
+str(Beta)+")\n for CFL="+ str(CFL)\
+ " and at t/T="+str(T);

pyplot.title(title_a);
pyplot.xlabel('X');
pyplot.ylabel('u(x)');

pyplot.xlim(-1.0,1.0);
pyplot.ylim(min(ylim_0)*1.05,max(ylim_1)*1.05);

xlabels=[-1.0,-0.5,0,0.5,1.0];
xlocs=[-1.0,-0.5,0,0.5,1.0];
pyplot.xticks(xlocs, xlabels);

ylabels=arange(0,1.2,0.2);
ylocs=arange(0,1.2,0.2);
pyplot.yticks(ylocs,ylabels);

pyplot.show()










