import numpy                 #loading our favorite library
from matplotlib import pyplot    #and the useful plotting library
from numpy import sin,cos,pi,linspace,ones,zeros,abs,min,max\
,exp, shape, empty_like , size , arange, log, log10, loadtxt, savetxt
import argparse
from decimal import Decimal
import csv

pyplot.rc('legend',**{'loc':'upper left'});
pyplot.rcParams[u'legend.fontsize'] = 18
pyplot.rcParams[u'legend.edgecolor']='white'
pyplot.rcParams[u'legend.facecolor']='0.8'
pyplot.rcParams[u'font.weight']='normal'
pyplot.rcParams[u'xtick.labelsize']=16
pyplot.rcParams[u'ytick.labelsize']=16
pyplot.rcParams[u'axes.titlesize']=18
pyplot.rcParams[u'axes.labelsize']=20
pyplot.rcParams[u'axes.spines.right']='false';
pyplot.rcParams[u'axes.spines.top']='false';
pyplot.rcParams[u'lines.linewidth'] = 1.5;
pyplot.rcParams[u'lines.markersize'] = 8;

#===================================
# Loading Data 
#===================================
parser = argparse.ArgumentParser(description='python_DG_argument_parsing');

parser.add_argument('-f', type=str, dest='python_input');

args = parser.parse_args();

with open(args.python_input) as file: 
    
    reader=csv.reader(file, delimiter=':');

    for row in reader:
        if row[0]=='mode':
            mode=str(row[1]);
        elif row[0] == 'wave':
            wave = str(row[1]);
        elif row[0]=='CFL':    
            CFL=Decimal(row[1]);
        elif row[0] == 'CFL_upw':
            CFL_upw = Decimal(row[1]);
        elif row[0] == 'dt':
            dt_=str(row[1]);
        elif row[0]=='DGp':
            DG=str(row[1]);
        elif row[0] == 'RK':    
            RK=str(row[1]);
        elif row[0] == 'Nelem':    
            Nelem_=str(int(row[1]));
        elif row[0] == 'T':     
            T=Decimal(row[1]);
        elif row[0] =='dir':    
            dir1 = str(row[1]);
        elif row[0] == 'errors': 
            errors = str(row[1]);
        elif row[0]=='Beta':
            Beta = Decimal(row[1]);

CFL = Decimal(CFL.quantize(Decimal('.001')));
CFL_upw = Decimal(CFL_upw.quantize(Decimal('.001')));
T=Decimal(T.quantize(Decimal('.001')));
Beta=Decimal(Beta.quantize(Decimal('.01')));
dt = float(dt_);

if mode=='dt_const':
    fname_errors = dir1+errors+str("_dt")+dt_+"_Beta"+str(Beta)\
                                  +'_'+str(T)+str("T.dat");
    fname_upw = dir1+errors+str("_dt")+dt_+"_Beta1.00"\
                               +'_'+str(T)+str("T.dat");
elif mode=='CFL_const':
    fname_errors = dir1 + errors + str("_CFL") + str(CFL) \
                   + "_Beta" + str(Beta) + '_' + str(T) + str("T.dat");
    fname_upw = dir1 + errors + str("_CFL") + str(CFL_upw) \
                   + "_Beta1.00" + '_' + str(T) + str("T.dat");

data_errors = loadtxt(fname_errors);
data_upw = loadtxt(fname_upw);

Nelem = data_errors[:,0]; 
nDOF = Nelem * (int(DG)+1);

L1_proj= data_errors[:,1];
L1_aver= data_errors[:,2];
L2_proj= data_errors[:,3];
L2_aver= data_errors[:,4];

L1_proj_upw = data_upw[:, 1];
L1_aver_upw = data_upw[:, 2];
L2_proj_upw = data_upw[:, 3];
L2_aver_upw = data_upw[:, 4];

# theoritcal curve:
theoretical_curve = exp(-(int(DG)+1) * log(nDOF))  ;
shift =  -0.7*log10(L2_proj[0]) + log10(theoretical_curve[0]);
theoretical_curve = theoretical_curve / 10**shift;

#===================================
# Order Calculations:
#===================================
order_L1_proj = zeros(size(L1_proj)-1);
order_L1_aver = zeros(size(L1_aver)-1);
order_L2_proj = zeros(size(L2_proj)-1);
order_L2_aver = zeros(size(L2_aver)-1);
order_exact = zeros(size(theoretical_curve)-1);

for i in range(1,size(L2_proj)):
    order_L1_proj[i-1] = log10(L1_proj[i-1]/L1_proj[i])/log10(2);
    order_L1_aver[i-1] = log10(L1_aver[i-1]/L1_aver[i])/log10(2);
    order_L2_proj[i-1] = log10(L2_proj[i-1]/L2_proj[i])/log10(2);
    order_L2_aver[i-1] = log10(L2_aver[i-1]/L2_aver[i])/log10(2);
    order_exact[i-1] = \
               log10(theoretical_curve[i-1]/theoretical_curve[i])/log10(2);

Nelem_ = zeros(size(Nelem)-1);
for i in range(1,size(Nelem)):
	Nelem_[i-1] = int(Nelem[i]);

order_print = numpy.transpose([Nelem_,order_L1_proj,order_L1_aver\
                               ,order_L2_proj,order_L2_aver]);

order_L1_print = numpy.transpose([Nelem_,order_L1_proj,order_L1_aver]);
order_L2_print = numpy.transpose([Nelem_,order_L2_proj,order_L2_aver]);

print('L1_order',order_L1_print);
print('L2_order',order_L2_print);

if mode=='CFL_const':
    order_out_name = dir1+'OAA_DGp'+DG+'_RK'+RK+'_Beta'\
    +str(Beta)+'_CFL'+str(CFL)+'.dat';

elif mode=='dt_const':
    order_out_name = dir1+'OAA_DGp'+DG+'_RK'+RK+'_Beta'\
    +str(Beta)+'_dt'+dt_+'.dat';

savetxt(order_out_name, order_print\
, fmt="%02d"+" %1.4f"+" %1.4f"+" %1.4f"+" %1.4f"\
,header="Nelem, order_L1_proj, order_L1_aver, order_L2_proj, order_L2_aver"\
,delimiter=' ');


#===================================
#! Plotting The figure:
#===================================


if mode=='dt_const':
    title_a = str("DGp") + DG + " RK" + RK \
              + " with dt=" + '{:1.2e}'.format(dt) \
              + ' at t/T=' + str(Decimal(T.quantize(Decimal('.1'))))

elif mode=='CFL_const':

	title_a = str("DGp")+ DG + " RK"+ RK\
              + ' at t/T=' + str(Decimal(T.quantize(Decimal('.1'))));

fig, ax = pyplot.subplots(figsize=(9.0, 7.5));

pyplot.plot(nDOF, L2_proj, '-or' \
            , label=r'hybrid, $\beta$=' + str(Beta)  + ', CFL='+str(CFL));
#pyplot.plot(nDOF, L2_aver, '-^b' \
#            , label=r'L$_{2}$ of averages, $\beta$=' + str(Beta));
pyplot.plot(nDOF, L2_proj_upw \
            , '--sb', label=r'upwind, $\beta=$1.00' + ', CFL='+str(CFL_upw));
#pyplot.plot(nDOF, L2_aver_upw, '--vb', label=r'L$_{2}$ of averages, upwind');
pyplot.plot(nDOF, theoretical_curve, ':k', linewidth=1.5,label='p+1 slope');

pyplot.legend(loc='upper right');
pyplot.yscale('log');
pyplot.xscale('log');

#pyplot.title(title_a)
pyplot.xlabel('nDOFs')
pyplot.ylabel(r'L$_{2}$ error of solution')

# pyplot.ylim(1.0e-2,6.5e-2)
from matplotlib import ticker

# ax.yaxis.set_major_formatter(ticker.LogFormatter())
ax.xaxis.set_major_locator(ticker.LogLocator(base=10.0, numticks=15))

ax.grid(which='minor', linestyle=':');

pyplot.show()


