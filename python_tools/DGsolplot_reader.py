
import argparse
from decimal import Decimal
import csv

parser = argparse.ArgumentParser(description='python_DG_argument_parsing');

parser.add_argument('-f', type=str, dest='python_input');

args = parser.parse_args();

with open(args.python_input) as file: 
    
    reader=csv.reader(file, delimiter=':');

    for row in reader:
        if row[0]=='mode':
            mode=str(row[1]);
        elif row[0]=='Eqn_set':
            eqn_set=str(row[1]);
        elif row[0]=='Diffusion_scheme':
            diffus_scheme=str(row[1]);
        elif row[0]=='CFL':
            #CFL=Decimal(row[1]);
            CFL = str(row[1]);
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
        elif row[0] == 'Beta':
            Beta=Decimal(row[1]);
        elif row[0] == 'Epsilon':
            Epsilon=Decimal(row[1]);
        elif row[0] == 'dt':
            dt_=str(row[1]);

from DGsolplot import plot_diffus, plot_advect, plot_AdvecDiffus

if eqn_set=='Advection':
    a =plot_advect(mode, DG, RK, CFL, Nelem, T, dt_\
                   , Beta, dir1, aver, nodal_exact, nodal_comp, discont )
elif eqn_set=='Diffusion':
    a = plot_diffus(diffus_scheme, mode, DG, RK, CFL, Nelem, T, dt_\
                  , Epsilon, dir1, aver, nodal_exact, nodal_comp, discont )
elif eqn_set=='Advection_Diffusion':
    plot_AdvecDiffus(diffus_scheme, mode, DG, RK, CFL, Nelem, T, dt_ \
                     , Beta, Epsilon, dir1, aver, nodal_exact, nodal_comp, discont)







