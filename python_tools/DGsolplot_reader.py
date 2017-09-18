
import argparse
from decimal import Decimal
import csv

parser = argparse.ArgumentParser(description='python_DG_argument_parsing');

parser.add_argument('-f', type=str, dest='python_input');
parser.add_argument('-t', type=float, dest='burgers_plot_time');
parser.add_argument('-o', type=int, dest='burger_plot_flag');

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
        elif row[0] == 'cont_unsteady_num': 
            cont_num_time = str(row[1]);
        elif row[0] == 'disc_unsteady_num': 
            disc_num_time = str(row[1]);
        elif row[0] == 'Beta':
            Beta=Decimal(row[1]);
        elif row[0] == 'Epsilon':
            Epsilon=Decimal(row[1]);
        elif row[0] == 'dt':
            dt_=str(row[1]);

from DGsolplot import plot_diffus, plot_advect, plot_AdvecDiffus, plot_burgers_decay_turb

if eqn_set=='Advection':
    a =plot_advect(mode, DG, RK, CFL, Nelem, T, dt_\
                   , Beta, dir1, aver, nodal_exact, nodal_comp, discont )
elif eqn_set=='Diffusion':
    a = plot_diffus(diffus_scheme, mode, DG, RK, CFL, Nelem, T, dt_\
                  , Epsilon, dir1, aver, nodal_exact, nodal_comp, discont )

elif args.burger_plot_flag==1:
    if not(args.burgers_plot_time is None):
        tt_=Decimal(args.burgers_plot_time)
        tt_ =  Decimal(tt_.quantize(Decimal('.001')));
        plot_burgers_decay_turb(dir1, DG, RK, CFL, Nelem, tt_, dt_ \
                     , Beta, Epsilon, cont_num_time, disc_num_time)
    else:
        print('bad option for time to plot')

elif eqn_set=='Advection_Diffusion':
    plot_AdvecDiffus(diffus_scheme, mode, DG, RK, CFL, Nelem, T, dt_ \
                     , Beta, Epsilon, dir1, aver, nodal_exact, nodal_comp, discont)







