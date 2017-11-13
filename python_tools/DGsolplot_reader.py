
import argparse
from decimal import Decimal
import csv

from numpy import pi

from subprocess import call, run, PIPE, Popen
from sys_cmd_toolbox import system_process    # locally defined file

parser = argparse.ArgumentParser(description='python_DG_argument_parsing');

parser.add_argument('-f', type=str, dest='python_input');
parser.add_argument('-t', type=float, dest='plot_time');
parser.add_argument('-o', type=int, dest='burger_plot_flag');
parser.add_argument('-a', type=int, dest='comp_fft_flag');

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
            CFL=Decimal(row[1]);
            #CFL = str(row[1]);
        elif row[0]=='p':
            DG=str(row[1]);
        elif row[0] == 'RK':    
            RK=str(row[1]);
        elif row[0] == 'Nelem':    
            Nelem=str(int(row[1]));
        elif row[0] == 'N_disc_ppt':    
            N_disc_ppt=int(row[1])
        elif row[0] == 'T':     
            T=Decimal(row[1]);
        elif row[0] =='dir':    
            dir_input = str(row[1]);
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
            
    tt_   = Decimal(args.plot_time)
    tt_ =  Decimal(tt_.quantize(Decimal('.001')))
    Beta = Decimal(Beta.quantize(Decimal('.01')))
    
    Epsilon = None;
    if not(Epsilon is None):
        Epsilon = Decimal(Epsilon.quantize(Decimal('.01')))
    CFL = Decimal(CFL.quantize(Decimal('.0001')))
    
    cmd=['mkdir',dir_input+'tempfig/']
    cmd_out,cmd_err=system_process(cmd,1000)

from DGsolplot import plot_diffus, plot_advec, plot_AdvecDiffus, plot_burgers_decay_turb

if eqn_set=='Advection':
    a =plot_advec(dir_input, mode, DG, RK, CFL, Nelem, N_disc_ppt, tt_, dt_ \
                     , Beta, Epsilon, cont_num_time, disc_num_time)
                   
    if args.comp_fft_flag==1:
        from fft_toolbox_python_new import load_data, compute_fft, plot_fft
        
        if (mode=='test') | (mode =='normal') | (mode =='CFL_const'):
            fname = dir_input + cont_sol + str("_N") + str(Nelem) \
                             + str("_CFL") + str(CFL) + str("_Beta") + str(Beta) \
                             + str("_") + str(tt_) + str("t.dat")  
        else:
            fname = dir_input + cont_sol + str("_N") + str(Nelem) \
                             + str("_dt") + dt_ + str("_Beta") + str(Beta) \
                             + str("_Eps") + str(Epsilon) \
                             + str("_") + str(tt_) + str("t.dat")
                             
        x_data_, u_data_ = load_data(fname)
        k_freq, u_amp, KEnerg = compute_fft(u_data_)
        
        plot_fft(k_freq*2*2*pi/80,u_amp)
        
elif eqn_set=='Diffusion':
    a = plot_diffus(diffus_scheme, mode, DG, RK, CFL, Nelem, T, dt_\
                  , Epsilon, dir1, aver, nodal_exact, nodal_comp, discont )

elif args.burger_plot_flag==1:
    if not(args.burgers_plot_time is None):
        plot_burgers_decay_turb(dir_input, mode, DG, RK, CFL, Nelem, tt_, dt_ \
                     , Beta, Epsilon, cont_num_time, disc_num_time)
    else:
        print('bad option for time to plot')

elif eqn_set=='Advection_Diffusion':
    plot_AdvecDiffus(diffus_scheme, mode, DG, RK, CFL, Nelem, T, dt_ \
                     , Beta, Epsilon, dir1, aver, nodal_exact, nodal_comp, discont)







