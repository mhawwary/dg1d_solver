
import argparse
from decimal import Decimal
import csv

from numpy import pi

from subprocess import call, run, PIPE, Popen
from sys_cmd_toolbox import system_process    # locally defined file

parser = argparse.ArgumentParser(description='python_DG_argument_parsing');

parser.add_argument('-i', type=str, dest='python_input');
parser.add_argument('-t', type=float, dest='plot_time');
parser.add_argument('-b', type=int, dest='burger_plot_flag');
parser.add_argument('-dEdt', type=int, dest='compute_dEdt_flag');
parser.add_argument('-K_den', type=float, dest='K_den');  # K = pi/k_den

args = parser.parse_args();

with open(args.python_input) as file: 
    
    reader=csv.reader(file, delimiter=':');

    for row in reader:
        if row[0]=='mode':
            mode=str(row[1]);
        elif row[0]=='Eqn_set':
            eqn_set=str(row[1]);
        elif row[0]=='Eqn_type':
            eqn_type=str(row[1]);
        elif row[0]=='wave_form':
            wave_form_=str(row[1]);
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
        elif row[0] == 'thermal_diffus':
            gamma_=str(row[1]);
        elif row[0] == 'dt':
            dt_=str(row[1]);
            
    if eqn_set == 'Advection':
        Epsilon = None;
    elif eqn_set == 'Diffusion':
        Beta = None;

    if not(Epsilon is None):
        Epsilon = Decimal(Epsilon.quantize(Decimal('.01')))
    if not(Beta is None):    
        Beta = Decimal(Beta.quantize(Decimal('.01')))
    
    tt_   = Decimal(args.plot_time) 
    if eqn_set=='Diffusion':
        tt_ =  Decimal(tt_.quantize(Decimal('.0001')))
    else:
        tt_ =  Decimal(tt_.quantize(Decimal('.001')))
    T =  Decimal(T.quantize(Decimal('.001')))
    CFL = Decimal(CFL.quantize(Decimal('.0001')))
    
    cmd=['mkdir',dir_input+'tempfig/']
    cmd_out,cmd_err=system_process(cmd,1000)
    cmd=['mkdir',dir_input+'tempfig/eps/']
    cmd_out,cmd_err=system_process(cmd,1000)
    cmd=['mkdir',dir_input+'tempfig/png/']
    cmd_out,cmd_err=system_process(cmd,1000)
    cmd=['mkdir',dir_input+'tempdata/']
    cmd_out,cmd_err=system_process(cmd,1000)

K_den = args.K_den    

from DGsolplot import plot_diffus, plot_advec, plot_AdvecDiffus\
                     , plot_burgers_decay_turb, compare_diffus_schemes\
                     , plot_dissipation_rates, plot_initial_proj_diffus

if eqn_set=='Advection':
    a =plot_advec(dir_input, mode, DG, RK, CFL, Nelem, N_disc_ppt, tt_, dt_ \
                     , Beta, Epsilon, cont_num_time, disc_num_time,T, K_den)
        
elif eqn_set=='Diffusion':
    if tt_>0.0000:
        a = plot_diffus(dir_input, mode, DG, RK, CFL, Nelem, N_disc_ppt, tt_, dt_ \
        , Epsilon, gamma_, diffus_scheme, cont_num_time\
        , disc_num_time,T, K_den, wave_form_)
    else:
        a = plot_initial_proj_diffus(dir_input, mode, DG, RK, CFL, Nelem, N_disc_ppt, tt_, dt_ \
                 , Epsilon, gamma_, diffus_scheme, cont_num_time\
                 , disc_num_time,T, K_den, wave_form_)
                     
    #a = compare_diffus_schemes(dir_input, mode, DG, RK, CFL, Nelem, N_disc_ppt, tt_, dt_ \
     #                , Epsilon, gamma_, diffus_scheme, cont_num_time, disc_num_time,T, K_den, wave_form_)

elif eqn_set=='Advection_Diffusion':
    if args.burger_plot_flag==1:
        if not(args.plot_time is None):
            plot_burgers_decay_turb(dir_input, mode, DG, RK, CFL, Nelem\
                         ,N_disc_ppt,tt_,dt_,Beta,Epsilon,gamma_,diffus_scheme\
                         ,cont_num_time, disc_num_time)
        else:
            print('bad option for time to plot')

        if args.compute_dEdt_flag==1:
            plot_dissipation_rates(dir_input, mode, DG, RK, CFL, Nelem\
                           ,N_disc_ppt,tt_,dt_,Beta,Epsilon,gamma_,diffus_scheme\
                           ,cont_num_time, disc_num_time)
    else:
        plot_AdvecDiffus(diffus_scheme, mode, DG, RK, CFL, Nelem, T, dt_ \
                        , Beta, Epsilon, gamma_, dir1, aver, nodal_exact, \
                        nodal_comp, discont)







