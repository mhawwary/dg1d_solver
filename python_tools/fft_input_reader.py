
import argparse
from decimal import Decimal
import csv
from numpy import zeros, vstack, append
from subprocess import call, run, PIPE, Popen

parser = argparse.ArgumentParser(description='python_FD_argument_parsing');

parser.add_argument('-f', type=str, dest='fd_python_input')
parser.add_argument('-g', type=str, dest='dg_python_input')

args = parser.parse_args();

if not(args.dg_python_input is None ):
    from fft_toolbox_python import DG_fft_input_reader
    dir_result_dg_, dir_input_dg_, DG, RK, Beta, Epsilon, Nelem_dg, CFL_dg\
    , dt_dg, T, tt_, N_aver_fft, dg_sol = DG_input_reader(args.dg_python_input)
elif not(args.fd_python_input is None ):
    from fft_toolbox_python import FD_fft_input_reader
#==============================================================================================================#
#                                     Compare FD and DG 
#==============================================================================================================#
if args.type_of_analysis == 'compare_FD_DG_DNS':
    from fft_toolbox_python import compare_fft_FD_DG4, compare_fft_FD_DG3, DG_input_reader, FD_input_reader
    
    dir_result_fd_, dir_input_fd_,FDOA, RK, Nelem_fd, CFL_fd, dt_fd, T, tt_\
    , N_aver_fft, fd_sol = FD_input_reader(args.fd_python_input)
    #compare_fft_FD_DG(dir_input_fd_, dir_result_fd_, FDOA, fd_sol, Nelem_fd\
    #, dt_fd, dir_input_dg_, dir_result_dg_, DG, dg_sol, Nelem_dg, Beta, Epsilon\
    #, RK, dt_dg, tt_, N_aver_fft)
    
    #dir_input_dg_2 = args.dg_python_input_dir;
    #compare_fft_FD_DG2(dir_input_fd_, dir_result_fd_, FDOA, fd_sol, Nelem_fd\
    #, dt_fd, dir_input_dg_, dir_input_dg_2, dir_result_dg_, DG, dg_sol, Nelem_dg, Beta, Epsilon\
    #, RK, dt_dg, tt_, N_aver_fft)
    
    compare_fft_FD_DG3(dir_input_fd_, dir_result_fd_, FDOA, fd_sol, Nelem_fd, CFL_fd\
    , dir_input_dg_, dir_result_dg_, DG, dg_sol, Nelem_dg, Beta, Epsilon, RK, CFL_dg\
    , tt_, N_aver_fft)
    
    compare_fft_FD_DG4(dir_input_fd_, dir_result_fd_, FDOA, fd_sol, Nelem_fd, dt_fd\
    , dir_input_dg_, dir_result_dg_, DG, dg_sol, Nelem_dg, Beta, Epsilon, RK, dt_dg\
    , tt_, N_aver_fft)

#==============================================================================================================#
#                                             DG input Reader 
#==============================================================================================================# 
elif not(args.dg_python_input is None ):    # DG only anlaysis
    with open(args.dg_python_input) as file: 
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
                Nelem_dg=str(int(row[1]));
            elif row[0] == 'T':     
                T=Decimal(row[1]);
            elif row[0] =='dir_result':    
                dir_result = str(row[1]);
            elif row[0] =='dir_input':    
                dir_input = str(row[1]);
            elif row[0] =='aver':   
                aver = str(row[1]);
            elif row[0] == 'cont_exact':
                nodal_exact = str(row[1]);
            elif row[0] == 'cont_num':
                nodal_comp = str(row[1]);
            elif row[0] == 'cont_unsteady_num':
                dg_sol = str(row[1]);
            elif row[0] == 'discont': 
                discont = str(row[1]);
            elif row[0] == 'Beta':
                Beta=Decimal(row[1]);
            elif row[0] == 'Epsilon':
                Epsilon=Decimal(row[1]);
            elif row[0] == 'dt':
                dt_=str(row[1]);
            elif row[0] == 'time_for_fft_compute':
                tt_ = Decimal(row[1])
            elif row[0] == 'N_ensembles_fft':
                N_aver_fft = str(row[1])
            elif row[0] == 'thermal_diffus':
                thermal_diffus = float(row[1])
                
    call(['mkdir','Results'])
    dir_dg = dir_result + str('p') + DG + str('_RK') + RK + str('/')
    call(['mkdir',dir_dg])

    if args.type_of_analysis== 'mesh_conv':    # mesh_conv analysis
        # loading DG fft_toolbox_functions:
        from fft_toolbox_python import compute_fft_mesh_refine_DG
        Nelem_dg_mesh_refine = [64,128,256,512,1024,2048,4096,8192]  # p3
        #Nelem_mesh_refine = [51,102,204,409]  # p5
        compute_fft_mesh_refine_DG(dir_dg,DG, RK, dg_sol, Nelem_dg_mesh_refine, Beta, Epsilon, dt_, tt_, N_aver_fft)
        
    elif args.type_of_analysis == 'dt_conv':    # dt_conv analysis
        # DG dt_ convergence analysis:
        from fft_toolbox_python import perform_dt_convergence_study_DG
        dt_ = [str('3.000e-03'), str('5.000e-06')]
        perform_dt_convergence_study_DG(dir_input, dir_dg, DG, RK, dg_sol, Nelem_dg, Beta, Epsilon, dt_, tt_, N_aver_fft)  
        
    elif args.type_of_analysis == 'compare_DG':  # compare 2 DG solutions
        if not(args.dg_python_input_dir is None ):
            dir_input2 = args.dg_python_input_dir;
        from fft_toolbox_python import compare_fft_DG
        compare_fft_comapre_DG(dir_input, dir_input2, dir_dg, dg_sol, Nelem_dg, Beta, Epsilon, RK, dt_, tt_, N_aver_fft) 
    

#==============================================================================================================#
#                                             FD input Reader 
#==============================================================================================================#
elif not(args.fd_python_input is None ):        # FD only anlaysis    
    with open(args.fd_python_input) as file: 
        
        reader=csv.reader(file, delimiter=':');

        for row in reader:
            if row[0]=='mode':
                mode=str(row[1]);
            elif row[0]=='Eqn_set':
                eqn_set=str(row[1]);
            elif row[0]=='CFL':
                #CFL=Decimal(row[1]);
                CFL = str(row[1]);
            elif row[0]=='FDOA':
                FD=str(row[1]);
            elif row[0] == 'RK':    
                RK=str(row[1]);
            elif row[0] == 'Nelem':    
                Nelem_fd=str(int(row[1]));
            elif row[0] == 'Nexact':    
                Nexact=str(int(row[1]));
            elif row[0] == 'T':     
                T=Decimal(row[1]);
            elif row[0] =='dir_result':    
                dir_result = str(row[1]);
            elif row[0] =='dir_input':    
                dir_input = str(row[1]);
            elif row[0] == 'dt':
                dt_=str(row[1]);
            elif row[0] == 'time_for_fft_compute':
                tt_ = Decimal(row[1])
            elif row[0] == 'N_ensembles_fft':
                N_aver_fft = str(row[1])
            elif row[0] == 'exact': 
                sol_exact = str(row[1]);
            elif row[0] == 'numerical': 
                fd_sol = str(row[1]);

    if FD==2:
        FDOA = FD + str('nd')
    elif FD==3:
        FDOA = FD + str('rd')
    else:
        FDOA = FD + str('th')
     
    call(['mkdir','Results'])
    dir_fd = dir_result + str('FD') + FD +str('th')+ str('_RK') + RK + str('/')
    call(['mkdir',dir_fd])

    if args.type_of_analysis== 'mesh_conv':  # mesh_conv analysis
        # FD Mesh Refinement Analysis:
        from fft_toolbox_python import compute_fft_mesh_refine_FD
        Nelem_fd_mesh_refine = [255,511,1023,2047,4095,8191,16383,32767]  # OA4
        compute_fft_mesh_refine_FD(dir_fd,FD, FDOA, RK, fd_sol, Nelem_fd_mesh_refine, dt_, tt_, N_aver_fft) 
    elif args.type_of_analysis == 'dt_conv': 
        # FD dt_ convergence analysis:       # dt_conv analysis
        from fft_toolbox_python import perform_dt_convergence_study_FD
        #dt_ = [str('1.000e-06'), str('2.500e-06'),str('5.000e-06'),str('1.000e-05'),str('5.000e-05')]
        dt_ = [str('1.000e-02'), str('5.000e-03'),str('5.000e-06')]
        perform_dt_convergence_study_FD(dir_fd, FDOA, RK, fd_sol, Nelem_fd, dt_, tt_, N_aver_fft)
        
    






