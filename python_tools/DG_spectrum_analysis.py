import argparse

parser = argparse.ArgumentParser(description='python_DG_argument_parsing');

parser.add_argument('-f', type=str, dest='dg_python_input')
parser.add_argument('-o', type=str, dest='type_of_analysis')  # dt_conv / mesh_conv 

args = parser.parse_args();

if not(args.dg_python_input is None ):
    from fft_toolbox_python_new import DG_fft_input_reader
    dir_result_dg, dir_input_dg, mode, poly_order, RK, Beta, Epsilon, Nelem_dg, CFL\
    , dt, T, tt_, N_aver_fft, dg_cont_sol, dg_disc_sol = DG_fft_input_reader(args.dg_python_input)
    
    if not(args.type_of_analysis is None ):
        if args.type_of_analysis== 'mesh_conv':    # mesh_conv analysis
            # loading DG fft_toolbox_functions:
            from fft_toolbox_python_new import compute_fft
            Nelem = [40] # p3
            #Nelem_mesh_refine = [51,102,204,409]  # p5
            compute_cont_fft(dir_input_dg, dir_result_dg, mode, poly_order\
            , RK, dg_cont_sol, Nelem_dg_mesh_refine, Beta, Epsilon, CFL, dt, tt_, N_aver_fft)
            
        elif args.type_of_analysis == 'dt_conv':    # dt_conv analysis
            # DG dt_ convergence analysis:
            from fft_toolbox_python_new import perform_dt_convergence_study_DG
            dt_ = [str('1.000e-05')]
            perform_dt_convergence_study_DG(dir_input, dir_dg, poly_order, RK, dg_sol, Nelem_dg, Beta, Epsilon, dt_, tt_, N_aver_fft)  
            
        elif args.type_of_analysis == 'compare_DG':  # compare 2 DG solutions
            dir_input2 = './input/....';
            from fft_toolbox_python_new import compare_fft_DG
            compare_fft_DG(dir_input, dir_input2, dir_dg, dg_sol, Nelem_dg, Beta, Epsilon, RK, dt_, tt_, N_aver_fft)
    else:
        print('error type of analysis is undefined') 

else:
    print('error input file is not specified')
    
