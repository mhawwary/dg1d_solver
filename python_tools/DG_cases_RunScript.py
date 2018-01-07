import argparse
from decimal import Decimal
import csv
from time import sleep
#from subprocess import call, run, PIPE, Popen
import fileinput
from numpy import size, pi

from sys_cmd_toolbox import system_process

parser = argparse.ArgumentParser(description='DG_argparsing')
parser.add_argument('-f', type=str, dest='inputfname')
parser.add_argument('-s', type=str, dest='spectrum_dir')
args = parser.parse_args()

fname = args.inputfname
spectrum_dir = args.spectrum_dir;

elem_num = 200  
Um = [0.0,75.0,100.0]
um_name = [200]
CFL_max = 0.073;
dx = 2*pi /elem_num;

# Computing and Running the first no. of elements :
#--------------------------------------------------

for lines in fileinput.input(fname, inplace=True):
    a = lines.split('=')
    b = a[0].strip()
    if b == 'spectrum_restart_flag':
       print('   {} = {}'.format(a[0].strip(),1))
    elif b == 'num_elements':
        print(' {} = {}'.format(a[0].strip(),elem_num))
    elif b == 'title':
        sim_title = a[1].strip()
        print(' {} = {}'.format(b,a[1].strip()))
    elif b == 'polynomial_order':
        P = a[1].strip()
        print('  {} = {}'.format(b,a[1].strip()))
    elif b == 'RK_order':
        Prk = a[1].strip()
        print('   {} = {}'.format(b,a[1].strip()))
    elif b == 'N_uniform_pts_per_elem':
        n_pts = a[1].strip()
        print('   {} = {}'.format(b,a[1].strip()))
    else:
       print(lines.rstrip())
    
res_dir = './Results_AdvecDiffus/'+ str(sim_title)+'/' 
log_name_dir = res_dir + 'p' + str(P) + 'RK' + str(Prk) +'/case'
solve_dir =  res_dir + 'p' + str(P) + 'RK' + str(Prk) +'/'

print('input fname:  ',fname,' .............')
print('title: ',sim_title,'    .............')
print('results_dir: ',res_dir,' .............')
print('log_name_dir: ',log_name_dir,' .............')
print('solve_dir: ',solve_dir,' .............')


print('\n==============================================\nnum of elements: ',elem_num,'\n==============================================')  

for k in range(0,size(Um)):
    print('\n===============================     Um = ',Um[k],'      ===============================\n')
    cmd=['cp','-r',spectrum_dir, res_dir]
    outputs,errors=system_process(cmd,2000)
    cmd=['mv', res_dir+'energ_spectrum', solve_dir];
    outputs,errors=system_process(cmd,2000)
    
    print('\n=========================       CFL_const 90% simulation      =========================\n')  
    CFL = 0.9 * CFL_max;  
    dtt_ = CFL * dx / (Um[k]+1.5) 
    iter_p = round(0.01/dtt_);
    print('CFL: ', CFL,'  dt: ',dtt_,'  iter_p: ',iter_p)
    
    for lines in fileinput.input(fname, inplace=True):
        a = lines.split('=')
        b = a[0].strip()
        if b == 'velocity_mean':
            print('   {} = {}'.format(a[0].strip(),float(Um[k])))
        elif b == 'mode':
            print(' {} = {}'.format(a[0].strip(),'CFL_const'))
        elif b == 'CFL_no':
            print(' {} = {}'.format(a[0].strip(),CFL))
        elif b == 'unsteady_data_print_flag':
            print(' {} = {}'.format(a[0].strip(),2))
        elif b == 'unsteady_data_print_iter':
            #if k==0:
            print(' {} = {}'.format(a[0].strip(),iter_p))
            #elif k==1:
            #    print(' {} = {}'.format(a[0].strip(),50))
        elif b == 'calculate_dt_flag':
            print(' {} = {}'.format(a[0].strip(),1))
        elif b == 'final_time':
            print(' {} = {}'.format(a[0].strip(),'0.100'))
        else:
            print(lines.rstrip())
        
    for i in range(1,65):
        case_index = i   
        if case_index < 10:
            case_index_s = str('0')+str(i)
        else:
            case_index_s = str(i)
         
        for lines in fileinput.input(fname, inplace=True):
            a = lines.split('=')
            b = a[0].strip()
            if b == 'case_no':
                print(' {} = {}'.format(a[0].strip(),case_index_s))
            else:
                print(lines.rstrip())

        cmd = ['./bin/DG1DFlow.exe', fname]
        outputs,errors=system_process(cmd,5000)

        log_name = log_name_dir+case_index_s+str('/log_case')+ case_index_s + str('.out')
        if not(outputs is None):
            log = open(log_name,'a+')
            log.writelines(outputs)
        
        print('case_no: ',case_index_s)  
        
    #--------------------------------------------------------------------------------------------------------------------------------------------#
    print('\n=========================       CFL_const 50% simulation      =========================\n') 
    CFL = 0.5 * CFL_max;  
    dtt_ = CFL * dx / (Um[k]+1.5) 
    iter_p = round(0.01/dtt_);
    print('CFL: ', CFL,'  dt: ',dtt_,'  iter_p: ',iter_p)
    
    for lines in fileinput.input(fname, inplace=True):
        a = lines.split('=')
        b = a[0].strip()
        if b == 'CFL_no':
            print(' {} = {}'.format(a[0].strip(),CFL))
        elif b == 'unsteady_data_print_flag':
            print(' {} = {}'.format(a[0].strip(),2))
        elif b == 'unsteady_data_print_iter':
            print(' {} = {}'.format(a[0].strip(),iter_p))
            #if k==0:
            #    print(' {} = {}'.format(a[0].strip(),5))
            #elif k==1:
            #    print(' {} = {}'.format(a[0].strip(),50))
        elif b == 'calculate_dt_flag':
            print(' {} = {}'.format(a[0].strip(),1))
        elif b == 'final_time':
            print(' {} = {}'.format(a[0].strip(),'0.100'))
        else:
            print(lines.rstrip())
                
    for i in range(1,65):
        case_index = i   
        if case_index < 10:
            case_index_s = str('0')+str(i)
        else:
            case_index_s = str(i)
             
        for lines in fileinput.input(fname, inplace=True):
            a = lines.split('=')
            b = a[0].strip()
            if b == 'case_no':
                print(' {} = {}'.format(a[0].strip(),case_index_s))
            else:
                print(lines.rstrip())

        cmd = ['./bin/DG1DFlow.exe', fname]
        outputs,errors=system_process(cmd,5000)

        log_name = log_name_dir+case_index_s+str('/log_case')+ case_index_s + str('.out')
        if not(outputs is None):
            log = open(log_name,'a+')
            log.writelines(outputs)
           
        print('case_no: ',case_index_s)   

    #--------------------------------------------------------------------------------------------------------------------------------------------#
    print('\n=========================       dt_const, dt = 1.0e-05      =========================\n') 
    dtt_ = 1.000e-05;
    CFL = (Um[k]+1.5) * dtt_ / dx 
    iter_p = round(0.01/dtt_);
    print('CFL: ', CFL,'  dt: ',dtt_,'  iter_p: ',iter_p)

    for lines in fileinput.input(fname, inplace=True):
        a = lines.split('=')
        b = a[0].strip()
        if b == 'mode':
            print(' {} = {}'.format(a[0].strip(),'dt_const'))
        elif b == 'unsteady_data_print_flag':
            print(' {} = {}'.format(a[0].strip(),0))
        elif b == 'unsteady_data_print_iter':
            print(' {} = {}'.format(a[0].strip(),500))
        elif b == 'calculate_dt_flag':
            print(' {} = {}'.format(a[0].strip(),0))
        elif b == 'dt':
            print(' {} = {}'.format(a[0].strip(),dtt_))
        elif b == 'final_time':
            print(' {} = {}'.format(a[0].strip(),'0.250'))
        else:
            print(lines.rstrip())
                
    for i in range(1,65):
        case_index = i   
        if case_index < 10:
            case_index_s = str('0')+str(i)
        else:
            case_index_s = str(i)
             
        for lines in fileinput.input(fname, inplace=True):
            a = lines.split('=')
            b = a[0].strip()
            if b == 'case_no':
                print(' {} = {}'.format(a[0].strip(),case_index_s))
            else:
                print(lines.rstrip())

        cmd = ['./bin/DG1DFlow.exe', fname]
        outputs,errors=system_process(cmd,5000)

        log_name = log_name_dir+case_index_s+str('/log_case')+ case_index_s + str('.out')
        if not(outputs is None):
            log = open(log_name,'a+')
            log.writelines(outputs)
            
        print('case_no: ',case_index_s)
        
    new_dir = res_dir + 'p' + str(P) + 'RK' + str(Prk) +'_' \
        + str(n_pts) +'pts_um' + str(um_name[k]) + '/'
    cmd = ['mv',solve_dir,new_dir]
    outputs,errors=system_process(cmd,5000)
    
