import os
from time import sleep
import sys
from subprocess import call, PIPE, Popen
import fileinput
from numpy import size

from sys_cmd_toolbox import system_process

elem_num = 256  # OA4
#elem_num = [51,102,204,409,819,1638,3276,6553]  # p5

P=6;
n_pts = [4*P,3*P,2*P,P+1]
Um=[10.0,25.0,50.0,100.0]
um_name = [10,25,50,100]
dir_temp = './Results_AdvecDiffus/decay_burg_turb2/DGp3_RK3_s/'
solve_dir = './Results_AdvecDiffus/decay_burg_turb2/DGp3_RK3/'

fname = './input/burgers_turb_case_input.in'
log_name_dir = './Results_AdvecDiffus/decay_burg_turb2/DGp3_RK3/case'

# Computing and Running the first no. of elements :
#--------------------------------------------------

for lines in fileinput.input(fname, inplace=True):
    a = lines.split('=')
    b = a[0].strip()
    if b == 'spectrum_restart_flag':
       print('   {} = {}'.format(a[0].strip(),1))
    elif b == 'num_elements':
        print(' {} = {}'.format(a[0].strip(),elem_num))
    elif b == 'mode':
        print(' {} = {}'.format(a[0].strip(),'CFL_const'))
    elif b == 'CFL_no':
        print(' {} = {}'.format(a[0].strip(),0.1170))
    elif b == 'unsteady_data_print_flag':
        print(' {} = {}'.format(a[0].strip(),1))
    elif b == 'calculate_dt_flag':
        print(' {} = {}'.format(a[0].strip(),1))
    else:
        print(lines.rstrip())

print('\n==============================================\nnum of elements: ',elem_num,'\n==============================================')  

for k in range(0,size(Um)):
    for lines in fileinput.input(fname, inplace=True):
        a = lines.split('=')
        b = a[0].strip()
        if b == 'velocity_mean':
            print('   {} = {}'.format(a[0].strip(),Um[k]))
        else:
            print(lines.rstrip())
    
    for j in range(0,size(n_pts)):
        for lines in fileinput.input(fname, inplace=True):
            a = lines.split('=')
            b = a[0].strip()
            if b == 'N_uniform_pts_per_elem':
                print('   {} = {}'.format(a[0].strip(),n_pts[j]))
            else:
                print(lines.rstrip())
            
        print('\n========================= Um=',Um[k],', n_pts=',n_pts[j],' =========================\n') 
        
        cmd=['cp','-r',dir_temp, solve_dir];
        output=system_process(cmd,1000)
        
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
            output=system_process(cmd,5000)

            log_name = log_name_dir+case_index_s+str('/log_case')+ case_index_s + str('.out')
            if not(output is None):
                log = open(log_name,'a+')
                log.writelines(output)
            
            print('case_no: ',case_index_s)
        
        new_dir = './Results_AdvecDiffus/decay_burg_turb2/p3_RK3_'+str(n_pts[j])+'pts_um'+ str(Um[k]) + '/'
        cmd = ['mv',solve_dir,new_dir]
        output=system_process(cmd,1000)



