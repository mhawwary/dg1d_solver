import os
from time import sleep
import sys
from subprocess import call, PIPE, Popen
import fileinput
from numpy import size

elem_num = 256  # OA4
#elem_num = [51,102,204,409,819,1638,3276,6553]  # p5

fname = './input/burgers_turb_case_input.in'
log_name_dir = './Results_AdvecDiffus/decay_burg_turb/DGp3_RK3/case'

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
        print(' {} = {}'.format(a[0].strip(),'dt_const'))
    elif b == 'unsteady_data_print_flag':
        print(' {} = {}'.format(a[0].strip(),0))
    elif b == 'calculate_dt_flag':
        print(' {} = {}'.format(a[0].strip(),0))
    else:
       print(lines.rstrip())

print('\n==============================================\nnum of elements: ',elem_num,'\n==============================================')  

print('\n=========================       dt_const simulation      =========================\n')  
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
    process = Popen(cmd, stdin=PIPE, stdout=PIPE,start_new_session=True)

    try:
        outs, errs = process.communicate(timeout=5000)
        if not(outs is None):
            output = outs.decode("utf-8")
            #print(output)

    except TimeoutExpired:
        process.kill()
        outs, errs = process.communicate()
        if not(outs is None):
            output = outs.decode("utf-8")
          #  print(output)
           
    process.terminate()

    #if not(errs is None):
     #   errors = errs.decode("utf-8")
      #  print('errors,\n',errors)
       # err_f = open('errs.out','w')
        #err_f.writelines(errors)

    log_name = log_name_dir+case_index_s+str('/log_case')+ case_index_s + str('.out')
    if not(outs is None):
        log = open(log_name,'a+')
        log.writelines(output)
    
    print('case_no: ',case_index_s)
    
#--------------------------------------------------------------------------------------------------------------------------------------------#

for lines in fileinput.input(fname, inplace=True):
        a = lines.split('=')
        b = a[0].strip()
        if b == 'mode':
            print(' {} = {}'.format(a[0].strip(),'CFL_const'))
        elif b == 'CFL_no':
            print(' {} = {}'.format(a[0].strip(),0.1170))
        elif b == 'unsteady_data_print_flag':
            print(' {} = {}'.format(a[0].strip(),1))
        elif b == 'calculate_dt_flag':
            print(' {} = {}'.format(a[0].strip(),1))
        else:
            print(lines.rstrip())
            
print('\n=========================       CFL_const 90% simulation      =========================\n') 
             
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
    process = Popen(cmd, stdin=PIPE, stdout=PIPE,start_new_session=True)

    try:
        outs, errs = process.communicate(timeout=5000)
        if not(outs is None):
            output = outs.decode("utf-8")
            #print(output)

    except TimeoutExpired:
        process.kill()
        outs, errs = process.communicate()
        if not(outs is None):
            output = outs.decode("utf-8")
          #  print(output)
           
    process.terminate()

    #if not(errs is None):
     #   errors = errs.decode("utf-8")
      #  print('errors,\n',errors)
       # err_f = open('errs.out','w')
        #err_f.writelines(errors)

    log_name = log_name_dir+case_index_s+str('/log_case')+ case_index_s + str('.out')
    if not(outs is None):
        log = open(log_name,'a+')
        log.writelines(output)
    
    print('case_no: ',case_index_s)


#--------------------------------------------------------------------------------------------------------------------------------------------#

for lines in fileinput.input(fname, inplace=True):
        a = lines.split('=')
        b = a[0].strip()
        if b == 'CFL_no':
            print(' {} = {}'.format(a[0].strip(),0.0650))
        else:
            print(lines.rstrip())
            
print('\n=========================       CFL_const 50% simulation      =========================\n') 
             
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
    process = Popen(cmd, stdin=PIPE, stdout=PIPE,start_new_session=True)

    try:
        outs, errs = process.communicate(timeout=5000)
        if not(outs is None):
            output = outs.decode("utf-8")
            #print(output)

    except TimeoutExpired:
        process.kill()
        outs, errs = process.communicate()
        if not(outs is None):
            output = outs.decode("utf-8")
          #  print(output)
           
    process.terminate()

    #if not(errs is None):
     #   errors = errs.decode("utf-8")
      #  print('errors,\n',errors)
       # err_f = open('errs.out','w')
        #err_f.writelines(errors)

    log_name = log_name_dir+case_index_s+str('/log_case')+ case_index_s + str('.out')
    if not(outs is None):
        log = open(log_name,'a+')
        log.writelines(output)
    
    print('case_no: ',case_index_s)

