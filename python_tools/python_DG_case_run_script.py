
from time import sleep
from subprocess import call, run, PIPE, Popen
import fileinput
from numpy import size

elem_num = [256,512,1024,2048,4096,8192,16384,32768]  # OA4
#elem_num = [51,102,204,409,819,1638,3276,6553]  # p5

fname = './input/burgers_decay_turb_case_input.in'
log_name_dir = './Results_AdvecDiffus/Decaying_Burgers_turb/DGp3_RK3_test/case'

# Computing and Running the first no. of elements :
#--------------------------------------------------

for lines in fileinput.input(fname, inplace=True):
    a = lines.split('=')
    b = a[0].strip()
    if b == 'spectrum_restart_flag':
       print(' {} = {}'.format(a[0].strip(),0))
    elif b == 'num_elements':
        print(' {} = {}'.format(a[0].strip(),elem_num[0]))
    else:
       print(lines.rstrip())

print('\n==============================================\nnum of elements: ',elem_num[0],'\n==============================================')  
        
for i in range(1,33):
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

    cmd = ['./bin/FD1DFlow.exe', fname]
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

# Computing and Running other no. of elements :
#--------------------------------------------------
for lines in fileinput.input(fname, inplace=True):
    a = lines.split('=')
    b = a[0].strip()
    if b == 'spectrum_restart_flag':
       print(' {} = {}'.format(a[0].strip(),1))
    else:
       print(lines.rstrip())
       
for j in range(1,size(elem_num)):

    for lines in fileinput.input(fname, inplace=True):
        a = lines.split('=')
        b = a[0].strip()
        if b == 'num_elements':
            print(' {} = {}'.format(a[0].strip(),elem_num[j]))
        else:
            print(lines.rstrip())
    
    print('\n==============================================\nnum of elements: ',elem_num[j],'\n==============================================')             
    
    for i in range(1,33):
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

        cmd = ['./bin/FD1DFlow.exe', fname]
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


