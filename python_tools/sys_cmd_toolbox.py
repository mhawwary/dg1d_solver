import os
from time import sleep
import sys
from subprocess import call, PIPE, Popen
import fileinput
from numpy import size

def system_process(cmd,time_to_wait):

    process = Popen(cmd, stdin=PIPE, stdout=PIPE , start_new_session=True)
    try:
        outs, errs = process.communicate(timeout=time_to_wait)

    except TimeoutExpired:
        process.kill()
        outs, errs = process.communicate()
        
    if not(outs is None):
        output = outs.decode("utf-8")
    else:
        output = None
    
    if not(errs is None):
        errors = errs.decode("utf-8")
       #err_f = open('errs.out','w')
       #err_f.writelines(errors)
    else:
        errors = None
           
    process.terminate()
    
    return output,errors 
