import os
from time import sleep
import sys
from subprocess import call, PIPE, Popen
import fileinput
from numpy import size

elem_num = 256  # OA4
#elem_num = [51,102,204,409,819,1638,3276,6553]  # p5

n_pts = [4,5,15,20];
Um=[0,1,10,25,50,75,100];

dir_result = './Results_AdvecDiffus/decay_burg_turb2/DGp3_RK3_s/'
solve_dir = './Results_AdvecDiffus/decay_burg_turb2/DGp3_RK3/'
new_dir = './Results_AdvecDiffus/decay_burg_turb2/p3_RK3_'+str(n_pts[1])+'pts_um'+ str(Um[5]) + '/'
cmd = [str('cp -r ')+dir_result+new_dir]
call(['cp','-r',dir_result, solve_dir]);
call(['mv',solve_dir,new_dir])
