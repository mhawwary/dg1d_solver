from matplotlib import pyplot ,ticker   #and the useful plotting library
from numpy import sin,cos,pi,linspace,ones,zeros,abs,min,max,exp, shape\
    , empty_like , size, loadtxt, savetxt, arange , log, sqrt, real, transpose

from numpy import random, fft, array, append, hstack, vstack, sort, argsort, argmin, argmax,sum, logspace,log10

ko = 10.0

n_pts=2048

kmax=2048

data = loadtxt('spectrum_data_N2048.dat')

k = data[:,0]
basic_epsi  = data[:,1]
E_dumped = data[:,2]
#k = linspace(1,kmax,n_pts)

del data

data = loadtxt('u_cont_N4096_dt5.000e-06_Beta1.00_Eps1.00_0.000t.dat')
x_init = data[:,0]
u_init = data[:,1]

n_pts = size(k)

A = 2 * ko**-5 / (3 * sqrt(pi) )

E = A * k**4 * exp(-(k/ko)**2)

#basic_epsi = random.random_sample(n_pts)
u_basic_epsi = sqrt(2*E) * exp (1j * 2*pi * basic_epsi)
u_max_k = sqrt(2*max(E))
u_tot = hstack((0, u_basic_epsi))
u_x = fft.irfft(u_tot)*n_pts
x1 = 2*pi*arange(0,1.,1./(2*kmax))

u_f = fft.rfft(u_x)/n_pts
E_f = 0.5*abs(u_f)**2

#====== Use a cosine function =========#
x= 2*pi*arange(0,1.,1./(2*kmax))
N_pts = size(x)
ux=zeros(N_pts)
for i in range(0,N_pts):
    temp = sqrt(2*E) * cos(k*x[i]+2*pi*basic_epsi)
    ux[i]= sum(temp)

uf = fft.rfft(ux)/(0.5*N_pts)
Ef = 0.5*abs(uf)**2
print(size(E_f))

pyplot.figure();

pyplot.plot(k,E,'--b',label='original')
pyplot.plot(k,Ef[1:n_pts+1],'-r',label='recover_cosine')
pyplot.plot(k,E_f[1:],'-k',label='recover')
pyplot.plot(k,E_dumped,'-ok',label='dumped')

pyplot.xscale('log')
pyplot.yscale('log')

pyplot.grid()
pyplot.legend()

pyplot.xlabel('k', labelpad=10);
pyplot.ylabel('E(k)', labelpad=10);

pyplot.ylim(10**-12,10**0)

pyplot.title('E using cosines')

pyplot.show()


#======================================#

output = transpose([x, ux])

#savetxt('u_init.dat',output,fmt="%1.5f"+" %1.5f",delimiter=' ')


pyplot.figure();

pyplot.plot(x1,u_x,'-b',label='numerical')
pyplot.plot(x,ux,'--r',label='cosines')
pyplot.plot(x_init,u_init,'--k',label='initial_proj')

pyplot.grid()

pyplot.legend()

pyplot.title('u(x) comparison')

pyplot.xlabel('x', labelpad=10);
pyplot.ylabel('u(x)', labelpad=10);


pyplot.show()

