import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from time import time
import sys
itime=time()
fig = plt.figure()
ax = fig.add_subplot(111)
cmap = plt.get_cmap("gnuplot")
color=lambda a:cmap(float(a))
l=2
N=18
Tmin=0.09
Tmax=0.1
Tsteps=10
T=np.linspace(Tmin,Tmax,Tsteps)
realizations=40
samples=100
repeats=1000
steps=300
SS=np.zeros((samples,steps,6),dtype=np.float_)
ene=np.zeros((samples,steps),dtype=np.float_)
ene2=np.zeros((samples,steps),dtype=np.float_)
norm=np.zeros((samples,steps),dtype=np.float_)
C=np.zeros((realizations,Tsteps),dtype=np.float_)
std_C=np.zeros((realizations,Tsteps),dtype=np.float_)
for realization in np.arange(realizations):
    print(time()-itime)
    for sample in np.arange(samples):
        SS[sample,:,:]=np.genfromtxt("zvo{}_SS_rand{}.dat".format(realization,sample))
        norm[sample,:]=np.cumprod(np.genfromtxt("zvo{}_Norm_rand{}.dat".format(realization,sample))[:,1])**2
        ene[sample,:]=SS[sample,:,1]*norm[sample,:]
        ene2[sample,:]=SS[sample,:,2]*norm[sample,:]
    index=0
    for t in T:
        C_t=np.zeros(repeats,dtype=np.float_)
        beta=1/t
        for repeat in np.arange(repeats):
            cnorm=0
            cene2=0
            cene=0
            id_b=np.random.choice(samples,samples,replace=True)
            factor=np.exp(-N*beta*l,dtype=np.float_)
            max_cnorm_k=0
            max_cene_k=0
            max_cene2_k=0
            for k in np.arange(steps):
                cnorm_k=((1+l*beta*N/(2*k+1))*norm[id_b,k]-beta/(2*k+1)*ene[id_b,k])*factor
                cene_k=((1+l*beta*N/(2*k+1))*ene[id_b,k]-beta/(2*k+1)*ene2[id_b,k])*factor
                cene2_k=((1-l*beta*N/(2*k+1))*ene2[id_b,k]+(l*l*beta*N*N/(2*k+1))*ene[id_b,k]-beta*N*N/(2*k+1)*ene[id_b,k+1])*factor
                #cene2_k=(ene2[id_b,k]+beta*N/(2*k+1)*(l*l*l*N*N*norm[id_b,k]-l*l*N*ene[id_b,k]-l*N*N*norm[id_b,k+1]-N*ene[id_b,k+1]))*factor
                max_cnorm_k=np.max(np.hstack((max_cnorm_k,np.abs(cnorm_k))))
                max_cene_k=np.max(np.hstack((max_cene_k,np.abs(cene_k))))
                max_cene2_k=np.max(np.hstack((max_cene2_k,np.abs(cene2_k))))
                if np.max(np.abs(cnorm_k)) < max_cnorm_k*sys.float_info.epsilon and np.max(np.abs(cene_k)) < max_cene_k*sys.float_info.epsilon and np.max(np.abs(cene2_k)) < max_cene2_k*sys.float_info.epsilon:
                    print("realization",realization,"repeat",repeat,"tmp",t,"cutoff",k)
                    break
                cnorm+=cnorm_k
                cene+=cene_k
                cene2+=cene2_k
                factor*=(N*beta)**2/(4*k*k+6*k+2)
            if k == steps-1:
                print("realization",realization,"repeat",repeat,"tmp",t,"unconveged!!")
                sys.exit()
            C_t[repeat]=np.mean(cene2)/np.mean(cnorm)-np.mean(cene)**2/(np.mean(cnorm)**2)
        C[realization,index]=beta*beta*np.mean(C_t)
        std_C[realization,index]=beta*beta*np.std(C_t,ddof=1)
        index+=1
print("finish",time()-itime)
yerr=np.sqrt(np.mean(std_C**2,axis=0)/realizations)/N
y=np.mean(C,axis=0)/N
np.save("CSHN18GC",np.vstack((T,y,yerr)).T)
ax.errorbar(T,y,yerr=yerr, fmt='x')
sc = ax.scatter(T, y, s=25, marker='x',linewidths=2,label="bootstrap")
ax.tick_params(labelsize=16)
ax.set_xlim(0,1.05*np.max(T))
ax.set_ylim(0,1.15*np.max(y))
ax.set_title(r'Specific Heat (N=18,$\Delta=1$)' ,size=16)
ax.set_xlabel(r'$T$', size=16)
ax.set_ylabel(r'$C$', size=16)
plt.tight_layout()
plt.show()
