import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from time import time
#from scipy.misc import factorial
import sys
itime=time()
fig = plt.figure()
ax = fig.add_subplot(111)
cmap = plt.get_cmap("gnuplot")
color=lambda a:cmap(float(a))
l=2
N=18
Tmin=1
Tmax=10
Tsteps=10
T=np.linspace(Tmin,Tmax,Tsteps)
realizations=1
samples=10
repeats=10
#realizations=40
#samples=100
#repeats=500
steps=300
SS=np.zeros((samples,steps,6))
ene=np.zeros((samples,realizations,steps))
ene2=np.zeros((samples,realizations,steps))
norm=np.zeros((samples,realizations,steps))
pnorm=np.zeros((samples,realizations,steps))
C=np.zeros((realizations,Tsteps))
std_C=np.zeros((realizations,Tsteps))
for realization in np.arange(realizations):
    print(time()-itime)
    for sample in np.arange(samples):
        SS[sample,:,:]=np.genfromtxt("zvo{}_SS_rand{}.dat".format(realization,sample))#,skip_header=1)
        pnorm[sample,realization,:]=np.genfromtxt("zvo{}_Norm_rand{}.dat".format(realization,sample))[:,1]#,skip_header=1)
        norm[sample,realization,0]=pnorm[sample,realization,0]
        for i in np.arange(1,steps):
            norm[sample,realization,i]=norm[sample,realization,i-1]*pnorm[sample,realization,i]
            print(norm[sample,realization,i])
        ene[sample,realization,:]=SS[sample,:,1]*norm[sample,realization,:]
        ene2[sample,realization,:]=SS[sample,:,2]*norm[sample,realization,:]
    index=0

    for t in T:
        C_t=np.zeros(repeats)
        beta=1/t
        cnorm=0
        cene2=0
        cene=0
        for repeat in np.arange(repeats):
            id_b=np.random.choice(samples,samples,replace=True)
            #factor=(beta*N)**2/2
            factor=1
            for k in np.arange(steps-1):
                if factor >sys.float_info.epsilon:
                    print("cene2",cene2)
                    cnorm+=factor*((1+l*beta*N/(2*k+1))*norm[id_b,realization,k]-beta/(2*k+1)*ene[id_b,realization,k])
                    cene+=factor*((1+l*beta*N/(2*k+1))*ene[id_b,realization,k]-beta/(2*k+1)*ene2[id_b,realization,k])
                    cene2+=factor*((1-l*beta*N/(2*k+1))*ene2[id_b,realization,k]+(l*l*beta*N*N/(2*k+1))*ene[id_b,realization,k]-beta*N*N/(2*k+1)*ene[id_b,realization,k+1])
                factor*=(N*beta)**2/(4*k*k+6*k+2)
            C_t[repeat]=np.mean(cene2)/np.mean(cnorm)-np.mean(cene)**2/(np.mean(cnorm)**2)
            #print(C_t[repeat])
        C[realization,index]=beta*beta*np.mean(C_t)
        std_C[realization,index]=beta*beta*np.std(C_t,ddof=1)
        index+=1
#print(C)
print("finish")
print(time()-itime)
yerr=np.sqrt(np.sum(std_C**2,axis=0)/realizations)/realizations
y=np.mean(C,axis=0)
np.save("CSHN18GC",np.vstack((T,y,yerr)).T)
ax.errorbar(T, y,yerr=yerr, fmt='x')
sc = ax.scatter(T, y, s=25, marker='x',linewidths=2,label="bootstrap")
#x=1/np.mean(np.mean(beta,axis=0),axis=0)
#y=np.mean(np.mean(beta**2*(tene2-tene**2),axis=0),axis=0)/N
#sc = ax.scatter(x, y, s=25, marker='x',linewidths=2,label="standard")
#plt.legend()
ax.tick_params(labelsize=16)
#ax.set_xlim(0,1.05*np.max(T))
#ax.set_ylim(0,1.15*np.max(y))
#ax.set_ylim(0,0.1)
#sub_axes.set_ylim(0,0.02)
# タイトル
ax.set_title(r'Specific Heat (N=18,$\Delta=1$)' ,size=16)
#ax.set_xlabel(r'$1/N$', size=16)
ax.set_xlabel(r'$T$', size=16)
ax.set_ylabel(r'$C$', size=16)
plt.tight_layout()
plt.show()
