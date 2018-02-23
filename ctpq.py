import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from time import time
from scipy.misc import factorial
import sys
itime=time()
fig = plt.figure()
ax = fig.add_subplot(111)
cmap = plt.get_cmap("gnuplot")
color=lambda a:cmap(float(a))
l=2
N=18
Tmin=1
Tmax=20
Tsteps=100
T=np.linspace(Tmin,Tmax,Tsteps)
#realizations=1
#samples=1
#repeats=1
realizations=40
samples=100
repeats=1000
steps=300
SS=np.zeros((samples,steps,6),dtype=np.float_)
ene=np.zeros((samples,realizations,steps),dtype=np.float_)
ene2=np.zeros((samples,realizations,steps),dtype=np.float_)
norm=np.zeros((samples,realizations,steps),dtype=np.float_)
pnorm=np.zeros((samples,realizations,steps),dtype=np.float_)
C=np.zeros((realizations,Tsteps),dtype=np.float_)
std_C=np.zeros((realizations,Tsteps),dtype=np.float_)
#old_settings = np.seterr(all='raise')
for realization in np.arange(realizations):
    print(time()-itime)
    for sample in np.arange(samples):
        SS[sample,:,:]=np.genfromtxt("zvo{}_SS_rand{}.dat".format(realization,sample))#,skip_header=1)
        norm[sample,realization,:]=np.cumprod(np.genfromtxt("zvo{}_Norm_rand{}.dat".format(realization,sample))[:,1])**2
        ene[sample,realization,:]=SS[sample,:,1]*(norm[sample,realization,:])
        ene2[sample,realization,:]=SS[sample,:,2]*(norm[sample,realization,:])
        #print(ene[sample,realization,:])
        #ene[sample,realization,:]=SS[sample,:,1]
        #ene2[sample,realization,:]=SS[sample,:,2]
    index=0
    for t in T:
        C_t=np.zeros(repeats,dtype=np.float_)
        beta=1/t
        for repeat in np.arange(repeats):
            cnorm=0
            cene2=0
            cene=0
            id_b=np.random.choice(samples,samples,replace=True)
            #print(id_b)
            #factor=(beta*N)**2/2
            factor=np.exp(-N*beta*l,dtype=np.float_)
            max_cnorm_k=0
            max_cene_k=0
            max_cene2_k=0
            for k in np.arange(steps-1):
                #print(l*l*N*ene[id_b,realization,k]-l*ene2[id_b,realization,k]-N*ene[id_b,realization,k+1])
                #print("cene",cene)
                #print("cene2",cene2)
                #print("cnorm",cnorm)
                cnorm_k=((1+l*beta*N/(2*k+1))*norm[id_b,realization,k]-beta/(2*k+1)*ene[id_b,realization,k])*factor
                cene_k=((1+l*beta*N/(2*k+1))*ene[id_b,realization,k]-beta/(2*k+1)*ene2[id_b,realization,k])*factor
                cene2_k=((1-l*beta*N/(2*k+1))*ene2[id_b,realization,k]+(l*l*beta*N*N/(2*k+1))*ene[id_b,realization,k]-beta*N*N/(2*k+1)*ene[id_b,realization,k+1])*factor
                #cene2_k=(ene2[id_b,realization,k]+beta*N/(2*k+1)*(l*l*l*N*N*norm[id_b,realization,k]-l*l*N*ene[id_b,realization,k]-l*N*N*norm[id_b,realization,k+1]-N*ene[id_b,realization,k+1]))*factor
                #print(k,cnorm_k,cene_k,cene2_k)
                #print(k,cene2_k,norm[id_b,realization,k],ene[id_b,realization,k],norm[id_b,realization,k+1],ene[id_b,realization,k+1],factor)
                max_cnorm_k=np.max(np.hstack((max_cnorm_k,np.abs(cnorm_k))))
                max_cene_k=np.max(np.hstack((max_cene_k,np.abs(cene_k))))
                max_cene2_k=np.max(np.hstack((max_cene2_k,np.abs(cene2_k))))
                #print("max",max_cnorm_k*sys.float_info.epsilon,max_cene_k*sys.float_info.epsilon,max_cene2_k*sys.float_info.epsilon)
                if np.max(np.abs(cnorm_k)) < max_cnorm_k*sys.float_info.epsilon and np.max(np.abs(cene_k)) < max_cene_k*sys.float_info.epsilon and np.max(np.abs(cene2_k)) < max_cene2_k*sys.float_info.epsilon:
                    break
                cnorm+=cnorm_k
                cene+=cene_k
                cene2+=cene2_k
                factor*=(N*beta)**2/(4*k*k+6*k+2)
            C_t[repeat]=np.mean(cene2)/np.mean(cnorm)-np.mean(cene)**2/(np.mean(cnorm)**2)
            #print(C_t[repeat])
        C[realization,index]=beta*beta*np.mean(C_t)
        std_C[realization,index]=beta*beta*np.std(C_t,ddof=1)
        index+=1
#print(C)
print("finish")
print(time()-itime)
yerr=np.sqrt(np.sum(std_C**2,axis=0)/realizations)/realizations/N
y=np.mean(C,axis=0)/N
np.save("CSHN18GC",np.vstack((T,y,yerr)).T)
ax.errorbar(T,y,yerr=yerr, fmt='x')
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
