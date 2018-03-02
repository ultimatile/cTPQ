l=2;N=18
Tmin=1;Tmax=2;Tsteps=100
T=linspace(Tmin,Tmax,Tsteps)
itime=time()
realizations=1;samples=100;repeats=1000;steps=300
#realizations=1;samples=100;repeats=2;steps=300
#repeats=1 cause NaN of std(C_t)
SS=zeros(samples,steps,6);ene=zeros(samples,steps);
ene2=zeros(samples,steps);norm=zeros(samples,steps)
C=zeros(realizations,Tsteps);std_C=zeros(realizations,Tsteps)
for realization in 1:realizations
    for sample in 1:samples
        SS[sample,:,:]=readdlm("zvo$(realization-1)_SS_rand$(sample-1).dat")
        norm[sample,:]=cumprod(readdlm("zvo$(realization-1)_Norm_rand$(sample-1).dat")[:,2]).^2
        ene[sample,:]=SS[sample,:,2].*norm[sample,:]
        ene2[sample,:]=SS[sample,:,3].*norm[sample,:]
    end #for sample in 0:samples
    for Tindex in 1:Tsteps
        C_t=zeros(repeats)
        beta=1/T[Tindex]
        for repeat in 1:repeats
            cnorm=0;cene2=0;cene=0
            id_b= rand(1:samples,samples)
            factor=exp(-N*beta*l)
            max_cnorm_k=0;max_cene_k=0;max_cene2_k=0
            for k in 0:steps-2
                cnorm_k=((1+l*beta*N/(2*k+1))*norm[id_b,k+1]-beta/(2*k+1)*ene[id_b,k+1])*factor
                cene_k=((1+l*beta*N/(2*k+1))*ene[id_b,k+1]-beta/(2*k+1)*ene2[id_b,k+1])*factor
                cene2_k=((1-l*beta*N/(2*k+1))*ene2[id_b,k+1]+(l*l*beta*N*N/(2*k+1))*ene[id_b,k+1]-beta*N*N/(2*k+1)*ene[id_b,k+2])*factor
                #cene2_k=(ene2[id_b,k]+beta*N/(2*k+1)*(l*l*l*N*N*norm[id_b,k]-l*l*N*ene[id_b,k]-l*N*N*norm[id_b,k+1]-N*ene[id_b,k+1]))*factor
                tmp_max_cnorm_k=maximum(abs.(cnorm_k))
                tmp_max_cene_k=maximum(abs.(cene_k))
                tmp_max_cene2_k=maximum(abs.(cene2_k))
                if max_cnorm_k<tmp_max_cnorm_k;max_cnorm_k=tmp_max_cnorm_k;end
                if max_cene_k<tmp_max_cene_k;max_cene_k=tmp_max_cene_k;end
                if max_cene2_k<tmp_max_cene2_k;max_cene2_k=tmp_max_cene2_k;end
                if tmp_max_cnorm_k< max_cnorm_k*eps(Float64) && tmp_max_cene_k< max_cene_k*eps(Float64) && tmp_max_cene2_k< max_cene2_k*eps(Float64)
                    println("realization:",realization," repeat:",repeat," t:",1/beta," cutoff:",k)
                    break
                end
                #println("$(max_cnorm_k) $(max_cene_k) $(max_cene2_k)")
                cnorm+=cnorm_k
                cene+=cene_k
                cene2+=cene2_k
                factor*=(N*beta)^2/(4*k*k+6*k+2)
            end #for k in 0:steps-1
            C_t[repeat]=mean(cene2)/mean(cnorm)-mean(cene)^2/(mean(cnorm)^2)
        end#repeat in 1:repeats
        mean_C_t=mean(C_t)
        C[realization,Tindex]=beta*beta*mean_C_t
        std_C[realization,Tindex]=beta*beta*std(C_t,corrected=true,mean=mean_C_t)
    end#for Tindex in 1:Tsteps
end#for realization=0:realizations
yerr=sqrt.(mean(std_C.^2,1)/realizations)/N
y=mean(C,1)/N
println("finish:$(time()-itime)")
#println(yerr)
writedlm("SHN$(N)T$(Tmin)_$(Tmax).dat", [T y' yerr'])
