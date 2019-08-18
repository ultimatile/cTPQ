struct mp
    l :: Float64
    N :: Int64
    realizations :: Int64
    samples :: Int64
    repeats :: Int64
    steps :: Int64
    Tmin :: Float64
    Tmax :: Float64
    Tsteps :: Int64
end

function main(mp)
    #repeats = 1 cause NaN of std(C_t)
    SS = zeros(mp.samples,mp.steps,6);ene = zeros(mp.samples,mp.steps);
    ene2 = zeros(mp.samples,mp.steps);norm = zeros(mp.samples,mp.steps)
    C = zeros(mp.realizations,mp.Tsteps);std_C = zeros(mp.realizations,mp.Tsteps)
    #eneGS=readdlm("energy.dat")
    eneGS = zeros(realizations)#
    T = range(Tmin,Tmax,length=Tsteps)
    for realization in 1:mp.realizations
        for sample in 1:mp.samples
            SS[sample,:,:] = readdlm("zvo$(realization-1)_SS_rand$(sample-1).dat",Float64,comments=true)
            norm[sample,:] = cumprod(readdlm("zvo$(realization-1)_Norm_rand$(sample-1).dat",Float64,comments=true)[:,2]).^2
            ene[sample,:] = SS[sample,:,2].*norm[sample,:]
            ene2[sample,:] = SS[sample,:,3].*norm[sample,:]
        end #for sample in 1:mp.samples
        for Tindex in 1:mp.Tsteps
            beta = 1/T[Tindex]
            C_t = zeros(mp.repeats)
            for repeat in 1:mp.repeats
                cnorm=zeros(mp.samples);cene2=zeros(mp.samples);cene=zeros(mp.samples)
                id_b = rand(1:mp.samples,mp.samples)
                factor = exp(BigFloat(-N*beta*l+beta*eneGS[realization]))
                max_cnorm_k = BigFloat(0);max_cene_k = BigFloat(0);max_cene2_k = BigFloat(0)
                ret=0
                for k in 0:mp.steps-2
                    cnorm_k = ((1+l*beta*mp.N/(2*k+1))*norm[id_b,k+1]-beta/(2*k+1)*ene[id_b,k+1])*factor
                    cene_k = ((1+l*beta*mp.N/(2*k+1))*ene[id_b,k+1]-beta/(2*k+1)*ene2[id_b,k+1])*factor
                    cene2_k = ((1-l*beta*mp.N/(2*k+1))*ene2[id_b,k+1]+(l*l*beta*mp.N*mp.N/(2*k+1))*ene[id_b,k+1]-beta*mp.N*mp.N/(2*k+1)*ene[id_b,k+2])*factor
                    #cene2_k=(ene2[id_b,k]+beta*N/(2*k+1)*(l*l*l*N*N*norm[id_b,k]-l*l*N*ene[id_b,k]-l*N*N*norm[id_b,k+1]-N*ene[id_b,k+1]))*factor
                    tmp_max_cnorm_k = maximum(abs.(cnorm_k))
                    tmp_max_cene_k = maximum(abs.(cene_k))
                    tmp_max_cene2_k = maximum(abs.(cene2_k))
                    if max_cnorm_k < tmp_max_cnorm_k;max_cnorm_k = tmp_max_cnorm_k;end
                    if max_cene_k < tmp_max_cene_k;max_cene_k = tmp_max_cene_k;end
                    if max_cene2_k < tmp_max_cene2_k;max_cene2_k = tmp_max_cene2_k;end
                    if tmp_max_cnorm_k < max_cnorm_k*eps(Float64) && tmp_max_cene_k < max_cene_k*eps(Float64) && tmp_max_cene2_k < max_cene2_k*eps(Float64)
                        println("realization:",realization," repeat:",repeat," t:",1/beta," cutoff:",k)
                        ret = 1
                        break
                    end
                    #println("$(max_cnorm_k) $(max_cene_k) $(max_cene2_k)")
                    cnorm += cnorm_k
                    cene += cene_k
                    cene2 += cene2_k
                    factor *= (mp.N*beta)^2/(4*k*k+6*k+2)
                end #for k in 0:steps-1
                #exit()
                if ret == 0;println("realization:",realization," repeat:",repeat," t:",1/beta," unconverged!!");end
                C_t[repeat] = mean(cene2)/mean(cnorm)-mean(cene)^2/(mean(cnorm)^2)
                @assert C_t[repeat] > 0 "C = $(C_t[repeat]) is negative!"
            end #repeat in 1:repeats
            mean_C_t = mean(C_t)
            C[realization,Tindex] = beta*beta*mean_C_t
            std_C[realization,Tindex] = beta*beta*std(C_t,corrected=true,mean=mean_C_t)
        end #for Tindex in 1:Tsteps
    end #for realization=0:realizations
    yerr = sqrt.(mean(std_C.^2,dims=1)/mp.realizations)/mp.N
    y = mean(C,dims=1)/mp.N
    #println(yerr)
    writedlm("CSHN$(mp.N)T$(mp.Tmin)_$(mp.Tmax).dat", [T y' yerr'])
end

using DelimitedFiles
using Statistics
const l=2
const N=18
const realizations=40
const samples=100
const repeats=5
const steps=300
const Tmin=0.2
const Tmax=20
const Tsteps=2
modpara=mp(l,N,realizations,samples,repeats,steps,Tmin,Tmax,Tsteps)
@time main(modpara)
