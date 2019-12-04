using ArgParse
using DelimitedFiles
using Statistics

function parse_commandline()
    s = ArgParseSettings()
    s.description = "This is for observation of cTPQ by mTPQ data obtained by HPhi."
    @add_arg_table s begin
    "--calibration", "-c"
        help = "an option for determine the minimum temperature"
        action = :store_true
    "--Tmin"
        default = 0.5
        arg_type = Float64
    "--Tmax"
        default = 2.0
        arg_type = Float64
    "--energy", "-e"
        help = "with GS energy"
        action = :store_true
    end
    return parse_args(s)
end

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
    flag_car :: Bool
    flag_ene :: Bool
end

function main(mp)
    #repeats = 1 cause NaN of std(C_t)
    SS = zeros(mp.samples, mp.steps, 6); ene = zeros(mp.samples, mp.steps)
    ene2 = zeros(mp.samples, mp.steps); norm = zeros(mp.samples, mp.steps)
    C = zeros(mp.realizations, mp.Tsteps); std_C = zeros(mp.realizations, mp.Tsteps)
    if mp.flag_ene
        eneGS = readdlm("energy.dat")
    else
        eneGS = zeros(mp.realizations)
    end
    T = range(mp.Tmin, mp.Tmax, length = mp.Tsteps)
    for realization in 1:mp.realizations
        for sample in 1:mp.samples
            SS[sample, :, :] = readdlm("zvo$(realization - 1)_SS_rand$(sample - 1).dat", Float64, comments = true)
            norm[sample, :] = cumprod(readdlm("zvo$(realization - 1)_Norm_rand$(sample - 1).dat", Float64, comments = true)[:, 2]) .^ 2
            ene[sample, :] = SS[sample, :, 2] .* norm[sample, :]
            ene2[sample, :] = SS[sample, :, 3] .* norm[sample, :]
        end #for sample in 1:mp.samples
        for Tindex in 1:mp.Tsteps
            beta = 1 / T[Tindex]
            C_t = zeros(mp.repeats)
            for repeat in 1:mp.repeats
                cnorm = zeros(mp.samples); cene2 = zeros(mp.samples); cene = zeros(mp.samples)
                id_b = rand(1:mp.samples, mp.samples)
                factor = exp(BigFloat(-mp.N * beta * mp.l + beta * eneGS[realization]))
                max_cnorm_k = BigFloat(0); max_cene_k = BigFloat(0); max_cene2_k = BigFloat(0)
                ret = 0
                for k in 0:mp.steps - 2
                    cnorm_k = ((1 + mp.l * beta * mp.N / (2 * k + 1)) * norm[id_b, k + 1] - beta/(2 * k + 1) * ene[id_b, k + 1]) * factor
                    cene_k = ((1 + mp.l * beta * mp.N / (2 * k + 1)) * ene[id_b, k + 1] - beta / (2 * k + 1) * ene2[id_b, k + 1]) * factor
                    cene2_k = ((1 - mp.l * beta * mp.N / (2 * k + 1)) * ene2[id_b, k + 1] + (mp.l * mp.l * beta * mp.N * mp.N / (2 * k + 1)) 
                    * ene[id_b, k + 1] - beta * mp.N * mp.N / (2 * k + 1) * ene[id_b, k + 2]) * factor
                    #cene2_k=(ene2[id_b,k]+beta*N/(2*k+1)*(l*l*l*N*N*norm[id_b,k]-l*l*N*ene[id_b,k]-l*N*N*norm[id_b,k+1]-N*ene[id_b,k+1]))*factor
                    tmp_max_cnorm_k = maximum(abs.(cnorm_k))
                    tmp_max_cene_k = maximum(abs.(cene_k))
                    tmp_max_cene2_k = maximum(abs.(cene2_k))
                    if max_cnorm_k < tmp_max_cnorm_k; max_cnorm_k = tmp_max_cnorm_k; end
                    if max_cene_k < tmp_max_cene_k; max_cene_k = tmp_max_cene_k; end
                    if max_cene2_k < tmp_max_cene2_k; max_cene2_k = tmp_max_cene2_k; end
                    if tmp_max_cnorm_k < max_cnorm_k * eps(Float64) && tmp_max_cene_k < max_cene_k * eps(Float64) && tmp_max_cene2_k < max_cene2_k * eps(Float64)
                        println("realization:", realization, " repeat:", repeat, " t:", 1/beta, " cutoff:", k)
                        ret = 1
                        break
                    end
                    #println("$(max_cnorm_k) $(max_cene_k) $(max_cene2_k)")
                    cnorm += cnorm_k
                    cene += cene_k
                    cene2 += cene2_k
                    factor *= (mp.N * beta) ^ 2 / (4 * k * k + 6 * k + 2)
                end #for k in 0:steps-1
                if ret == 0; println("realization:", realization, " repeat:", repeat, " t:", 1 / beta, " unconverged!!"); end
                C_t[repeat] = mean(cene2) / mean(cnorm) - mean(cene) ^ 2 / (mean(cnorm) ^ 2)
                @assert C_t[repeat] > 0 "C = $(C_t[repeat]) is negative!"
            end #repeat in 1:repeats
            mean_C_t = mean(C_t)
            C[realization, Tindex] = beta * beta * mean_C_t
            std_C[realization, Tindex] = beta * beta * std(C_t, corrected = true, mean = mean_C_t)
        end #for Tindex in 1:Tsteps
    end #for realization=0:realizations
    yerr = sqrt.(mean(std_C .^ 2, dims = 1) / mp.realizations) / mp.N
    y = mean(C, dims = 1) / mp.N
    if mp.flag_car; return; end
    writedlm("CSHN$(mp.N)T$(mp.Tmin)_$(mp.Tmax).dat", [T y' yerr'])
end


#l = 0.894202698 #L30J1.0g0.6
l =  0.8932157447#L18J1.0g1.0
#l =  0.8937822841#L18J1.0g0.6
#l =  0.9623002067#L18J1.0g0.0

N = 18
realizations = 40
samples = 100
steps = 1000

parsed_args = parse_commandline()
if parsed_args["calibration"]
    Tsteps = 2
    repeats = 10
else
    Tsteps = 100
    repeats = 250
end

modpara = mp(l, N, realizations, samples, repeats, steps, parsed_args["Tmin"], parsed_args["Tmax"], Tsteps, parsed_args["calibration"], parsed_args["energy"])
@time main(modpara)
