using ArgParse
using DelimitedFiles
using Statistics

function parse_commandline()
    s = ArgParseSettings()
    s.description = "This is for observation of cTPQ by mTPQ data obtained by HPhi."
    @add_arg_table! s begin
    "--calibration", "-c"
        help = "an option for determine the minimum temperature"
        action = :store_true
    "--Tmin"
        default = 0.5
        arg_type = Float64
    "--Tmax"
        default = 2.0
        arg_type = Float64
    "--hmin"
        default = 0.0
        arg_type = Float64
    "--hmax"
        default = 5.0
        arg_type = Float64
    "--T", "-T"
        default = 2.0
        arg_type = Float64
    "--energy", "-e"
        help = "with GS energy"
        action = :store_true
    "--Sz"
        help = "Sz divided calculation"
        action = :store_true
    "--mag", "-m"
        help = "magnetizaiton calculation"
        action = :store_true
    end
    return parse_args(s)
end

struct mp
    l :: Float64
    N :: Int64
    twoS :: Int64
    realizations :: Int64
    samples :: Int64
    repeats :: Int64
    steps :: Int64
    Tmin :: Float64
    Tmax :: Float64
    Tsteps :: Int64
    flag_cal :: Bool
    flag_ene :: Bool
    flag_Sz :: Bool
end

struct mp_mag
    l :: Float64
    N :: Int64
    twoS :: Int64
    T :: Float64
    realizations :: Int64
    samples :: Int64
    repeats :: Int64
    steps :: Int64
    hmin :: Float64
    hmax :: Float64
    hsteps :: Int64
    flag_cal :: Bool
    flag_ene :: Bool
end

mutable struct Intermediates
    norm :: Array{Float64}
    ene :: Array{Float64}
    ene2 :: Array{Float64}
    function Intermediates(mp, realization)
        ene = zeros(mp.samples, mp.steps)
        ene2 = zeros(mp.samples, mp.steps)
        norm = zeros(mp.samples, mp.steps)
        for sample in 1:mp.samples
            if mp.realizations == 1
                try
                    SS = readdlm("SS_rand$(sample - 1).dat", Float64, comments = true)
                    norm[sample, :] = cumprod(readdlm("Norm_rand$(sample - 1).dat", Float64, comments = true)[:, 2]) .^ 2
                catch
                    SS = readdlm("zvo$(realization - 1)_SS_rand$(sample - 1).dat", Float64, comments = true)
                    norm[sample, :] = cumprod(readdlm("zvo$(realization - 1)_Norm_rand$(sample - 1).dat", Float64, comments = true)[:, 2]) .^ 2
                end
            else
                SS = readdlm("zvo$(realization - 1)_SS_rand$(sample - 1).dat", Float64, comments = true)
                norm[sample, :] = cumprod(readdlm("zvo$(realization - 1)_Norm_rand$(sample - 1).dat", Float64, comments = true)[:, 2]) .^ 2
            end
            ene[sample, :] = SS[:, 2] .* norm[sample, :]
            ene2[sample, :] = SS[:, 3] .* norm[sample, :]
        end #for sample in 1:mp.samples
        new(norm, ene, ene2)
    end
end

mutable struct Intermediates_Sz
    norm :: Array{Float64}
    ene :: Array{Float64}
    ene2 :: Array{Float64}
    function Intermediates_Sz(mp, realization)
        ene = zeros(div(mp.N * mp.twoS, 2) + 1, mp.samples, mp.steps)
        ene2 = zeros(div(mp.N * mp.twoS, 2) + 1, mp.samples, mp.steps)
        norm = zeros(div(mp.N * mp.twoS, 2) + 1, mp.samples, mp.steps)
        iSz = 1
        for twoSz in mod(mp.twoS * mp.N, 2):2:mp.twoS * mp.N
            for sample in 1:mp.samples
                if mp.realizations == 1
                    try
                        SS = readdlm("SS_rand$(sample - 1)_twoSz$(twoSz).dat", Float64, comments = true)
                        norm[iSz, sample, :] = cumprod(readdlm("Norm_rand$(sample - 1)_twoSz$(twoSz).dat", Float64, comments = true)[:, 2]) .^ 2
                    catch
                        SS = readdlm("zvo$(realization - 1)_SS_rand$(sample - 1)_twoSz$(twoSz).dat", Float64, comments = true)
                        norm[iSz, sample, :] .= cumprod(readdlm("zvo$(realization - 1)_Norm_rand$(sample - 1)_twoSz$(twoSz).dat", Float64, comments = true)[:, 2]) .^ 2
                    end
                else
                    SS = readdlm("zvo$(realization - 1)_SS_rand$(sample - 1)_twoSz$(twoSz).dat", Float64, comments = true)
                    norm[iSz, sample, :] = cumprod(readdlm("zvo$(realization - 1)_Norm_rand$(sample - 1)_twoSz$(twoSz).dat", Float64, comments = true)[:, 2]) .^ 2
                end
                ene[iSz, sample, :] = SS[:, 2] .* norm[iSz, sample, :]
                ene2[iSz, sample, :] = SS[:, 3] .* norm[iSz, sample, :]
            end #for sample in 1:mp.samples
            iSz += 1
        end# for twoSz in -twoS*N:2:twoS*N
        new(norm, ene, ene2)
    end
end


function calc_canon(mp, intr, eneGS, beta)
    cnorm = zeros(mp.samples); cene2 = zeros(mp.samples); cene = zeros(mp.samples)
    id_b = rand(1:mp.samples, mp.samples)
    factor = exp(BigFloat(-mp.N * beta * mp.l + beta * eneGS))
    max_cnorm_k = BigFloat(0); max_cene_k = BigFloat(0); max_cene2_k = BigFloat(0)
    for k in 0:mp.steps - 2
        cnorm_k = ((2k + 1 + mp.l * beta * mp.N) * intr.norm[id_b, k + 1] - beta * intr.ene[id_b, k + 1]) * factor
        cene_k = ((2k + 1 + mp.l * beta * mp.N) * intr.ene[id_b, k + 1] - beta * intr.ene2[id_b, k + 1]) * factor
        cene2_k = ((2k + 1 - mp.l * beta * mp.N) * intr.ene2[id_b, k + 1] + (mp.l * mp.l * beta * mp.N * mp.N)
        * intr.ene[id_b, k + 1] - beta * mp.N * mp.N * intr.ene[id_b, k + 2]) * factor
        #cene2_k=(ene2[id_b,k]+beta*N/(2*k+1)*(l*l*l*N*N*norm[id_b,k]-l*l*N*ene[id_b,k]-l*N*N*norm[id_b,k+1]-N*ene[id_b,k+1]))*factor
        tmp_max_cnorm_k = maximum(abs.(cnorm_k))
        tmp_max_cene_k = maximum(abs.(cene_k))
        tmp_max_cene2_k = maximum(abs.(cene2_k))
        if max_cnorm_k < tmp_max_cnorm_k; max_cnorm_k = tmp_max_cnorm_k; end
        if max_cene_k < tmp_max_cene_k; max_cene_k = tmp_max_cene_k; end
        if max_cene2_k < tmp_max_cene2_k; max_cene2_k = tmp_max_cene2_k; end
        if tmp_max_cnorm_k < max_cnorm_k * eps(Float64) && tmp_max_cene_k < max_cene_k * eps(Float64) && tmp_max_cene2_k < max_cene2_k * eps(Float64)
            return cnorm, cene, cene2, k, true
        end
        #println("$(max_cnorm_k) $(max_cene_k) $(max_cene2_k)")
        cnorm += cnorm_k
        cene += cene_k
        cene2 += cene2_k
        factor *= (mp.N * beta) ^ 2 / (4k ^ 2 + 6k + 2)
    end #for k in 0:steps-1
    return cnorm, cene, cene2, k, false
end

function calc_canon(mp, intr, eneGS, beta, iSz)
    cnorm = zeros(mp.samples); cene2 = zeros(mp.samples); cene = zeros(mp.samples)
    id_b = rand(1:mp.samples, mp.samples)
    factor = exp(BigFloat(-mp.N * beta * mp.l + beta * eneGS))
    max_cnorm_k = BigFloat(0); max_cene_k = BigFloat(0); max_cene2_k = BigFloat(0)
    for k in 0:mp.steps - 2
        cnorm_k = ((2k + 1 + mp.l * beta * mp.N ) * intr.norm[iSz, id_b, k + 1] - beta * intr.ene[iSz, id_b, k + 1]) * factor
        cene_k = ((2k + 1 + mp.l * beta * mp.N) * intr.ene[iSz, id_b, k + 1] - beta * intr.ene2[iSz, id_b, k + 1]) * factor
        cene2_k = ((2k + 1 - mp.l * beta * mp.N) * intr.ene2[iSz, id_b, k + 1] + (mp.l * mp.l * beta * mp.N * mp.N)
        * intr.ene[iSz, id_b, k + 1] - beta * mp.N * mp.N * intr.ene[iSz, id_b, k + 2]) * factor
        tmp_max_cnorm_k = maximum(abs.(cnorm_k))
        tmp_max_cene_k = maximum(abs.(cene_k))
        tmp_max_cene2_k = maximum(abs.(cene2_k))
        if max_cnorm_k < tmp_max_cnorm_k; max_cnorm_k = tmp_max_cnorm_k; end
        if max_cene_k < tmp_max_cene_k; max_cene_k = tmp_max_cene_k; end
        if max_cene2_k < tmp_max_cene2_k; max_cene2_k = tmp_max_cene2_k; end
        if tmp_max_cnorm_k < max_cnorm_k * eps(Float64) && tmp_max_cene_k < max_cene_k * eps(Float64) && tmp_max_cene2_k < max_cene2_k * eps(Float64)
            return cnorm, cene, cene2, k, true
        end
        #println("$(max_cnorm_k) $(max_cene_k) $(max_cene2_k)")
        cnorm += cnorm_k
        cene += cene_k
        cene2 += cene2_k
        factor *= (mp.N * beta) ^ 2 / (4 * k * k + 6 * k + 2)
    end #
    return cnorm, cene, cene2, k, false
end

function calc_canon_mag(mp, intr, eneGS, beta, iSz)
    cmag = zeros(mp.samples)
    id_b = rand(1:mp.samples, mp.samples)
    factor = exp(BigFloat(-mp.N * beta * mp.l + beta * eneGS))
    max_cmag_k = BigFloat(0)
    for k in 0:mp.steps - 2
        cmag_k = ((1 + mp.l * beta * mp.N / (2k + 1)) * intr.norm[iSz, id_b, k + 1] - beta / (2k + 1) * intr.ene[iSz, id_b, k + 1]) * factor
        tmp_max_cmag_k = maximum(abs.(cmag_k))
        if max_cmag_k < tmp_max_cmag_k; max_cmag_k = tmp_max_cmag_k; end
        if tmp_max_cmag_k < max_cmag_k * eps(Float64)
            return cmag, k, true
        end
        #println("$(max_cnorm_k) $(max_cene_k) $(max_cene2_k)")
        cmag += cmag_k
        factor *= (mp.N * beta) ^ 2 / (4 * k * k + 6 * k + 2)
    end #
    return cmag, k, false
end

function check_converge(ret, realization, repeat, beta, k)
    if ret
        println("realization:", realization, " repeat:", repeat, " t:", 1 / beta, " cutoff:", k)
    else
        println("realization:", realization, " repeat:", repeat, " t:", 1 / beta, " unconverged!!")
    end
end

function main_Sz(mp)
    #repeats = 1 cause NaN of std(C_t)
    C = zeros(mp.realizations, mp.Tsteps); std_C = zeros(mp.realizations, mp.Tsteps)
    χ = zeros(mp.realizations, mp.Tsteps); std_χ = zeros(mp.realizations, mp.Tsteps)
    if mp.flag_ene
        eneGS = readdlm("energy.dat")
    else
        eneGS = zeros(mp.realizations)
    end
    
    T = range(mp.Tmin, mp.Tmax, length = mp.Tsteps)
    
    for realization in 1:mp.realizations
        intr = Intermediates_Sz(mp, realization)    
        for Tindex in 1:mp.Tsteps
            beta = 1 / T[Tindex]
            C_t = zeros(mp.repeats)
            χ_t = zeros(mp.repeats)
            for repeat in 1:mp.repeats
                iSz = 1
                cnorm = zeros(mp.samples); cene2 = zeros(mp.samples); cene = zeros(mp.samples)
                for twoSz in mod(mp.twoS * mp.N, 2):2:mp.twoS * mp.N
                    tmp_cnorm, tmp_cene, tmp_cene2, k, ret = calc_canon(mp, intr, eneGS[realization], beta, iSz)
                    check_converge(ret, realization, repeat, beta, k)
                    mean_tmp_cnorm = mean(tmp_cnorm)
                    if twoSz == 0
                        cnorm += tmp_cnorm
                        cene += tmp_cene
                        cene2 += tmp_cene2
                    else
                        cnorm += 2 * tmp_cnorm
                        cene += 2 * tmp_cene
                        cene2 += 2 * tmp_cene2
                        χ_t[repeat] += twoSz * twoSz * mean_tmp_cnorm
                        @assert χ_t[repeat] >= 0 "χ = $(χ_t[repeat]) is negative!"
                    end
                    iSz += 1
                end
                mean_cnorm = mean(cnorm)
                C_t[repeat] = mean(cene2) / mean_cnorm - mean(cene) ^ 2 / (mean_cnorm ^ 2)
                @assert C_t[repeat] > 0 "C = $(C_t[repeat]) is negative!"
            end #repeat in 1:repeats
            mean_C_t = mean(C_t)
            C[realization, Tindex] = beta * beta * mean_C_t
            std_C[realization, Tindex] = beta * beta * std(C_t, corrected = true, mean = mean_C_t)
            χ_t　/= 4
            mean_χ_t = mean(χ_t)
            χ[realization, Tindex] = beta * mean_χ_t
            std_χ[realization, Tindex] = beta * std(χ_t, corrected = true, mean = mean_χ_t)
        end #for Tindex in 1:Tsteps
    end #for realization=0:realizations
    if mp.flag_cal; return; end
    Cave = mean(C, dims = 1) / mp.N 
    Cerr = sqrt.(mean(std_C .^ 2, dims = 1) / mp.realizations) / mp.N
    χave = mean(χ, dims = 1) / mp.N
    χerr = sqrt.(mean(std_χ .^ 2, dims = 1) / mp.realizations) / mp.N
    writedlm("CN$(mp.N)T$(mp.Tmin)_$(mp.Tmax)_Sz.dat", [T C' Cerr'])
    writedlm("chiN$(mp.N)T$(mp.Tmin)_$(mp.Tmax)_Sz.dat", [T χ' χerr'])
end

function main_mag(mp)
    #repeats = 1 cause NaN of std(C_t)
    m = zeros(mp.realizations, mp.Tsteps); std_m = zeros(mp.realizations, mp.Tsteps)
    if mp.flag_ene
        eneGS = readdlm("energy.dat")
    else
        eneGS = zeros(mp.realizations)
    end
    
    h = range(mp.hmin, mp.hmax, length = mp.hsteps)
    
    for realization in 1:mp.realizations
        intr = Intermediates_Sz(mp, realization)    
        for hindex in 1:mp.hsteps
            current_h =  h[hindex]
            beta = 1 / mp.T
            m_t = zeros(mp.repeats)
            for repeat in 1:mp.repeats
                iSz = 1
                cnorm = zeros(mp.samples); cene = zeros(mp.samples)
                for twoSz in mod(mp.twoS * mp.N, 2):2:mp.twoS * mp.N
                    tmp_cnorm, tmp_cene, k, ret = calc_canon_mag(mp, intr, eneGS[realization], beta, iSz)
                    check_converge(ret, realization, repeat, beta, k)
                    mean_tmp_cnorm = mean(tmp_cnorm)
                    if twoSz == 0
                        cnorm += tmp_cnorm
                        cene += tmp_cene
                    else
                        cnorm += 2 * tmp_cnorm
                        cene += 2 * tmp_cene
                        m_t[repeat] += twoSz * mean_tmp_cnorm
                        @assert χ_t[repeat] >= 0 "χ = $(χ_t[repeat]) is negative!"
                    end
                    iSz += 1
                end
                mean_cnorm = mean(cnorm)
                C_t[repeat] = mean(cene2) / mean_cnorm - mean(cene) ^ 2 / (mean_cnorm ^ 2)
                @assert C_t[repeat] > 0 "C = $(C_t[repeat]) is negative!"
            end #repeat in 1:repeats
            mean_C_t = mean(C_t)
            C[realization, Tindex] = beta * beta * mean_C_t
            std_C[realization, Tindex] = beta * beta * std(C_t, corrected = true, mean = mean_C_t)
            χ_t　/= 4
            mean_χ_t = mean(χ_t)
            χ[realization, Tindex] = beta * mean_χ_t
            std_χ[realization, Tindex] = beta * std(χ_t, corrected = true, mean = mean_χ_t)
        end #for Tindex in 1:Tsteps
    end #for realization=0:realizations
    if mp.flag_cal; return; end
    Cave = mean(C, dims = 1) / mp.N 
    Cerr = sqrt.(mean(std_C .^ 2, dims = 1) / mp.realizations) / mp.N
    χave = mean(χ, dims = 1) / mp.N
    χerr = sqrt.(mean(std_χ .^ 2, dims = 1) / mp.realizations) / mp.N
    writedlm("CN$(mp.N)T$(mp.Tmin)_$(mp.Tmax)_Sz.dat", [T C' Cerr'])
    writedlm("chiN$(mp.N)T$(mp.Tmin)_$(mp.Tmax)_Sz.dat", [T χ' χerr'])
end

function main_fullSz(mp)
    #repeats = 1 cause NaN of std(C_t)
    C = zeros(mp.realizations, mp.Tsteps); std_C = zeros(mp.realizations, mp.Tsteps)

    if mp.flag_ene
        eneGS = readdlm("energy.dat")
    else
        eneGS = zeros(mp.realizations)
    end

    T = range(mp.Tmin, mp.Tmax, length = mp.Tsteps)

    for realization in 1:mp.realizations
        intr = Intermediates(mp, realization)

        for Tindex in 1:mp.Tsteps
            beta = 1 / T[Tindex]
            C_t = zeros(mp.repeats)
            for repeat in 1:mp.repeats
                cnorm, cene, cene2, k, ret = calc_canon(mp, intr, eneGS[realization], beta)
                check_converge(ret, realization, repeat, beta, k)
                mean_cnorm = mean(cnorm)
                C_t[repeat] = mean(cene2) / mean_cnorm - mean(cene) ^ 2 / (mean_cnorm ^ 2)
                @assert C_t[repeat] > 0 "C = $(C_t[repeat]) is negative!"
            end #repeat in 1:repeats
            mean_C_t = mean(C_t)
            C[realization, Tindex] = beta * beta * mean_C_t
            std_C[realization, Tindex] = beta * beta * std(C_t, corrected = true, mean = mean_C_t)
        end #for Tindex in 1:Tsteps
    end #for realization=0:realizations
    if mp.flag_cal; return; end
    Cave = mean(C, dims = 1) / mp.N
    Cerr = sqrt.(mean(std_C .^ 2, dims = 1) / mp.realizations) / mp.N
    writedlm("CN$(mp.N)T$(mp.Tmin)_$(mp.Tmax).dat", [T C' Cerr'])
end

function main(mp)
    if mp.flag_Sz
        main_Sz(mp)
    else
        main_fullSz(mp)
    end
end

#l = 0.894202698 #L6J1.0g0.6
#l = 0.894202698 #L30J1.0g0.6
#l =  0.8932157447#L18J1.0g1.0
#l =  0.8937822841#L18J1.0g0.6
#l =  0.9623002067#L18J1.0g0.0

l = 2.00001
N = 6
realizations = 1
samples = 100
steps = 1000
twoS = 2

parsed_args = parse_commandline()

if parsed_args["mag"]
    if parsed_args["calibration"]
        hsteps = 2
        repeats = 10
    else
        hsteps = 50
        repeats = 250
    end
    modpara = mp_mag(l, N, twoS, T, realizations, samples, repeats, steps, parsed_args["hmin"], parsed_args["hmax"], hsteps, parsed_args["calibration"], parsed_args["energy"])
    @time main_mag(modpara)
else
    if parsed_args["calibration"]
        Tsteps = 2
        repeats = 10
    else
        Tsteps = 50
        repeats = 250
    end
    modpara = mp(l, N, twoS, realizations, samples, repeats, steps, parsed_args["Tmin"], parsed_args["Tmax"], Tsteps, parsed_args["calibration"], parsed_args["energy"], parsed_args["Sz"])
    @time main(modpara)
end