function get_generator_data(data; use_gpu=false)
    ngen = length(data.generators)

    if use_gpu
        pgmin = CuArray{Float64}(undef, ngen)
        pgmax = CuArray{Float64}(undef, ngen)
        qgmin = CuArray{Float64}(undef, ngen)
        qgmax = CuArray{Float64}(undef, ngen)
        c2 = CuArray{Float64}(undef, ngen)
        c1 = CuArray{Float64}(undef, ngen)
        c0 = CuArray{Float64}(undef, ngen)
    else
        pgmin = Array{Float64}(undef, ngen)
        pgmax = Array{Float64}(undef, ngen)
        qgmin = Array{Float64}(undef, ngen)
        qgmax = Array{Float64}(undef, ngen)
        c2 = Array{Float64}(undef, ngen)
        c1 = Array{Float64}(undef, ngen)
        c0 = Array{Float64}(undef, ngen)
    end

    copyto!(pgmin, data.genvec.Pmin)
    copyto!(pgmax, data.genvec.Pmax)
    copyto!(qgmin, data.genvec.Qmin)
    copyto!(qgmax, data.genvec.Qmax)
    copyto!(c0, data.genvec.coeff0)
    copyto!(c1, data.genvec.coeff1)
    copyto!(c2, data.genvec.coeff2)

    return pgmin,pgmax,qgmin,qgmax,c2,c1,c0
end

function get_bus_data(data; use_gpu=false)
    ngen = length(data.generators)
    nbus = length(data.buses)
    nline = length(data.lines)

    FrIdx = [l for b=1:nbus for l in data.FromLines[b]]
    ToIdx = [l for b=1:nbus for l in data.ToLines[b]]
    GenIdx = [g for b=1:nbus for g in data.BusGenerators[b]]
    FrStart = accumulate(+, vcat([1], [length(data.FromLines[b]) for b=1:nbus]))
    ToStart = accumulate(+, vcat([1], [length(data.ToLines[b]) for b=1:nbus]))
    GenStart = accumulate(+, vcat([1], [length(data.BusGenerators[b]) for b=1:nbus]))

    Pd = [data.buses[i].Pd for i=1:nbus]
    Qd = [data.buses[i].Qd for i=1:nbus]

    if use_gpu
        cuFrIdx = CuArray{Int}(undef, length(FrIdx))
        cuToIdx = CuArray{Int}(undef, length(ToIdx))
        cuGenIdx = CuArray{Int}(undef, length(GenIdx))
        cuFrStart = CuArray{Int}(undef, length(FrStart))
        cuToStart = CuArray{Int}(undef, length(ToStart))
        cuGenStart = CuArray{Int}(undef, length(GenStart))
        cuPd = CuArray{Float64}(undef, nbus)
        cuQd = CuArray{Float64}(undef, nbus)

        copyto!(cuFrIdx, FrIdx)
        copyto!(cuToIdx, ToIdx)
        copyto!(cuGenIdx, GenIdx)
        copyto!(cuFrStart, FrStart)
        copyto!(cuToStart, ToStart)
        copyto!(cuGenStart, GenStart)
        copyto!(cuPd, Pd)
        copyto!(cuQd, Qd)

        return cuFrStart,cuFrIdx,cuToStart,cuToIdx,cuGenStart,cuGenIdx,cuPd,cuQd
    else
        return FrStart,FrIdx,ToStart,ToIdx,GenStart,GenIdx,Pd,Qd
    end
end

function get_branch_data(data; use_gpu=false)
    buses = data.buses
    lines = data.lines
    BusIdx = data.BusIdx
    nline = length(data.lines)
    ybus = Ybus{Array{Float64}}(computeAdmitances(data.lines, data.buses, data.baseMVA; VI=Array{Int}, VD=Array{Float64})...)
    frBound = [ x for l=1:nline for x in (buses[BusIdx[lines[l].from]].Vmin^2, buses[BusIdx[lines[l].from]].Vmax^2) ]
    toBound = [ x for l=1:nline for x in (buses[BusIdx[lines[l].to]].Vmin^2, buses[BusIdx[lines[l].to]].Vmax^2) ]

    if use_gpu
        cuYshR = CuArray{Float64}(undef, length(ybus.YshR))
        cuYshI = CuArray{Float64}(undef, length(ybus.YshI))
        cuYffR = CuArray{Float64}(undef, nline)
        cuYffI = CuArray{Float64}(undef, nline)
        cuYftR = CuArray{Float64}(undef, nline)
        cuYftI = CuArray{Float64}(undef, nline)
        cuYttR = CuArray{Float64}(undef, nline)
        cuYttI = CuArray{Float64}(undef, nline)
        cuYtfR = CuArray{Float64}(undef, nline)
        cuYtfI = CuArray{Float64}(undef, nline)
        cuFrBound = CuArray{Float64}(undef, 2*nline)
        cuToBound = CuArray{Float64}(undef, 2*nline)
        copyto!(cuYshR, ybus.YshR)
        copyto!(cuYshI, ybus.YshI)
        copyto!(cuYffR, ybus.YffR)
        copyto!(cuYffI, ybus.YffI)
        copyto!(cuYftR, ybus.YftR)
        copyto!(cuYftI, ybus.YftI)
        copyto!(cuYttR, ybus.YttR)
        copyto!(cuYttI, ybus.YttI)
        copyto!(cuYtfR, ybus.YtfR)
        copyto!(cuYtfI, ybus.YtfI)
        copyto!(cuFrBound, frBound)
        copyto!(cuToBound, toBound)

        return cuYshR, cuYshI, cuYffR, cuYffI, cuYftR, cuYftI,
               cuYttR, cuYttI, cuYtfR, cuYtfI, cuFrBound, cuToBound
    else
        return ybus.YshR, ybus.YshI, ybus.YffR, ybus.YffI, ybus.YftR, ybus.YftI,
               ybus.YttR, ybus.YttI, ybus.YtfR, ybus.YtfI, frBound, toBound
    end
end

function init_values(data, ybus, pg_start, qg_start, pij_start, qij_start,
                     pji_start, qji_start, wi_i_ij_start, wi_j_ji_start,
                     rho_pq, rho_va, u_curr, v_curr, l_curr, rho, wRIij)
    lines = data.lines
    buses = data.buses
    BusIdx = data.BusIdx
    ngen = length(data.generators)
    nline = length(data.lines)

    YffR = ybus.YffR; YffI = ybus.YffI
    YttR = ybus.YttR; YttI = ybus.YttI
    YftR = ybus.YftR; YftI = ybus.YftI
    YtfR = ybus.YtfR; YtfI = ybus.YtfI
    YshR = ybus.YshR; YshI = ybus.YshI

    for g=1:ngen
        v_curr[pg_start+g] = 0.5*(data.genvec.Pmin[g] + data.genvec.Pmax[g])
        v_curr[qg_start+g] = 0.5*(data.genvec.Qmin[g] + data.genvec.Qmax[g])
    end

    for l=1:nline
        wij0 = (buses[BusIdx[lines[l].from]].Vmax^2 + buses[BusIdx[lines[l].from]].Vmin^2) / 2
        wji0 = (buses[BusIdx[lines[l].to]].Vmax^2 + buses[BusIdx[lines[l].to]].Vmin^2) / 2

        wR0 = sqrt(wij0 * wji0)
        u_curr[pij_start+l] = YffR[l] * wij0 + YftR[l] * wR0
        u_curr[qij_start+l] = -YffI[l] * wij0 - YftI[l] * wR0
        u_curr[wi_i_ij_start+l] = wij0
        u_curr[pji_start+l] = YttR[l] * wji0 + YtfR[l] * wR0
        u_curr[qji_start+l] = -YttI[l] * wji0 - YtfI[l] * wR0
        u_curr[wi_j_ji_start+l] = wji0
        #wRIij[2*(l-1)+1] = wR0
        #wRIij[2*l] = 0.0

        v_curr[wi_i_ij_start+l] = 1.0
        v_curr[wi_j_ji_start+l] = 1.0
    end

    l_curr .= 0

    rho[1:2*ngen+4*nline] .= rho_pq
    rho[2*ngen+4*nline+1:end] .= rho_va

    return
end

function admm_rect_gpu(case; iterlim=800, rho_pq=400.0, rho_va=40000.0, use_gpu=false)
    data = opf_loaddata(case)

    ngen = length(data.generators)
    nline = length(data.lines)
    nbus = length(data.buses)
    nvar = 2*ngen + 6*nline

    baseMVA = data.baseMVA
    n = 8
    mu_max = 1e8
    rho_max = 1e6
    rho_min_pq = 5.0
    rho_min_w = 5.0
    eps_rp = 1e-4
    eps_rp_min = 1e-5
    rt_inc = 2
    rt_dec = 2
    eta = 0.99
    Kf = 100
    Kf_mean = 10

    ybus = Ybus{Array{Float64}}(computeAdmitances(data.lines, data.buses, data.baseMVA; VI=Array{Int}, VD=Array{Float64})...)

    pgmin, pgmax, qgmin, qgmax, c2, c1, c0 = get_generator_data(data)
    YshR, YshI, YffR, YffI, YftR, YftI, YttR, YttI, YtfR, YtfI, FrBound, ToBound = get_branch_data(data)
    FrStart, FrIdx, ToStart, ToIdx, GenStart, GenIdx, Pd, Qd = get_bus_data(data)

    cu_pgmin, cu_pgmax, cu_qgmin, cu_qgmax, cu_c2, cu_c1, cu_c0 = get_generator_data(data; use_gpu=true)
    cuYshR, cuYshI, cuYffR, cuYffI, cuYftR, cuYftI, cuYttR, cuYttI, cuYtfR, cuYtfI, cuFrBound, cuToBound = get_branch_data(data; use_gpu=true)
    cu_FrStart, cu_FrIdx, cu_ToStart, cu_ToIdx, cu_GenStart, cu_GenIdx, cu_Pd, cu_Qd = get_bus_data(data; use_gpu=true)

    pg_start = 0
    qg_start = ngen
    pij_start = 2*ngen
    qij_start = 2*ngen + nline
    pji_start = 2*ngen + 2*nline
    qji_start = 2*ngen + 3*nline
    wi_i_ij_start = 2*ngen + 4*nline
    wi_j_ji_start = 2*ngen + 5*nline

    u_curr = zeros(nvar)
    v_curr = zeros(nvar)
    l_curr = zeros(nvar)
    u_prev = zeros(nvar)
    v_prev = zeros(nvar)
    l_prev = zeros(nvar)
    rho = zeros(nvar)
    rd = zeros(nvar)
    rp = zeros(nvar)
    rp_old = zeros(nvar)
    rp_k0 = zeros(nvar)
    param = zeros(24, nline)
    wRIij = zeros(2*nline)

    init_values(data, ybus, pg_start, qg_start, pij_start, qij_start,
                pji_start, qji_start, wi_i_ij_start, wi_j_ji_start,
                rho_pq, rho_va, u_curr, v_curr, l_curr, rho, wRIij)

    cu_u_curr = CuArray{Float64}(undef, nvar)
    cu_v_curr = CuArray{Float64}(undef, nvar)
    cu_l_curr = CuArray{Float64}(undef, nvar)
    cu_u_prev = CuArray{Float64}(undef, nvar)
    cu_v_prev = CuArray{Float64}(undef, nvar)
    cu_l_prev = CuArray{Float64}(undef, nvar)
    cu_rho = CuArray{Float64}(undef, nvar)
    cu_rd = CuArray{Float64}(undef, nvar)
    cu_rp = CuArray{Float64}(undef, nvar)
    cu_rp_old = CuArray{Float64}(undef, nvar)
    cu_rp_k0 = CuArray{Float64}(undef, nvar)
    cuParam = CuArray{Float64}(undef, (24, nline))
    cuWRIij = CuArray{Float64}(undef, 2*nline)

    copyto!(cu_u_curr, u_curr)
    copyto!(cu_v_curr, v_curr)
    copyto!(cu_l_curr, l_curr)
    copyto!(cu_rho, rho)
    copyto!(cuParam, param)
    copyto!(cuWRIij, wRIij)

    max_auglag = 50

    nblk_gen = div(ngen, 32, RoundUp)
    nblk_br = nline
    nblk_bus = div(nbus, 32, RoundUp)

    ABSTOL = 1e-6
    RELTOL = 1e-5

    it = 0
    time_gen = time_br = time_bus = 0

    h_u_curr = zeros(nvar)
    h_param = zeros(24, nline)
    h_wRIij = zeros(2*nline)

    while it < iterlim
        it += 1

        if !use_gpu
            u_prev .= u_curr
            v_prev .= v_curr
            l_prev .= l_curr

            generator_kernel_cpu(baseMVA, ngen, pg_start, qg_start, u_curr, v_curr, l_curr, rho,
                                pgmin, pgmax, qgmin, qgmax, c2, c1, c0)
            auglag_it = auglag_kernel_cpu(n, nline, it, max_auglag, pij_start, qij_start, pji_start, qji_start,
                                        wi_i_ij_start, wi_j_ji_start, mu_max,
                                        u_curr, v_curr, l_curr, rho,
                                        wRIij, param, YffR, YffI, YftR, YftI,
                                        YttR, YttI, YtfR, YtfI, FrBound, ToBound)
            bus_kernel_cpu(baseMVA, nbus, pg_start, qg_start, pij_start, qij_start,
                                pji_start, qji_start, wi_i_ij_start, wi_j_ji_start,
                                FrStart, FrIdx, ToStart, ToIdx, GenStart,
                                GenIdx, Pd, Qd, u_curr, v_curr, l_curr, rho, YshR, YshI)

            l_curr .+= rho .* (u_curr .- v_curr)
            rd .= -rho .* (v_curr .- v_prev)
            rp .= u_curr .- v_curr
            rp_old .= u_prev .- v_prev

            primres = norm(rp)
            dualres = norm(rd)

            eps_pri = sqrt(length(l_curr))*ABSTOL + RELTOL*max(norm(u_curr), norm(-v_curr))
            eps_dual = sqrt(length(u_curr))*ABSTOL + RELTOL*norm(l_curr)

            @printf("[CPU] %10d  %.6e  %.6e  %.6e  %.6e  %.2f\n", it, primres, dualres, eps_pri, eps_dual, auglag_it)
        else
            cu_u_prev .= cu_u_curr
            cu_v_prev .= cu_v_curr
            cu_l_prev .= cu_l_curr

            tgpu = CUDA.@timed @cuda threads=32 blocks=nblk_gen generator_kernel(baseMVA, ngen, pg_start, qg_start,
                                                                                    cu_u_curr, cu_v_curr, cu_l_curr, cu_rho,
                                                                                    cu_pgmin, cu_pgmax, cu_qgmin, cu_qgmax, cu_c2, cu_c1, cu_c0)
            time_gen += tgpu.time
            tgpu = CUDA.@timed @cuda threads=(n,n) blocks=nblk_br shmem=(sizeof(Float64)*(14*n+3*n^2) + sizeof(Int)*(4*n)) auglag_kernel(n, it, max_auglag, pij_start, qij_start, pji_start, qji_start,
                                                                                                                        wi_i_ij_start, wi_j_ji_start, mu_max,
                                                                                                                        cu_u_curr, cu_v_curr, cu_l_curr, cu_rho,
                                                                                                                        cuWRIij, cuParam, cuYffR, cuYffI, cuYftR, cuYftI,
                                                                                                                        cuYttR, cuYttI, cuYtfR, cuYtfI, cuFrBound, cuToBound)
            time_br += tgpu.time
            tgpu = CUDA.@timed @cuda threads=32 blocks=nblk_bus bus_kernel(baseMVA, nbus, pg_start, qg_start, pij_start, qij_start,
                                                                           pji_start, qji_start, wi_i_ij_start, wi_j_ji_start,
                                                                           cu_FrStart, cu_FrIdx, cu_ToStart, cu_ToIdx, cu_GenStart,
                                                                           cu_GenIdx, cu_Pd, cu_Qd, cu_u_curr, cu_v_curr, cu_l_curr,
                                                                           cu_rho, cuYshR, cuYshI)
            time_bus += tgpu.time

            cu_l_curr .+= cu_rho .* (cu_u_curr .- cu_v_curr)
            cu_rd .= -cu_rho .* (cu_v_curr .- cu_v_prev)
            cu_rp .= cu_u_curr .- cu_v_curr
            cu_rp_old .= cu_u_prev .- cu_v_prev

            gpu_primres = CUDA.norm(cu_rp)
            gpu_dualres = CUDA.norm(cu_rd)

            gpu_eps_pri = sqrt(length(l_curr))*ABSTOL + RELTOL*max(CUDA.norm(cu_u_curr), CUDA.norm(-cu_v_curr))
            gpu_eps_dual = sqrt(length(u_curr))*ABSTOL + RELTOL*CUDA.norm(cu_l_curr)

            @printf("[GPU] %10d  %.6e  %.6e  %.6e  %.6e\n", it, gpu_primres, gpu_dualres, gpu_eps_pri, gpu_eps_dual)
        end
    end

    if use_gpu
        copyto!(u_curr, cu_u_curr)
        @printf(" ** [GPU] Kernel time\n")
        @printf("Generator = %.2f\n", time_gen)
        @printf("Branch    = %.2f\n", time_br)
        @printf("Bus       = %.2f\n", time_bus)
        @printf("Total     = %.2f\n", time_gen + time_br + time_bus)
    end

    objval = sum(data.generators[g].coeff[data.generators[g].n-2]*(baseMVA*u_curr[pg_start + g])^2 +
                 data.generators[g].coeff[data.generators[g].n-1]*(baseMVA*u_curr[pg_start + g]) +
                 data.generators[g].coeff[data.generators[g].n]
                 for g in 1:ngen)
    @printf("Objective value = %.6e\n", objval)

    return
end