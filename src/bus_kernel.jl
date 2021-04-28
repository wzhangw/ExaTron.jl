function bus_kernel(
    baseMVA, nbus, gen_start, line_start,
    FrStart, FrIdx, ToStart, ToIdx, GenStart, GenIdx,
    Pd, Qd, u, v, l, rho, YshR, YshI
)
    I = threadIdx().x + (blockDim().x * (blockIdx().x - 1))
    T = eltype(l)
    if I <= nbus
        common_wi = zero(T)
        common_ti = zero(T)
        inv_rhosum_pij_ji = zero(T)
        inv_rhosum_qij_ji = zero(T)
        rhosum_wi_ij_ji = zero(T)
        rhosum_ti_ij_ji = zero(T)

        @inbounds begin
            for k=FrStart[I]:FrStart[I+1]-1
                pij_idx = line_start + 8*(FrIdx[k]-1)
                common_wi += l[pij_idx+4] + rho[pij_idx+4]*u[pij_idx+4]
                common_ti += l[pij_idx+6] + rho[pij_idx+6]*u[pij_idx+6]
                inv_rhosum_pij_ji += one(T) / rho[pij_idx]
                inv_rhosum_qij_ji += one(T) / rho[pij_idx+1]
                rhosum_wi_ij_ji += rho[pij_idx+4]
                rhosum_ti_ij_ji += rho[pij_idx+6]
            end
            for k=ToStart[I]:ToStart[I+1]-1
                pij_idx = line_start + 8*(ToIdx[k]-1)
                common_wi += l[pij_idx+5] + rho[pij_idx+5]*u[pij_idx+5]
                common_ti += l[pij_idx+7] + rho[pij_idx+7]*u[pij_idx+7]
                inv_rhosum_pij_ji += one(T) / rho[pij_idx+2]
                inv_rhosum_qij_ji += one(T) / rho[pij_idx+3]
                rhosum_wi_ij_ji += rho[pij_idx+5]
                rhosum_ti_ij_ji += rho[pij_idx+7]
            end
        end

        common_wi /= rhosum_wi_ij_ji

        rhs1 = zero(T)
        rhs2 = zero(T)
        inv_rhosum_pg = zero(T)
        inv_rhosum_qg = zero(T)

        @inbounds begin
            for g=GenStart[I]:GenStart[I+1]-1
                pg_idx = gen_start + 2*(GenIdx[g]-1)
                rhs1 += u[pg_idx] + (l[pg_idx]/rho[pg_idx])
                rhs2 += u[pg_idx+1] + (l[pg_idx+1]/rho[pg_idx+1])
                inv_rhosum_pg += one(T) / rho[pg_idx]
                inv_rhosum_qg += one(T) / rho[pg_idx+1]
            end

            rhs1 -= (Pd[I] / baseMVA)
            rhs2 -= (Qd[I] / baseMVA)

            for k=FrStart[I]:FrStart[I+1]-1
                pij_idx = line_start + 8*(FrIdx[k]-1)
                rhs1 -= u[pij_idx] + (l[pij_idx]/rho[pij_idx])
                rhs2 -= u[pij_idx+1] + (l[pij_idx+1]/rho[pij_idx+1])
            end

            for k=ToStart[I]:ToStart[I+1]-1
                pij_idx = line_start + 8*(ToIdx[k]-1)
                rhs1 -= u[pij_idx+2] + (l[pij_idx+2]/rho[pij_idx+2])
                rhs2 -= u[pij_idx+3] + (l[pij_idx+3]/rho[pij_idx+3])
            end

            rhs1 -= YshR[I]*common_wi
            rhs2 += YshI[I]*common_wi

            A11 = (inv_rhosum_pg + inv_rhosum_pij_ji) + (YshR[I]^2 / rhosum_wi_ij_ji)
            A12 = -YshR[I]*(YshI[I] / rhosum_wi_ij_ji)
            A21 = A12
            A22 = (inv_rhosum_qg + inv_rhosum_qij_ji) + (YshI[I]^2 / rhosum_wi_ij_ji)
            mu2 = (rhs2 - (A21/A11)*rhs1) / (A22 - (A21/A11)*A12)
            mu1 = (rhs1 - A12*mu2) / A11
            #mu = A \ [rhs1 ; rhs2]
            wi = common_wi + ( (YshR[I]*mu1 - YshI[I]*mu2) / rhosum_wi_ij_ji )
            ti = common_ti / rhosum_ti_ij_ji

            for k=GenStart[I]:GenStart[I+1]-1
                pg_idx = gen_start + 2*(GenIdx[k]-1)
                v[pg_idx] = u[pg_idx] + (l[pg_idx] - mu1) / rho[pg_idx]
                v[pg_idx+1] = u[pg_idx+1] + (l[pg_idx+1] - mu2) / rho[pg_idx+1]
            end
            for j=FrStart[I]:FrStart[I+1]-1
                pij_idx = line_start + 8*(FrIdx[j]-1)
                v[pij_idx] = u[pij_idx] + (l[pij_idx] + mu1) / rho[pij_idx]
                v[pij_idx+1] = u[pij_idx+1] + (l[pij_idx+1] + mu2) / rho[pij_idx+1]
                v[pij_idx+4] = wi
                v[pij_idx+6] = ti
            end
            for j=ToStart[I]:ToStart[I+1]-1
                pij_idx = line_start + 8*(ToIdx[j]-1)
                v[pij_idx+2] = u[pij_idx+2] + (l[pij_idx+2] + mu1) / rho[pij_idx+2]
                v[pij_idx+3] = u[pij_idx+3] + (l[pij_idx+3] + mu2) / rho[pij_idx+3]
                v[pij_idx+5] = wi
                v[pij_idx+7] = ti
            end
        end
    end
end

function bus_kernel_cpu(
    baseMVA, nbus, gen_start, line_start,
    FrStart, FrIdx, ToStart, ToIdx, GenStart, GenIdx,
    Pd, Qd, u, v, l, rho, YshR, YshI
)
    T = eltype(l)
    Threads.@threads for I=1:nbus
        common_wi = T(0.0)
        common_ti = T(0.0)
        inv_rhosum_pij_ji = T(0.0)
        inv_rhosum_qij_ji = T(0.0)
        rhosum_wi_ij_ji = T(0.0)
        rhosum_ti_ij_ji = T(0.0)

        @inbounds begin
            if FrStart[I] < FrStart[I+1]
                for k=FrStart[I]:FrStart[I+1]-1
                    pij_idx = line_start + 8*(FrIdx[k]-1)
                    common_wi += l[pij_idx+4] + rho[pij_idx+4]*u[pij_idx+4]
                    common_ti += l[pij_idx+6] + rho[pij_idx+6]*u[pij_idx+6]
                    inv_rhosum_pij_ji += one(T) / rho[pij_idx]
                    inv_rhosum_qij_ji += one(T) / rho[pij_idx+1]
                    rhosum_wi_ij_ji += rho[pij_idx+4]
                    rhosum_ti_ij_ji += rho[pij_idx+6]
                end
            end

            if ToStart[I] < ToStart[I+1]
                for k=ToStart[I]:ToStart[I+1]-1
                    pij_idx = line_start + 8*(ToIdx[k]-1)
                    common_wi += l[pij_idx+5] + rho[pij_idx+5]*u[pij_idx+5]
                    common_ti += l[pij_idx+7] + rho[pij_idx+7]*u[pij_idx+7]
                    inv_rhosum_pij_ji += one(T) / rho[pij_idx+2]
                    inv_rhosum_qij_ji += one(T) / rho[pij_idx+3]
                    rhosum_wi_ij_ji += rho[pij_idx+5]
                    rhosum_ti_ij_ji += rho[pij_idx+7]
                end
            end
        end

        common_wi /= rhosum_wi_ij_ji

        rhs1 = zero(T)
        rhs2 = zero(T)
        inv_rhosum_pg = zero(T)
        inv_rhosum_qg = zero(T)

        @inbounds begin
            if GenStart[I] < GenStart[I+1]
                for g=GenStart[I]:GenStart[I+1]-1
                    pg_idx = gen_start + 2*(GenIdx[g]-1)
                    rhs1 += u[pg_idx] + (l[pg_idx]/rho[pg_idx])
                    rhs2 += u[pg_idx+1] + (l[pg_idx+1]/rho[pg_idx+1])
                    inv_rhosum_pg += one(T) / rho[pg_idx]
                    inv_rhosum_qg += one(T) / rho[pg_idx+1]
                end
            end

            rhs1 -= (Pd[I] / baseMVA)
            rhs2 -= (Qd[I] / baseMVA)

            if FrStart[I] < FrStart[I+1]
                for k=FrStart[I]:FrStart[I+1]-1
                    pij_idx = line_start + 8*(FrIdx[k]-1)
                    rhs1 -= u[pij_idx] + (l[pij_idx]/rho[pij_idx])
                    rhs2 -= u[pij_idx+1] + (l[pij_idx+1]/rho[pij_idx+1])
                end
            end

            if ToStart[I] < ToStart[I+1]
                for k=ToStart[I]:ToStart[I+1]-1
                    pij_idx = line_start + 8*(ToIdx[k]-1)
                    rhs1 -= u[pij_idx+2] + (l[pij_idx+2]/rho[pij_idx+2])
                    rhs2 -= u[pij_idx+3] + (l[pij_idx+3]/rho[pij_idx+3])
                end
            end

            rhs1 -= YshR[I]*common_wi
            rhs2 += YshI[I]*common_wi

            A11 = (inv_rhosum_pg + inv_rhosum_pij_ji) + (YshR[I]^2 / rhosum_wi_ij_ji)
            A12 = -YshR[I]*(YshI[I] / rhosum_wi_ij_ji)
            A21 = A12
            A22 = (inv_rhosum_qg + inv_rhosum_qij_ji) + (YshI[I]^2 / rhosum_wi_ij_ji)
            mu2 = (rhs2 - (A21/A11)*rhs1) / (A22 - (A21/A11)*A12)
            mu1 = (rhs1 - A12*mu2) / A11
            #mu = A \ [rhs1 ; rhs2]
            wi = common_wi + ( (YshR[I]*mu1 - YshI[I]*mu2) / rhosum_wi_ij_ji )
            ti = common_ti / rhosum_ti_ij_ji

            for k=GenStart[I]:GenStart[I+1]-1
                pg_idx = gen_start + 2*(GenIdx[k]-1)
                v[pg_idx] = u[pg_idx] + (l[pg_idx] - mu1) / rho[pg_idx]
                v[pg_idx+1] = u[pg_idx+1] + (l[pg_idx+1] - mu2) / rho[pg_idx+1]
            end
            for j=FrStart[I]:FrStart[I+1]-1
                pij_idx = line_start + 8*(FrIdx[j]-1)
                v[pij_idx] = u[pij_idx] + (l[pij_idx] + mu1) / rho[pij_idx]
                v[pij_idx+1] = u[pij_idx+1] + (l[pij_idx+1] + mu2) / rho[pij_idx+1]
                v[pij_idx+4] = wi
                v[pij_idx+6] = ti
            end
            for j=ToStart[I]:ToStart[I+1]-1
                pij_idx = line_start + 8*(ToIdx[j]-1)
                v[pij_idx+2] = u[pij_idx+2] + (l[pij_idx+2] + mu1) / rho[pij_idx+2]
                v[pij_idx+3] = u[pij_idx+3] + (l[pij_idx+3] + mu2) / rho[pij_idx+3]
                v[pij_idx+5] = wi
                v[pij_idx+7] = ti
            end
        end
    end
end
