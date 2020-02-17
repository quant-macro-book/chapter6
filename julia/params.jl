struct Params{TF <: AbstractFloat,TI <: Integer,TV <: Vector}
    # PARAMETERS
    α::TF
    δ::TF
    β::TF
    Nj::TI    # LIVE 60 YEARS MAX (80YRS OLD)
    meaJ::TV # AGE DISTRIBUTION (no population growth)
    θ::TV
    L::TF
    ss0::TV
    τ::TF
end


"""
    Params(ind_E)
Set parameters depending on experiments.

# Arguments

- `ind_E::Integer`: experiment types

# Return

- 'p::Params': parameter struct

"""
function Params(;ind_E = 1)

    α = 0.40
    δ = 0.08
    β = 0.98
    Nj = 61     # LIVE 60 YEARS MAX (80YRS OLD)
    NjW = 45    # WORKING YEARS (RETIRE AT NjW+1) : ENTER AT 21 AND RETIRE AT 65
    ρ = 0.5   # SS REPLACEMENT RATE (0.0 OR 0.5)
    meaJ = (1 / Nj) * ones(Nj) # AGE DISTRIBUTION (no population growth)

    θ = zeros(Nj)
    ss0 = zeros(Nj)

    θ[1:NjW] .= 1
    L = sum(meaJ[1:NjW] .* θ[1:NjW])
    ss0[NjW + 1:Nj] .= ρ
    τ = ρ * (Nj - NjW) / NjW

    if ind_E == 1

        return Params(α, δ, β, Nj, meaJ, θ, L, ss0, τ)

    elseif ind_E == 2
        ρ = 0.25
        ss0[NjW + 1:Nj] .= ρ
        τ = ρ * (Nj - NjW) / NjW

        return Params(α, δ, β, Nj, meaJ, θ, L, ss0, τ)

    elseif ind_E == 3
        NjW = 50

        θ = zeros(Nj)
        ss0 = zeros(Nj)

        θ[1:NjW] .= 1
        L = sum(meaJ[1:NjW] .* θ[1:NjW])
        ss0[NjW + 1:Nj] .= ρ
        τ = ρ * (Nj - NjW) / NjW

        return Params(α, δ, β, Nj, meaJ, θ, L, ss0, τ)

    elseif ind_E == 4

        Nj = 66
        meaJ = (1 / 61) * ones(Nj) # IN EXPERIMENT TO RAISE LONGEVITY

        θ = zeros(Nj)
        ss0 = zeros(Nj)

        θ[1:NjW] .= 1
        L = sum(meaJ[1:NjW] .* θ[1:NjW])
        ss0[NjW + 1:Nj] .= ρ
        τ = ρ * (Nj - NjW) / NjW

        return Params(α, δ, β, Nj, meaJ, θ, L, ss0, τ)

    else
        throw(DomainError(ind_E, "ind_E must be between 1 and 4"))
    end
end
