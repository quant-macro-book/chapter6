function policyA(p, aX)

    a = zeros(p.Nj + 1)
    a[2:p.Nj] = aX

    @views K = dot(a[1:p.Nj], p.meaJ)

    r = p.α * (K / p.L)^(p.α - 1) - p.δ # INTEREST RATE
    w = (1 - p.α) * (K / p.L)^p.α  # WAGE

    ss = p.ss0 * w
    f = zeros(p.Nj - 1)
    # EULER EQUATION (u(c)=log(c))
    @inbounds for j in 1:p.Nj - 1
        ff = p.β * (1.0 + r) * ((1.0 + r) * a[j] + w * p.θ[j] * (1.0 - p.τ) + ss[j] - a[j + 1]) # beta*(1+r)*c, where c=(1+r)a+w*theta*(1-tau)-a'
        f[j] = (1.0 + r) * a[j + 1] + w * p.θ[j + 1] * (1.0 - p.τ) + ss[j + 1] - a[j + 2] - ff   # 0=c'-ff
    end
    return f
end
