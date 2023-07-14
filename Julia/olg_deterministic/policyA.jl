function policyA(m,aX)
    """
    --------------------------------------------------------------
    === Returns Euler Equation Error given HH's Asset Profile ===
    --------------------------------------------------------------
    ※THIS SHOULD BE ZERO AT EQUILIBRIUM
    <input>
    ・m: model structure 
    ・aX: HH's asset profile(from age 2 to Nj -> (Nj-2+1) element vector)
    <output>
    ・f: Euler equation errors from age 1 to Nj-1 -> (Nj-1+1) element vector
    """

    a = zeros(m.Nj+1);
    # initial asset
    a[1] = 0.0;
    a[m.Nj+1] = 0.0;

    a[2:m.Nj] .= copy(aX);

    K = sum(a[1:end-1] .* m.meaJ); 

    r = m.alpha * (K/m.L)^(m.alpha-1) - m.del; # interest rate
    w = (1-m.alpha)*(K/m.L)^m.alpha; # wage

    ss = m.ss0.*w;

    # euler equation (u(c)=log(c))
    f = zeros(m.Nj-1)
    for j in 1:m.Nj-1
        # beta*(1+r)*c
        ff =  m.beta * (1+r) * ((1+r)*a[j] + w*m.theta[j]*(1-m.tau) + ss[j] - a[j+1]);
        # 0=c'-ff
        f[j] = (1+r)*a[j+1] + w*theta[j+1]*(1-m.tau) + ss[j+1] - a[j+2] - ff;
    end

    return f
end