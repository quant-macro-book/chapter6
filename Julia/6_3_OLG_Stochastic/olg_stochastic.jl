# ============================================================ #
#  STOCHASTIC OLG: STEADY STATE                                #
#  Quantiative Macro Textbook Chapter 7                        #
#  WRITTEN By Sagiri Kitao (Translated in Julia by Taiki Ono)  #
# ============================================================ #

using Plots

function main()
    
    # ================ #
    #  SET PARAMETERS  #
    # ================ #
    
    alpha = 0.40;   # capital share
    delta = 0.08;   # depreciation rate
    beta = 0.98;    # discount factor
    rho = 0.5;      # SS replacement rate

    # iteration parameters
    maxiter = 2000; # max number of iterations
    tol = 1e-3;     # error tolerance level
    adjK = 0.2;     # adjustment factor of capital in updating guess

    # grid construction
    Nj = 61;        # max age: 61 years max (age 20-80)
    Njw = 45;       # working years (retire at Njw+1: enter at 21 and retire at 65)
    Na = 201;       # asset state
    Na2 = 8001;     # asset choice
    Ne = 2;         # labor productivity (low/high)

    minA = 0;       # min asset grid
    maxA = 25;      # max asset grid
    curvA = 1.2;    # asset grid density (=1 eequal size)

    # labor productivity grid (mean=1)
    gride = zeros(Ne);
    dif = 0.3;
    gride[1] = 1.0 - dif;
    gride[2] = 1.0 + dif;

    # labor productivity transition matrix
    Pe = zeros(Ne,Ne);
    Pe[1,1] = 0.8;
    Pe[1,2] = 1.0 - Pe[1,1];

    Pe[2,2] = copy(Pe[1,1]);
    Pe[2,1] = 1.0 - Pe[2,2];


    # age distribution (assume no death and no pop growth. sum to 1)
    meaJ = 1/Nj .* ones(Nj);


    # aggregate labor supply
    L = sum(meaJ[1:Njw]);


    # asset state grid
    grida = zeros(Na);
    grida[1] = minA;
    for ac in 2:Na
        grida[ac] = grida[1] + (maxA-minA)*((ac-1)/(Na-1))^curvA;
    end

    # asset choice grid
    grida2 = zeros(Na2);
    grida2[1] = minA;
    for ac in 2:Na2
        grida2[ac] = grida2[1] + (maxA-minA)*((ac-1)/(Na2-1))^curvA;
    end


    # split grids in grida2 to nearby two grids in grida
    ac1vec = zeros(Int64,Na2);
    ac2vec = zeros(Int64,Na2);

    pra1vec = zeros(Na2);
    pra2vec = zeros(Na2);

    for ac in 1:Na2

        xx = grida2[ac];

        if xx >= grida[Na]
  
            ac1vec[ac] = Na;
            ac2vec[ac] = Na;

            pra1vec[ac] = 1.0;
            pra2vec[ac] = 0.0;

        else

            ind = 1;
            while xx > grida[ind+1]
                ind += 1;
                if ind+1 >= Na
                    break
                end
            end

            ac1vec[ac] = ind;

            if ind < Na

                ac2vec[ac] = ind+1;
                dA = (xx-grida[ind])/(grida[ind+1]-grida[ind]);
                pra1vec[ac] = 1.0-dA;
                pra2vec[ac] = dA;

            else

                ac2vec[ac] = ind;
                pra1vec[ac] = 1.0;
                pra2vec[ac] = 0.0;

            end
        end
    end



    # value function/solution (initialization)
    vfun = zeros(Nj,Ne,Na);
    afunG = zeros(Int64,Nj,Ne,Na);
    afun = zeros(Nj,Ne,Na);

    # income grid (initialization)
    yvec = zeros(Nj,Ne);

    # distribution (initialization)
    mea = zeros(Nj,Ne,Na)

    # equilibrium tax rate
    tau = rho*sum(meaJ[Njw+1:Nj])/sum(meaJ[1:Njw]);
    
    # initial guess of K
    K = 7.0;
    r = 0.0;

    for iter = 1:maxiter

        # compute r/w/SS
        r = alpha*(K/L)^(alpha-1) - delta;
        w = (1-alpha)*(K/L)^alpha;
        ss = rho*w
        
        # net income by age
        for jc in 1:Njw
            for ec in 1:Ne
                yvec[jc,ec] = (1-tau)*w*gride[ec]
            end
        end
        yvec[Njw+1:Nj,:] .= ss;


        # household problem: from last age Nj to 1 (backwards)
        
        # (1) age Nj (last age)
        jc = Nj;
        for ec in 1:Ne
            for ac in 1:Na
                
                c = yvec[jc,ec] + (1+r)*grida[ac];
                vfun[jc,ec,ac] = log(c);
                afunG[jc,ec,ac] = 1; # no saving
                afun[jc,ec,ac] = grida2[1];

            end
        end

        # (2) age Nj-1:1
        for jc in Nj-1:-1:1
            for ec in 1:Ne

                y = yvec[jc,ec];

                for ac in 1:Na

                    vtemp = -1000000 .* ones(Na2); # initialization (store values)
                    accmax = Na2;

                    for acc in 1:Na2

                        c = y + (1+r)*grida[ac] - grida2[acc];

                        if c <= 0.0
                            accmax = acc-1;
                            break
                        end

                        acc1 = ac1vec[acc];
                        acc2 = ac2vec[acc];

                        vpr = 0.0;
                        for ecc in 1:Ne
                            vpr += Pe[ec,ecc]*(pra1vec[acc]*vfun[jc+1,ecc,acc1] + pra2vec[acc]*vfun[jc+1,ecc,acc2]);
                        end

                        vtemp[acc] = log(c) + beta*vpr;

                    end

                    val,index = findmax(vtemp[1:accmax]);
                    vfun[jc,ec,ac] = val;
                    afunG[jc,ec,ac] = index; # grid from grida2
                    afun[jc,ec,ac] = grida2[index];

                end
            end
        end

        # compute distribution mea
        mea = zeros(Nj,Ne,Na);    # initialization
        mea[1,:,1] .= meaJ[1]/Ne; # zero asset at age 1

        for jc in 1:Nj-1
            for ec in 1:Ne
                for ac in 1:Na

                    mea0 = mea[jc,ec,ac];
                    acc = afunG[jc,ec,ac];
                    
                    acc1 = ac1vec[acc];
                    acc2 = ac2vec[acc];

                    pra1 = pra1vec[acc];
                    pra2 = pra2vec[acc];

                    for ecc in 1:Ne
                        
                        mea[jc+1,ecc,acc1] += Pe[ec,ecc]*pra1*mea0;
                        mea[jc+1,ecc,acc2] += Pe[ec,ecc]*pra2*mea0;

                    end
                end
            end
        end

        errm = abs(sum(mea)-1);
        if errm > 1e-4
            println("error in computation of distribution: $errm")
            break
        end

        mea_maxA = sum(mea[:,:,Na]);
        if mea_maxA > 1e-3
            println("measure at max asset grid is large: $mea_maxA")
        end

        # compute error in K=A
        A = sum(afun.*mea);
        errK = abs(K-A);

        if errK < tol
            println("K coverged: $iter")
            break
        end

        if iter > maxiter
            println("WARN: iter>maxiter: $iter, $errK")
        end

        # update Kitao
        K += adjK*(A-K);

        println("$iter errK = $errK")

    end


    # =============== #
    #  COMPUTE STATS  #
    # =============== #

    # asset by age
    afunJ = zeros(Nj);
    for jc in 1:Nj
        tempA = 0.0;
        for ac in 1:Na
            tempA += grida[ac] * sum(mea[jc,:,ac])
        end
        afunJ[jc] = tempA / meaJ[jc]
    end 

    # asset by age/skill
    afunJE = zeros(Nj,Ne);
    for jc in 1:Nj
        for ec in 1:Ne
            temp = 0.0;
            for ac in 1:Na
                temp += grida[ac] * mea[jc,ec,ac];
            end
            afunJE[jc,ec] = temp/sum(mea[jc,ec,:]);
        end
    end

    cfun = zeros(Nj,Ne,Na); # consumption function
    sfun = zeros(Nj,Ne,Na); # saving function
    srat = zeros(Nj,Ne,Na); # stores saving rate given state

    for jc in 1:Nj
        for ec in 1:Ne
            for ac in 1:Na

                acc = afunG[jc,ec,ac];
                y = yvec[jc,ec];
                c = y + (1+r)*grida[ac] - grida2[acc];
                inc = y + r*grida[ac];

                sfun[jc,ec,ac] = inc - c;
                cfun[jc,ec,ac] = c;
                srat[jc,ec,ac] = (inc-c)/inc;

            end
        end
    end

    # consumption by age
    cfunJ = zeros(Nj);
    for jc in 1:Nj
        cfunJ[jc] = sum(mea[jc,:,:] .* cfun[jc,:,:]) / sum(mea[jc,:,:]);
    end

    # saving and saving rate by age/skill
    sfunJE = zeros(Nj,Ne);
    sratJE = zeros(Nj,Ne);
    for jc in 1:Nj
        for ec in 1:Ne

            sfunJE[jc,ec] = sum(mea[jc,ec,:] .* sfun[jc,ec,:]) / sum(mea[jc,ec,:]);
            sratJE[jc,ec] = sum(mea[jc,ec,:] .* srat[jc,ec,:]) / sum(mea[jc,ec,:]);

        end
    end

    age = collect(1:Nj);
    age .+= 19;

    norm = 1/cfunJ[1];

    minJ = 20;
    maxJ = 80;


    plot(age,norm.*afunJE[:,2],ls=:solid,lc=:black,lw=3,label="high",legend=:topleft)
    plot!(age,norm.*afunJE[:,1],ls=:dashdot,lc=:black,lw=3,label="low")
    xlims!(minJ,maxJ)
    title!("Asset (Stock) By Age and Prod")
    xlabel!("Age")
    ylabel!("Asset")
    savefig("fig_olg2_a.pdf")

    
    # saving rate by asset: choose a particular age to plot
    jc = 21;
    sfunAE = zeros(Na,Ne);
    sratAE = zeros(Na,Ne);
    for ac in 1:Na
        for ec in 1:Ne

            sfunAE[ac,ec] = sfun[jc,ec,ac];
            sratAE[ac,ec] = srat[jc,ec,ac];

        end
    end

    minA = 0.0;
    maxA = 20.0;
    minY = -0.2;
    maxY = 0.5;
    zerovec = zeros(Na);

    plot(grida,sratAE[:,2],ls=:solid,lc=:black,lw=3,label="high")
    plot!(grida,sratAE[:,1],ls=:dashdot,lc=:black,lw=3,label="low")
    plot!(grida,zerovec,ls=:dot,lc=:black,lw=1,label="")
    xlims!(minA,maxA)
    title!("Saving Rate By Asset and Prod")
    xlabel!("Asset")
    ylabel!("Saving Rate")
    savefig("fig_olg2_s.pdf")


    plot(grida,sfunAE[:,2],ls=:solid,lc=:black,lw=3,label="high")
    plot!(grida,sfunAE[:,1],ls=:dashdot,lc=:black,lw=3,label="low")
    plot!(grida,zerovec,ls=:dot,lc=:black,lw=1,label="")
    xlims!(minA,maxA)
    title!("Saving Level By Asset and Prod")
    xlabel!("Asset")
    ylabel!("Saving (level)")
    savefig("fig_olg2_slevel.pdf")

    println(sratAE[:,1])

end

main()

