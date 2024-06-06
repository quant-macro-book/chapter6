# =========================================================== #
#  STOCHASTIC OLG: Transition Dynamics                        #
#  Quantitative Macro Textbook Chapter7                       #
#  WRITTEN By Sagiri Kitao (Translated in Julia by Taiki Ono) #
# =========================================================== #

# import libraries
using Plots
using JLD


function main()
    
    ind_TR = 4;
        # =1 COMPUTE INITIAL/FINAL SS + TRANSITION STARTING WITH GENERIC INITIAL GUESS
        # =2 COMPUTE INITIAL/FINAL SS + TRANSITION STARTING WITH SAVED GUESS
        # =3 PLOT TRANSITION RESULTS (AFTER 1/2)
        # =4 COMPUTE WELFARE EFFECTS (AFTER 3)


    # ====================== #
    #  SET PARAMETER VALUES  #
    # ====================== #

    alpha = 0.40;   # capital share
    delta = 0.08;   # depreciation rate
    beta = 0.98;    # discount factor

    rho0 = 0.5;     # SS replacement rate (initial SS)
    rho1 = 0.25;    # SS replacement rate (final SS) 

    NT = 100;       # transition period

    # iteration parameters
    maxiter = 2000; # max number of iterations
    tol = 1e-3;     # error tolerance level
    adjK = 0.2;     # adjustment factor of  capital in updating guess

    # grid construction
    Nj = 61;        # max age: 61 years max (age 20-80)
    Njw = 45;       # working years (retires at Njw+1): enter at 21 and retire at 65
    Na = 101;       # asset state
    Na2 = 2001;     # asset choice
    Ne = 2;         # labor productivity (low/high)

    minA = 0;       # min asset grid
    maxA  = 25;     # max asset grid
    curvA = 1.2;    # asset grid density (=1 equal size)

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

    # path of SS replacement rate
    rhoT = zeros(NT);
    TT = 25; # coverged to rho1 in TT years
    for tc in 1:TT
        rhoT[tc] = rho0 + ((rho1-rho0)/(TT-1))*(tc-1);
    end
    rhoT[TT+1:NT] .= rho1;

    # path of tax rate
    tauT = zeros(NT);
    for tc in 1:NT
        tauT[tc] = rhoT[tc]*sum(meaJ[Njw+1:Nj])/sum(meaJ[1:Njw]);
    end
    tau0 = copy(tauT[1]);
    tau1 = copy(tauT[NT]);

    
    if (ind_TR == 1) || (ind_TR == 2) # compute initial and final SS

        # ======================================================================== #
        # THIS IS TO AVOID USING GLOBAL VARIABLES FOR FASTER COMPUTATION IN JULIA

        K_SS0 = 0.0;
        mea_SS0 = zeros(Nj,Ne,Na);
        vfun_SS0 = zeros(Nj,Ne,Na);
        K_SS1 = 0.0;
        vfun_SS1 = zeros(Nj,Ne,Na);
        afun = zeros(Nj,Ne,Na);


        # ======================================================================== #

        for ind_SS in 1:2

            println("NOW COMPUTING SS $ind_SS")

            # ===================== #
            #  COMPUTING INTIAL SS  #
            # ===================== #

            # value function/solution (initialization)
            vfun = zeros(Nj,Ne,Na);
            afunG = zeros(Int64,Nj,Ne,Na);
            afun = zeros(Nj,Ne,Na);

            # income grid (initialization)
            yvec = zeros(Nj,Ne);

            # distribution (initialization)
            mea = zeros(Nj,Ne,Na);

            # SS replacement rate and tax rate
            if ind_SS == 1
                rho = copy(rho0);
                tau = copy(tau0);
            else
                rho = copy(rho1);
                tau = copy(tau1);
            end

            # initial guess of K
            K = 6.0;
            if ind_SS == 1
                K = 6.1676;
            else
                K = 6.8425;
            end

            for iter in 1:maxiter

                # compute r/w/SS
                r = alpha*(K/L)^(alpha-1) - delta;
                w = (1-alpha)*(K/L)^alpha;
                ss = rho*w;

                # net income by age
                for jc in 1:Njw
                    for ec in 1:Ne
                        yvec[jc,ec] = (1-tau)*w*gride[ec];
                    end
                end
                yvec[Njw+1:Nj,:] .= ss;

                
                # hosehold problem: from last age Nj to 1 (backwards)
                
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

                                acc1 = ac1vec[Int(acc)];
                                acc2 = ac2vec[Int(acc)];

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
                    
                            acc1 = ac1vec[Int(acc)];
                            acc2 = ac2vec[Int(acc)];

                            pra1 = pra1vec[Int(acc)];
                            pra2 = pra2vec[Int(acc)];

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
                    flush(stdout)
                    break
                end

                if iter > maxiter
                    println("WARN: iter>maxiter: $iter, $errK")
                end

                # update Kitao
                K += adjK*(A-K);
                println("$iter errK = $errK")
                
            end

            if ind_SS == 1

                K_SS0 = copy(K);
                mea_SS0 = copy(mea);
                vfun_SS0 = copy(vfun); # used for welfare computation only

            else

                K_SS1 = copy(K);
                vfun_SS1 = copy(vfun);

            end
        end

        save("SS_saver.jld","K_SS0",K_SS0,"K_SS1",K_SS1,"vfun_SS1",vfun_SS1,"mea_SS0",mea_SS0);
        save("SS_welf.jld","vfun_SS0",vfun_SS0,"mea_SS0",mea_SS0);
        
        
        println("saved")
        flush(stdout)


        # ========================= #
        #  INITIAL GUESS OF KT:KT0  #
        # ========================= #

        if ind_TR == 1
            # use generic initial guess
            KT0 = K_SS1 .* ones(NT);
            NT0 = 30;
            intK = (K_SS1-K_SS0)/(NT0-1);

            for tc in 1:NT0
                KT0[tc] = K_SS0 + intK*(tc-1)
            end

        elseif ind_TR == 2
            # use saved initial guess
            KT0 = copy(load("TR_saver.jld")["KT0"]); # read saved "KT0"
        end

        ind_PLOT = 1;
        if ind_PLOT == 1
            
            time = collect(1:NT);
            plot(time,KT0,labels="")
            title!("KT guess")
        
        end


        # ==================== #
        #  COMPUTE TRANSITION  #
        # ==================== #

        # new capital (initialization)
        KT1 = zeros(NT);

        # policy function (initialization)
        afunGT = zeros(Int64,NT,Nj,Ne,Na);

        # asset (level): (ADDED IN JULIA CODE)
        afunT = zeros(NT,Nj,Ne,Na);

        # value of cohort (used for welfare computation)
        vfun_TR = zeros(NT,Ne);
        vfun_TR0 = zeros(Nj,Ne,Na);

        maxiterTR = 300; 
        iterTR = 1;

        errK = 1.0;
        errKTol = 1e-4; # 0.001

        errKvec = zeros(maxiterTR)

        adjK = 0.05;

        while (errK > errKTol) && (iterTR < maxiterTR)

            # =================================================== #
            #  COMPUTE VALUE FUNCTION FROM t=NT to 1 (BACKWARDS)  #
            # =================================================== #

            vfun0 = copy(vfun_SS1) # value in the final SS (next period at time NT)

            for tc in NT:-1:1

                r = alpha * ((KT0[tc]/L)^(alpha-1)) - delta;
                w = (1-alpha) * (KT0[tc]/L)^alpha;
                rho = rhoT[tc];
                tau = tauT[tc];
                ss = rho*w;

                # net income by age
                yvec = zeros(Nj,Ne) # initialization (ADDED IN JULIA CODE)
                for jc in 1:Njw
                    for ec in 1:Ne
                        yvec[jc,ec] = (1-tau)*w*gride[ec];
                    end
                end
                yvec[Njw+1:Nj,:] .= ss;

                # initialization
                vfun1 = zeros(Nj,Ne,Na); # new value function (at time tc)
                afunG = zeros(Int64,Nj,Ne,Na); # solution grid
                afun = zeros(Nj,Ne,Na);  # solution level
                
                # hosehold problem from last age Nj to 1 (backwards) at time tc
                # age Nj (last age)
                jc = Nj;
                for ec in 1:Ne
                    for ac in 1:Na

                        c = yvec[jc,ec] + (1+r)*grida[ac];
                        vfun1[jc,ec,ac] = log(c); # note: vfun_SS1
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
                                    vpr += Pe[ec,ecc]*(pra1vec[acc]*vfun0[jc+1,ecc,acc1]+pra2vec[acc]*vfun0[jc+1,ecc,acc2]) # note: vfun0
                                end
                                vtemp[acc] = log(c) + beta*vpr;

                            end

                            val,index = findmax(vtemp[1:accmax]);
                            vfun1[jc,ec,ac] = val; # note: vfun1
                            afunG[jc,ec,ac] = index; # grid from grida2
                            afun[jc,ec,ac] = grida2[index];

                        end
                    end
                end

                # update vfun0 for next periiod (tc-1)
                vfun0 = copy(vfun1);

                # save policy function (solution grid)
                afunGT[tc,:,:,:] .= copy(afunG);

                # save asset (level)
                afunT[tc,:,:,:] .= copy(afun);
                
                println("iterTR = $iterTR, tc = $tc");
                flush(stdout)

                vfun_TR[tc,:] = copy(vfun1[1,:,1]);
                if tc == 1
                    vfun_TR0 = copy(vfun1);
                end

            end


            # ==================================================== #
            #  COMPUTE DISTRIBUTION meaT: FROM t=1 TO NT (FORWARD) #
            # ==================================================== #

            meaT = zeros(NT,Nj,Ne,Na) # initialization
            mea0 = zeros(Nj,Ne,Na);   # initialization
        

            meaT[1,:,:,:] .= copy(mea_SS0) # dist in the initial SS

            mea0 .= copy(mea_SS0);

            errm = sum(mea0) - 1;

            if errm > 1e-4
                error("error in computation of distribution mea1: errm = $errm, mea0 before computing transition")
            end

            for tc in 1:NT-1

                afunG = copy(afunGT[tc,:,:,:]);

                mea1 = zeros(Nj,Ne,Na);    # initialization
                mea1[1,:,1] .= meaJ[1]/Ne; # zero asset at age 1

                for jc in 1:Nj-1
                    for ec in 1:Ne
                        for ac in 1:Na

                            meaX = mea0[jc,ec,ac] # note meaX/mea0
                            
                            acc = afunG[jc,ec,ac];

                            acc1 = ac1vec[acc];
                            acc2 = ac2vec[acc];

                            pra1 = pra1vec[acc];
                            pra2 = pra2vec[acc];

                            for ecc in 1:Ne

                                mea1[jc+1,ecc,acc1] += Pe[ec,ecc]*pra1*meaX # note mea1/meaX
                                mea1[jc+1,ecc,acc2] += Pe[ec,ecc]*pra2*meaX # note mea1/meaX

                            end
                        end
                    end
                end

                errm = abs(sum(mea1)-1);
                if errm > 1e-4
                    error("error in computation of distribution mea1: errm = $errm, tc = $tc")
                end

                mea_maxA = sum(mea1[:,:,Na]);
                if mea_maxA > 1e-3
                    println("measure at max asset grid is large: mea_maxA = $mea_maxA, tc = $tc")
                    flush(stdout)
                end

                meaT[tc+1,:,:,:] .= copy(mea1);
                
                mea0 = copy(mea1)

            end


            # ============= #
            #  COMPUTE KT1  #
            # ============= # 

            errKT = zeros(NT);
            KT1[1] = copy.(KT0[1]); # predetermined
            errKT[1] = 0.0;
            

            for tc in 1:NT-1

                afun .= copy(afunT[tc,:,:,:]); # saving for the next period
                mea0 .= copy(meaT[tc,:,:,:]);

                KT1[tc+1] = sum(mea0.*afun); # saving at the begginig of next period

                errKT[tc+1] = abs(KT1[tc]-KT0[tc]);

            end

            errK = maximum(errKT);
            errKvec[iterTR] = errK;

            # update guess KT0
            if errK > errKTol

                # KT0[1] is predetermined
                for tc in 2:NT
                    KT0[tc] += adjK*(KT1[tc]-KT0[tc]);
                end

            else
                println("converged: iterTR = $iterTR, errK = $errK")
                flush(stdout)
            end

            println("iterTR = $iterTR, errK = $errK")
            flush(stdout)

            iterTR += 1;

            # save outcome (DO THIS IN EACH ITERATION JUST IN CASE)
            save("TR_saver.jld","KT0",KT0);
            save("TR_solution.jld","K_SS0",K_SS0,"K_SS1",K_SS1,"KT0",KT0);
            save("TR_welf.jld","vfun_TR",vfun_TR,"vfun_TR0",vfun_TR0);

        end

        if iterTR == maxiterTR
            println("iteration not converged: iterTR = $iterTR, errK = $errK")
            flush(stdout)
        end

    elseif ind_TR == 3
        
        # load saved arrays in ind_TR = 1
        TR_solution = load("TR_solution.jld");
        K_SS0 = copy(TR_solution["K_SS0"]);
        K_SS1 = copy(TR_solution["K_SS1"]);
        KT0 = copy(TR_solution["KT0"]);

        # plot transition dynamics
        rT = zeros(NT);
        wT = zeros(NT);

        for tc in 1:NT

            rT[tc] = alpha * ((KT0[tc]/L)^(alpha-1)) - delta;
            wT[tc] = (1-alpha) * (KT0[tc]/L)^alpha;

        end

        maxY = 100;
        norm = copy(K_SS0);
        gridT = collect(1:NT);

        plot([1],[K_SS0/norm], mc=:red, markershapes=:circle, lw=2, legend=false)
        plot!(gridT, KT0./norm, color=:blue, lw=2)
        plot!([maxY], [K_SS1/norm], mc=:red, markershapes=:circle, lw=2)
        title!("Capital")
        xlabel!("Time")
        xlims!(1-0.9,maxY+0.9)
        savefig("fig_olg2_tr_k.pdf")


        plot([1],[rT[1]], mc=:red, markershapes=:circle, lw=2, legend=false)
        plot!(gridT, rT, color=:blue, lw=2)
        plot!([maxY], [rT[NT]], mc=:red, markershapes=:circle, lw=2)
        title!("Interest Rate")
        xlabel!("Time")
        xlims!(1-0.9,maxY+0.9)
        savefig("fig_olg2_tr_r.pdf")

    elseif ind_TR == 4

        SS_welf = load("SS_welf.jld"); # vfun_SS0(Nj,Ne,Na),mea_SS0(Nj,Ne,Na)
        TR_welf = load("TR_welf.jld"); # vfun_TR(NT,Ne),vfun_TR0(Nj,Ne,Na)

        vfun_SS0 = copy(SS_welf["vfun_SS0"]);
        mea_SS0 = copy(SS_welf["mea_SS0"]);
        vfun_TR = copy(TR_welf["vfun_TR"]);
        vfun_TR0 = copy(TR_welf["vfun_TR0"]);

        betaJ = zeros(Nj);
        
        for jc in 1:Nj
            temp = 0.0;
            for ic in jc:Nj
                temp += beta^(ic-jc);
            end
            betaJ[jc] = temp;
        end 

        welf0 = zeros(Nj,Ne,Na);
        for jc in 1:Nj
            for ec in 1:Ne
                for ac in 1:Na
                    welf0[jc,ec,ac] = exp( (vfun_TR0[jc,ec,ac] - vfun_SS0[jc,ec,ac]) / betaJ[jc] ) - 1;
                end
            end
        end

        welf0_JE = zeros(Nj,Ne);
        for jc in 1:Nj
            for ec in 1:Ne
                welf0_JE[jc,ec] = sum(welf0[jc,ec,:].*mea_SS0[jc,ec,:])/sum(mea_SS0[jc,ec,:]);
            end
        end

        welfTR = zeros(NT,Ne);
        for tc in 1:NT
            for ec in 1:Ne
                welfTR[tc,ec] = exp( (vfun_TR[tc,ec]-vfun_SS0[1,ec,1]) / betaJ[1] ) - 1;
            end
        end 

        time = collect(1:NT);

        plot(time, welfTR[:,1], ls=:dash, c=:black, lw=3, label="low", legend=:bottomright)
        plot!(time, welfTR[:,2], ls=:solid, c=:black, lw=3, label="high")
        xlims!(1,60)
        ylims!(0.01,0.08)
        xlabel!("Cohort")
        ylabel!("CEV")
        savefig("fig_welf_tr.pdf")


        age = collect(1:Nj);
        age .+= 19;

        plot(age, welf0_JE[:,1], ls=:dash, c=:black, lw=3, label="low", legend=:bottomright)
        plot!(age, welf0_JE[:,2], ls=:solid, c=:black, lw=3, label="high")
        ylims!(-0.1,0.04)
        xlabel!("Age")
        ylabel!("CEV")
        savefig("fig_welf_0.pdf")


        jc = 21; # 40 yrs old
        welf0_AE = zeros(Na,Ne);
        for ac in 1:Na
            for ec in 1:Ne
                welf0_AE[ac,ec] = welf0[jc,ec,ac];
            end
        end

        plot(grida, welf0_AE[:,1], ls=:dash, c=:black, lw=3, label="low", legend=:bottomright)
        plot!(grida, welf0_AE[:,2], ls=:solid, c=:black, lw=3, label="high")
        xlabel!("Asset")
        ylabel!("CEV")
        savefig("fig_welf_aged_40.pdf")

    end
end

main()