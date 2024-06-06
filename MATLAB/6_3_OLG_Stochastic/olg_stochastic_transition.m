
%==========================================================================
% STOCHASTIC OLG: Transition Dynamics
% Quantiative Macro Textbook Chapter 7
% WRITTEN By Sagiri Kitao
% Comments welcome --> sagiri.kitao@gmail.com
%==========================================================================

clc; close all; clear all;

ind_TR=3;
    % =1 COMPUTE INITIAL/FINAL SS
    % =2 COMPUTE TRANSITION STARTING WITH GENERIC INITIAL GUESS (AFTER 1)
    % =3 COMPUTE TRANSITION STARTING WITH SAVED INITIAL GUESS (AFTER 2)
    % =4 PLOT TRANSITION RESULTS (AFTER 2/3)
    % =5 COMPUTE WELFARE EFFECTS


%==========================================================================
% SET PARAMETER VALUES
%==========================================================================

alpha=0.40;         % CAPITAL SHARE
delta=0.08;         % DEPRECIATION RATE
beta=0.98;          % DISCOUNT FACTOR

rho0=0.5;           % SS REPLACEMENT RATE (INITIAL SS)
rho1=0.25;          % SS REPLACEMENT RATE (FINAL SS)

NT=100;             % TRANSITION PERIOD

% ITERATION PARAMETERS
maxiter=2000;       % MAX NUMBER OF ITERATIONS
tol=1e-3;           % ERROR TOLERANCE LEVEL
adjK=0.2;           % ADJUSTMENT FACTOR OF CAPITAL IN UPDATING GUESS


% GRID CONSTRUCTION
Nj=61;              % MAX AGE: 61 YEARS MAX (AGE 20-80)
NjW=45;             % WORKING YEARS (RETIRE AT NjW+1) : ENTER AT 21 AND RETIRE AT 65
Na=101;             % ASSET STATE
Na2=2001;           % ASSET CHOICE
Ne=2;               % LABOR PRODUCTIVITY (LOW/HIGH)


minA=0;             % MIN ASSET GRID
maxA=25;            % MAX ASSET GRID
curvA=1.2;          % ASSET GRID DENSITY (=1 EQUAL SIZE)



% LABOR PRODUCTIVITY GRID (MEAN=1)
gride=zeros(Ne,1);
dif=0.3;
gride(1)=1-dif;
gride(2)=1+dif;


% LABOR PRODUCTIVITY TRANSITION MATRIX
Pe=zeros(Ne,Ne);
Pe(1,1)=0.8;
Pe(1,2)=1-Pe(1,1);

Pe(2,2)=Pe(1,1);
Pe(2,1)=1-Pe(2,2);



% AGE DISTRIBUTION (ASSUME NO DEATH AND NO POP GROWTH. SUM TO 1)
meaJ=(1/Nj)*ones(Nj,1);


% AGGREGATE LABOR SUPPLY
L=sum(meaJ(1:NjW));


% ASSET STATE GRID
grida=zeros(Na,1);
grida(1)=minA;
for ac=2:Na
    grida(ac)=grida(1)+(maxA-minA)*((ac-1)/(Na-1))^curvA;
end


% ASSET CHOICE GRID
grida2=zeros(Na2,1);
grida2(1)=minA;
for ac=2:Na2
    grida2(ac)=grida2(1)+(maxA-minA)*((ac-1)/(Na2-1))^curvA;
end


% SPLIT GRIDS IN grida2 TO NEARBY TWO GRIDS IN grida
ac1vec=zeros(Na2,1);
ac2vec=zeros(Na2,1);

pra1vec=zeros(Na2,1);
pra2vec=zeros(Na2,1);

for ac=1:Na2

    xx=grida2(ac);

    if xx>=grida(Na)
        ac1vec(ac)=Na;
        ac2vec(ac)=Na;

        pra1vec(ac)=1;
        pra2vec(ac)=0;
    else
        ind=1;
        while xx>grida(ind+1)
            ind=ind+1;
            if ind+1>=Na
                break
            end
        end

        ac1vec(ac)=ind;

        if ind<Na
            ac2vec(ac)=ind+1;

            dA=(xx-grida(ind))/(grida(ind+1)-grida(ind));
            pra1vec(ac)=1-dA;
            pra2vec(ac)=dA;
        else
            ac2vec(ac)=ind;

            pra1vec(ac)=1;
            pra2vec(ac)=0;
        end

    end

end


% PATH OF SS REPLACEMENT RATE
rhoT=zeros(NT,1);
TT=25; % CONVERGE TO rho1 IN TT YEARS
for tc=1:TT
    rhoT(tc)=rho0+((rho1-rho0)/(TT-1))*(tc-1);
end
rhoT(TT+1:NT)=rho1;


% PATH OF TAX RATE
tauT=zeros(NT,1);
for tc=1:NT
    tauT(tc)=rhoT(tc)*sum(meaJ(NjW+1:Nj))/sum(meaJ(1:NjW));
end
tau0=tauT(1);
tau1=tauT(NT);


if ind_TR==1 % COMPUTE INITIAL AND FINAL SS

    for ind_SS=1:2

        disp(['NOW COMPUTING SS ',num2str(ind_SS)])

        %=======================
        % COMPUTE INITIAL SS
        %=======================

        % VALUE FUNCTION/SOLUTION (INITIALIZATION)
        vfun=zeros(Nj,Ne,Na);
        afunG=zeros(Nj,Ne,Na);

        afun=zeros(Nj,Ne,Na);

        % INCOME GRID (INITIALIZATION)
        y=zeros(Nj,Ne);

        % DISTRIBUTION (INITIALIZATION)
        mea=zeros(Nj,Ne,Na);

        % SS REPLACEMENT RATE AND TAX RATE
        if ind_SS==1
            rho=rho0;
            tau=tau0;
        else
            rho=rho1;
            tau=tau1;
        end

        % INITIAL GUESS OF K
        K=6;
        if ind_SS==1
            K=6.1676;
        else
            K=6.8425;
        end

        for iter=1:maxiter

            % COMPUTE r/w/ss
            r=alpha*(K/L)^(alpha-1)-delta;
            w=(1-alpha)*(K/L)^alpha;
            ss=rho*w;

            % NET INCOME BY AGE
            for jc=1:NjW
                for ec=1:Ne
                    yvec(jc,ec)=(1-tau)*w*gride(ec);
                end
            end
            yvec(NjW+1:Nj,1:Ne)=ss;


            % HOUSEHOLD PROBLEM: FROM LAST AGE Nj TO 1 (BACKWARDS)

            % (1) AGE Nj (LAST AGE)
            jc=Nj;
            for ec=1:Ne
                for ac=1:Na
                    c=yvec(jc,ec)+(1+r)*grida(ac);
                    vfun(jc,ec,ac)=log(c);
                    afunG(jc,ec,ac)=1; % NO SAVING

                    afun(jc,ec,ac)=grida2(1);
                end
            end


            % (2) AGE Nj-1:1
            for jc=Nj-1:-1:1

                for ec=1:Ne

                    y=yvec(jc,ec);

                    for ac=1:Na

                        vtemp  = -1000000*ones(Na2,1); % INITIALIZATION (STORE VALUES)
                        accmax = Na2;

                        for acc=1:Na2

                            c=y+(1+r)*grida(ac)-grida2(acc);

                            if c<=0
                                accmax=acc-1;
                                break
                            end

                            acc1 = ac1vec(acc);
                            acc2 = ac2vec(acc);

                            vpr=0;
                            for ecc=1:Ne
                                vpr=vpr+Pe(ec,ecc)*(pra1vec(acc)*vfun(jc+1,ecc,acc1)+pra2vec(acc)*vfun(jc+1,ecc,acc2));
                            end
                            vtemp(acc)=log(c)+beta*vpr;

                        end % acc

                        [val,index]     = max(vtemp(1:accmax));

                        vfun(jc,ec,ac)  = val;
                        afunG(jc,ec,ac) = index; % GRID FROM grida2
                        afun(jc,ec,ac)  = grida2(index);

                    end % ac

                end % ec

            end % jc



            % COMPUTE DISTRIBUTION mea
            mea=zeros(Nj,Ne,Na);        % INITIALIZATION
            mea(1,:,1)=meaJ(1)/Ne;      % ZERO ASSET AT AGE 1

            for jc=1:Nj-1
                for ec=1:Ne
                    for ac=1:Na

                        meaX=mea(jc,ec,ac);

                        acc=afunG(jc,ec,ac);

                        acc1 = ac1vec(acc);
                        acc2 = ac2vec(acc);

                        pra1 = pra1vec(acc);
                        pra2 = pra2vec(acc);

                        for ecc=1:Ne
                            mea(jc+1,ecc,acc1)=mea(jc+1,ecc,acc1)+Pe(ec,ecc)*pra1*meaX;
                            mea(jc+1,ecc,acc2)=mea(jc+1,ecc,acc2)+Pe(ec,ecc)*pra2*meaX;
                        end

                    end
                end
            end

            errm=abs(sum(sum(sum(mea)))-1);
            if errm>1e-4
                disp(['error in computation of distribution: ind_SS=',num2str(ind_SS),'errm=',num2str(errm),' ind_SS=',num2str(ind_SS)])
                break
            end

            mea_maxA=sum(sum(mea(:,:,Na)));
            if mea_maxA > 1e-3
                disp(['measure at max asset grid is large',num2str(mea_maxA)])
            end

            % COMPUTE ERROR IN K=A
            A=sum(sum(sum(afun.*mea)));
            errK = abs(K-A);

            if errK<tol
                disp(['K converged: iter=',num2str(iter)])
                break
            end

            if iter>maxiter
                disp(['WARN: iter>maxiter: iter=',num2str(iter),' errK=',num2str(errK)])
            end

            % UPDATE K
            K=K+adjK*(A-K);
            disp([num2str(iter),'   errK   ',num2str(errK)])

        end % iter

        if ind_SS==1
            K_SS0=K;
            mea_SS0=mea;
            vfun_SS0=vfun; % USED FOR WELFARE COMPUTATION ONLY
        else
            K_SS1=K;
            vfun_SS1=vfun;
        end

    end % ind_SS=1 OR 2

    save('SS_saver','K_SS0','K_SS1','vfun_SS1','mea_SS0')
    save('SS_welf','vfun_SS0','mea_SS0')
    disp('saved')


elseif ind_TR==2 | ind_TR==3

    %==========================
    % INITIAL GUESS FOR KT: KT0
    %==========================

    load('SS_saver') % READ 'K_SS0','K_SS1','vfun_SS1','mea_SS0'

    if  ind_TR==2 % USE GENERIC INITIAL GUESS

        KT0=K_SS1*ones(NT,1);
        NT0=30;
        intK=(K_SS1-K_SS0)/(NT0-1);

        for tc=1:NT0
            KT0(tc)=K_SS0+intK*(tc-1);
        end

    elseif ind_TR==3 % USE SAVED INITIAL GUESS

        load('TR_saver') % READ SAVED 'KT0'

    end % ind_TR


    ind_PLOT=1;
    if ind_PLOT==1
        time=1:NT;
        figure('Name','KT guess')
        plot(time,KT0)
    end


    %=============================
    % COMPUTE TRANSITION
    %=============================

    % NEW CAPITAL (INITIALIZATION)
    KT1=zeros(NT,1);

    % POLICY FUNCTION (INITIALIZATION)
    afunGT=zeros(NT,Nj,Ne,Na);

    % VALUE OF COHORT (USED FOR WELFARE COMPUTATION)
    vfun_TR=zeros(NT,Ne);
    vfun_TR0=zeros(Nj,Ne,Na);

    maxiterTR = 2; %50; % 10;
    iterTR    = 1;

    errK     = 1;
    errKTol  = 1e-4; %0.001;

    errKvec=zeros(maxiterTR,1); 

    adjK=0.05;

    while (errK>errKTol) & (iterTR<maxiterTR)

        %=====================================================
        % COMPUTE VALUE FUNCTION FROM t=NT to 1 (BACKWARDS)
        %=====================================================

        vfun0 = vfun_SS1; % VALUE IN THE FINAL SS (NEXT PERIOD AT TIME NT)

        for tc=NT:-1:1

            r=alpha*((KT0(tc)/L)^(alpha-1))-delta;
            w=(1-alpha)*(KT0(tc)/L)^alpha;
            rho=rhoT(tc);
            tau=tauT(tc);
            ss=rho*w;

            % NET INCOME BY AGE
            for jc=1:NjW
                for ec=1:Ne
                    yvec(jc,ec)=(1-tau)*w*gride(ec);
                end
            end
            yvec(NjW+1:Nj,1:Ne)=ss;


            %  INITIALIZATION
            vfun1  = zeros(Nj,Ne,Na); % NEW VALUE FUNCTION (AT TIME tc)
            afunG  = zeros(Nj,Ne,Na); % SOLUTION GRID
            afun   = zeros(Nj,Ne,Na); % SOLUTION LEVEL


            % HOUSEHOLD PROBLEM: FROM LAST AGE Nj TO 1 (BACKWARDS) AT time tc

            % (1) AGE Nj (LAST AGE)
            jc=Nj;
            for ec=1:Ne
                for ac=1:Na
                    c=yvec(jc,ec)+(1+r)*grida(ac);
                    vfun1(jc,ec,ac)=log(c); % NOTE: vfun1
                    afunG(jc,ec,ac)=1; % NO SAVING
                    afun(jc,ec,ac)=grida2(1);
                end
            end


            % (2) AGE Nj-1:1
            for jc=Nj-1:-1:1

                for ec=1:Ne

                    y=yvec(jc,ec);

                    for ac=1:Na

                        vtemp  = -1000000*ones(Na2,1); % INITIALIZATION (STORE VALUES)
                        accmax = Na2;

                        for acc=1:Na2

                            c=y+(1+r)*grida(ac)-grida2(acc);

                            if c<=0
                                accmax=acc-1;
                                break
                            end

                            acc1 = ac1vec(acc);
                            acc2 = ac2vec(acc);

                            vpr=0;
                            for ecc=1:Ne
                                vpr=vpr+Pe(ec,ecc)*(pra1vec(acc)*vfun0(jc+1,ecc,acc1)+pra2vec(acc)*vfun0(jc+1,ecc,acc2)); % NOTE: vfun0
                            end
                            vtemp(acc)=log(c)+beta*vpr;

                        end % acc

                        [val,index]     = max(vtemp(1:accmax));

                        vfun1(jc,ec,ac) = val; % NOTE: vfun1
                        afunG(jc,ec,ac) = index; % GRID FROM grida2
                        afun(jc,ec,ac)  = grida2(index);

                    end % ac

                end % ec

            end % jc

            % UPDATE vfun0 FOR NEXT PERIOD (tc-1)
            vfun0 = vfun1;

            % SAVE POLICY FUNCTION (SOLUTION GRID)
            afunGT(tc,:,:,:)=afunG;

            % SAVE ASSET (LEVEL)
            afunT(tc,:,:,:)=afun;

            disp(['iterTR=',num2str(iterTR),'  tc=',num2str(tc)])


            vfun_TR(tc,:)=vfun1(1,:,1);
            if tc==1
                vfun_TR0=vfun1;
            end

        end % tc



        %=====================================================
        % COMPUTE DISTRIBUTION meaT: FROM t=1 TO NT (FORWARD)
        %=====================================================

        meaT=zeros(NT,Nj,Ne,Na); % INITIALIZATION

        meaT(1,:,:,:)=mea_SS0; % DIST IN THE INITIAL SS

        mea0=mea_SS0;

        errm=abs(sum(sum(sum(mea0)))-1);

        if errm>1e-4
            disp(['error in computation of distribution mea1: errm=',num2str(errm),'mea0 before computing transition'])
            stop
            %break
        end

        for tc=1:NT-1

            afunG(:,:,:)=afunGT(tc,:,:,:);

            mea1=zeros(Nj,Ne,Na);           % INITIALIZATION
            mea1(1,:,1)=meaJ(1)/Ne;         % ZERO ASSET AT AGE 1

            for jc=1:Nj-1
                for ec=1:Ne
                    for ac=1:Na

                        meaX=mea0(jc,ec,ac); % NOTE meaX/mea0

                        acc=afunG(jc,ec,ac);

                        acc1 = ac1vec(acc);
                        acc2 = ac2vec(acc);

                        pra1 = pra1vec(acc);
                        pra2 = pra2vec(acc);

                        for ecc=1:Ne
                            mea1(jc+1,ecc,acc1)=mea1(jc+1,ecc,acc1)+Pe(ec,ecc)*pra1*meaX; % NOTE mea1/meaX
                            mea1(jc+1,ecc,acc2)=mea1(jc+1,ecc,acc2)+Pe(ec,ecc)*pra2*meaX; % NOTE mea1/meaX
                        end

                    end
                end
            end

            errm=abs(sum(sum(sum(mea1)))-1);
            if errm>1e-4
                disp(['error in computation of distribution mea1: errm=',num2str(errm),' tc=',num2str(tc)])
                stop
                %break
            end

            mea_maxA=sum(sum(mea1(:,:,Na)));
            if mea_maxA > 1e-3
                disp(['measure at max asset grid is large: mea_maxA=',num2str(mea_maxA),' tc=',num2str(tc)])
            end

            meaT(tc+1,:,:,:)=mea1;

            mea0=mea1;


        end % tc


        %========================================
        % COMPUTE KT1
        %========================================

        errKT=zeros(NT,1);

        KT1(1)=KT0(1);   % PREDETERMINED
        errKT(1)=0;

        for tc=1:NT-1

            afun(:,:,:)=afunT(tc,:,:,:); % SAVING FOR THE NEXT PERIOD
            mea0(:,:,:)=meaT(tc,:,:,:);

            KT1(tc+1)=sum(sum(sum(mea0.*afun))); % SAVING AT THE BEGINNING OF NEXT PERIOD

            errKT(tc+1)=abs(KT1(tc)-KT0(tc));

        end % tc

        errK=max(errKT);

        errKvec(iterTR)=errK;

        % UPDATE GUESS KT0
        if errK>errKTol

            % KT0(1) IS PREDETERMINED
            for tc=2:NT
                KT0(tc)=KT0(tc)+adjK*(KT1(tc)-KT0(tc));
            end

        else
            disp(['converged: iterTR= ',num2str(iterTR), ' errK= ',num2str(errK)])
        end

        disp(['iterTR= ',num2str(iterTR),'  errK=',num2str(errK)])

        iterTR=iterTR+1;

        % SAVE OUTCOME (DO THIS IN EACH ITERATION JUST IN CASE)
        save('TR_saver','KT0')
        save('TR_solution','K_SS0','K_SS1','KT0')
        save('TR_welf','vfun_TR','vfun_TR0')

    end % while

    if iterTR==maxiterTR
        disp(['iteration not converged: iterTR = ',num2str(iterTR),'errK = ',num2str(errK)])        
    end

elseif ind_TR==4

    load('TR_solution')  %'K_SS0','K_SS1','KT0'

    for tc=1:NT
        rT(tc)=alpha*((KT0(tc)/L)^(alpha-1))-delta;
        wT(tc)=(1-alpha)*(KT0(tc)/L)^alpha;
    end
    
    maxY=100;
    norm=K_SS0;
    gridT=1:NT;

    figure('Name','Capital')
    hold on
    plot(1,K_SS0/norm,'o r','LineWidth',2)
    plot(gridT,KT0/norm,'b','LineWidth',2)
    plot(maxY,K_SS1/norm,'o r','LineWidth',2)
    hold off
    xlabel('Time')
    xlim([1 maxY])
    xlabel('期間')
    box on
    grid on
    set(gca,'Fontsize',14)
    set(gca,'FontName','Times New Roman')
    set(gcf,'color','w')
    saveas(gcf,'fig_olg2_tr_k.eps','epsc2')


    figure('Name','Interest Rate')
    hold on
    plot(1,rT(1),'o r','LineWidth',2)
    plot(gridT,rT,'b','LineWidth',2)
    plot(maxY,rT(NT),'o r','LineWidth',2)
    hold off
    xlabel('Time')
    xlim([1 maxY])
    xlabel('期間')
    box on
    grid on
    set(gca,'Fontsize',14)
    set(gca,'FontName','Times New Roman')
    set(gcf,'color','w')
    saveas(gcf,'fig_olg2_tr_r.eps','epsc2')


elseif ind_TR==5

    load('SS_welf') % vfun_SS0(Nj,Ne,Na),mea_SS0(Nj,Ne,Na)
    load('TR_welf') % vfun_TR(NT,Ne),vfun_TR0(Nj,Ne,Na)
    

    betaJ=zeros(Nj,1);

    for jc=1:Nj
        temp=0;
        for ic=jc:Nj
            temp=temp+beta^(ic-jc);
        end
        betaJ(jc)=temp;
    end

    welf0=zeros(Nj,Ne,Na);
    for jc=1:Nj
        for ec=1:Ne
            for ac=1:Na
                welf0(jc,ec,ac)=exp( (vfun_TR0(jc,ec,ac)-vfun_SS0(jc,ec,ac)) / betaJ(jc))-1;
            end
        end
    end

    welf0_JE=zeros(Nj,Ne);
    for jc=1:Nj
        for ec=1:Ne
            welf0_JE(jc,ec)=sum(welf0(jc,ec,:).*mea_SS0(jc,ec,:))/sum(mea_SS0(jc,ec,:));
        end
    end

    welfTR=zeros(NT,Ne);
    for tc=1:NT
        for ec=1:Ne
            welfTR(tc,ec)=exp( (vfun_TR(tc,ec)-vfun_SS0(1,ec,1)) / betaJ(1) )-1;
        end
    end

    time=1:NT;
    figure('Name','Welf tr')
    hold on
    plot(time,welfTR(:,1),'--k','LineWidth',3)
    plot(time,welfTR(:,2),'-k','LineWidth',3)
    hold off
    legend('low','high','Location','SE')
    box on
    grid on
    xlim([1 60])
    xlabel('コホート生年')
    ylabel('CEV')
    set(gca,'Fontsize',16)
%    set(gca,'FontName','Times New Roman')
    set(gcf,'color','w')
    saveas(gcf,'fig_welf_tr.eps','epsc2')



    age=1:Nj;
    age=age+19;
    figure('Name','Welf initial')
    hold on
    plot(age,welf0_JE(:,1),'--k','LineWidth',3)
    plot(age,welf0_JE(:,2),'-k','LineWidth',3)
    hold off
    legend('low','high','Location','SE')
    box on 
    grid on
    xlabel('年齢')
    ylabel('CEV')
    set(gca,'Fontsize',16)
%    set(gca,'FontName','Times New Roman')
    set(gcf,'color','w')
    saveas(gcf,'fig_welf_0.eps','epsc2')


    jc=21; % 40 YRS OLD
    welf0_AE=zeros(Na,Ne);
    for ac=1:Na
        for ec=1:Ne
            welf0_AE(ac,ec)=welf0(jc,ec,ac);
        end
    end


    figure('Name','Welf initial by asset')
    hold on
    plot(grida,welf0_AE(:,1),'--k','LineWidth',3)
    plot(grida,welf0_AE(:,2),'-k','LineWidth',3)
    hold off
    legend('low','high','Location','SE')
    box on 
    grid on
    xlabel('Asset')
    ylabel('CEV')
    set(gca,'Fontsize',16)
%    set(gca,'FontName','Times New Roman')
    set(gcf,'color','w')

end
