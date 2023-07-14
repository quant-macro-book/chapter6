clear all; clc; close all;

%==========================================================================
% STOCHASTIC OLG: STEADY STATE
% Quantiative Macro Textbook Chapter 7
% WRITTEN By Sagiri Kitao
% Comments welcome --> sagiri.kitao@gmail.com
%==========================================================================

%==========================================================================
% SET PARAMETER VALUES
%==========================================================================

alpha=0.40;         % CAPITAL SHARE
delta=0.08;         % DEPRECIATION RATE
beta=0.98;          % DISCOUNT FACTOR
rho=0.5;            % SS REPLACEMENT RATE 


% ITERATION PARAMETERS
maxiter=2000;       % MAX NUMBER OF ITERATIONS
tol=1e-3;           % ERROR TOLERANCE LEVEL
adjK=0.2;           % ADJUSTMENT FACTOR OF CAPITAL IN UPDATING GUESS


% GRID CONSTRUCTION
Nj=61;              % MAX AGE: 61 YEARS MAX (AGE 20-80)
NjW=45;             % WORKING YEARS (RETIRE AT NjW+1) : ENTER AT 21 AND RETIRE AT 65
Na=201;             % ASSET STATE
Na2=8001;           % ASSET CHOICE
Ne=2;               % LABOR PRODUCTIVITY (LOW/HIGH)


minA=0;             % MIN ASSET GRID 
maxA=25;            % MAX ASSET GRID
curvA=1.2;          % ASSET GRID DENSITY (=1 EQUAL SIZE)


    
% LABOR PRODUCTIVITY GRID (MEAN=1)
gride=zeros(Ne,1);
dif=0.2; %0.3;
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


% VALUE FUNCTION/SOLUTION (INITIALIZATION) 
vfun=zeros(Nj,Na);
afunG=zeros(Nj,Ne,Na);
afun=zeros(Nj,Ne,Na);

% INCOME GRID (INITIALIZATION)
y=zeros(Nj,Ne);

% DISTRIBUTION (INITIALIZATION)
mea=zeros(Nj,Ne,Na);

% EQUILIBRIUM TAX RATE 
tau=rho*sum(meaJ(NjW+1:Nj))/sum(meaJ(1:NjW));

% INITIAL GUESS OF K
K=7;

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

                mea0=mea(jc,ec,ac);

                acc=afunG(jc,ec,ac);

                acc1 = ac1vec(acc);
                acc2 = ac2vec(acc);

                pra1 = pra1vec(acc);
                pra2 = pra2vec(acc);

                for ecc=1:Ne
                    mea(jc+1,ecc,acc1)=mea(jc+1,ecc,acc1)+Pe(ec,ecc)*pra1*mea0;
                    mea(jc+1,ecc,acc2)=mea(jc+1,ecc,acc2)+Pe(ec,ecc)*pra2*mea0;
                end

            end
        end
    end

    errm=abs(sum(sum(sum(mea)))-1);
    if errm>1e-4
        disp(['error in computation of distribution',num2str(errm)])
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
        disp(['K converged',num2str(iter)])
        iter
        break
    end

    if iter>maxiter
        disp(['WARN: iter>maxiter',num2str(iter),num2str(errK)])
    end

    % UPDATE K
    K=K+adjK*(A-K);

    disp([num2str(iter),'   errK   ',num2str(errK)])

end % iter


%==========================================================================
% COMPUTE STATS
%==========================================================================

% ASSET BY AGE
afunJ=zeros(Nj,1);
for jc=1:Nj    
    tempA=0;
    for ac=1:Na
        tempA=tempA+grida(ac)*sum(mea(jc,:,ac));
    end
    afunJ(jc)=tempA/meaJ(jc);
end


% ASSET BY AGE/SKILL
afunJE=zeros(Nj,Ne);
for jc=1:Nj    
    for ec=1:Ne
        temp=0;
        for ac=1:Na
            temp=temp+grida(ac)*mea(jc,ec,ac);
        end
        afunJE(jc,ec)=temp/sum(mea(jc,ec,:));
    end
end

cfun=zeros(Nj,Ne,Na);
sfun=zeros(Nj,Ne,Na);
srat=zeros(Nj,Ne,Na);
for jc=1:Nj
    for ec=1:Ne
        for ac=1:Na
            acc=afunG(jc,ec,ac);
            y=yvec(jc,ec);
            c=y+(1+r)*grida(ac)-grida2(acc);
            inc=y+r*grida(ac);
            
            sfun(jc,ec,ac)=inc-c;
            cfun(jc,ec,ac)=c;
            srat(jc,ec,ac)=(inc-c)/inc;
        end
    end
end


% CONSUMPTION BY AGE
cfunJ=zeros(Nj,1);
for jc=1:Nj
    cfunJ(jc)=sum(sum(mea(jc,:,:).*cfun(jc,:,:)))/sum(sum(mea(jc,:,:)));
end


sfunJE=zeros(Nj,Ne);
sratJE=zeros(Nj,Ne);
for jc=1:Nj    
    for ec=1:Ne
        sfunJE(jc,ec)=sum(mea(jc,ec,:).*sfun(jc,ec,:))/sum(mea(jc,ec,:));
        sratJE(jc,ec)=sum(mea(jc,ec,:).*srat(jc,ec,:))/sum(mea(jc,ec,:));
    end
end


age=1:Nj;
age=age+19;

norm=1/cfunJ(1);

minJ=20;
maxJ=80;
figure('Name','Asset (STOCK) BY AGE AND PROD')
hold on
plot(age,norm*afunJE(:,2),'k-','LineWidth',3)
plot(age,norm*afunJE(:,1),'k-.','LineWidth',3)
hold off
xlim([minJ maxJ])
legend('high','low','Location','NW')
grid on 
box on 
xlabel('年齢')
ylabel('資産')
set(gcf,'color','w')
set(gca,'Fontsize',16);
%set(gca,'Fontsize',16,'FontName','Times New Roman');
saveas(gcf,'fig_olg2_a.eps','epsc2')




% SAVING RATE BY ASSET: CHOOSE A PARTICULAR AGE TO PLOT
jc=21;
sfunAE=zeros(Na,Ne);
sratAE=zeros(Na,Ne);
for ac=1:Na    
    for ec=1:Ne
        sfunAE(ac,ec)=sfun(jc,ec,ac);
        sratAE(ac,ec)=srat(jc,ec,ac);
    end
end


minA=0;
maxA=20;
minY=-0.2;
maxY=0.5;
zerovec=zeros(Na,1);
figure('Name','SAVING RATE BY ASSET AND PROD')
hold on
plot(grida,sratAE(:,2),'- k','LineWidth',3)
plot(grida,sratAE(:,1),'-. k','LineWidth',3)
plot(grida,zerovec,': k','LineWidth',2)
hold off
xlim([minA maxA])
%ylim([minY maxY])
legend('high','low','Location','Best')
xlabel('資産')
ylabel('貯蓄率')
grid on 
box on 
set(gcf,'color','w')
set(gca,'Fontsize',16)
%set(gca,'Fontsize',16,'FontName','Times New Roman');
saveas(gcf,'fig_olg2_s.eps','epsc2')



minA=0;
maxA=20;
minY=-0.2;
maxY=0.5;
zerovec=zeros(Na,1);
figure('Name','SAVING LEVEL BY ASSET AND PROD')
hold on
plot(grida,sfunAE(:,2),'- k','LineWidth',3)
plot(grida,sfunAE(:,1),'-. k','LineWidth',3)
plot(grida,zerovec,': k','LineWidth',2)
hold off
xlim([minA maxA])
%ylim([minY maxY])
legend('high','low','Location','Best')
xlabel('資産')
ylabel('貯蓄（水準）')
grid on 
box on 
set(gcf,'color','w')
set(gca,'Fontsize',16)
%set(gca,'Fontsize',16,'FontName','Times New Roman');
%saveas(gcf,'fig_olg2_s.eps','epsc2')




