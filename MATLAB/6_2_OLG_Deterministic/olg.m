clear all; clc; close all;

global alpha del beta Nj L theta tau ss0 meaJ

%==========================================================================
% MATLAB CODE WRITTEN FOR CH7 OF KEIZAI SEMINAR
% WRITTEN BY SAGIRI KITAO
% 09/14/2019
%==========================================================================


ind_E=1;
    %=1 BASELINE
    %=2 rho 50% down
    %=3 raise NRA TO 70
    %=4 live up to 85

% PARAMETERS
alpha=0.40;
del=0.08;
beta=0.98; 
Nj=61;     % LIVE 60 YEARS MAX (80YRS OLD)
NjW=45;    % WORKING YEARS (RETIRE AT NjW+1) : ENTER AT 21 AND RETIRE AT 65
rho=0.5;   % SS REPLACEMENT RATE (0.0 OR 0.5)

meaJ=(1/Nj)*ones(Nj,1); % AGE DISTRIBUTION (no population growth) 

if ind_E==2
    rho=0.25
end

if ind_E==3
    NjW=50
end

if ind_E==4
   Nj=66
   meaJ=(1/61)*ones(Nj,1); % IN EXPERIMENT TO RAISE LONGEVITY ：人口も増加
end

theta=zeros(Nj,1);
theta(1:NjW)=1.0;

L=sum(meaJ(1:NjW).*theta(1:NjW));

ss0=zeros(Nj,1);
ss0(NjW+1:Nj)=rho; 

tau=rho*(Nj-NjW)/NjW;

% INITIAL GUESS FOR ASSET DECISIONS
a_g=zeros(Nj+1,1);
a_g=zeros(Nj-1,1); % a(2) to a(Nj) (NOTE a(1)=a(Nj+1)=0)
a_g(1)=0.01;
for j=1:Nj-2
    a_g(j+1)=a_g(j)+0.01;
end

% FIND SOLUTIONS : aX(1:Nj-1)
aX=fsolve(@policyA,a_g);

a=zeros(Nj+1,1);
a(2:Nj)=aX;

K=0;
for j=1:Nj
   K=K+a(j)*meaJ(j);
end

r=alpha*(K/L)^(alpha-1)-del; % INTEREST RATE
w=(1-alpha)*(K/L)^alpha; % WAGE

ss=ss0*w;

c=zeros(Nj,1);
for j=1:Nj
    c(j)=w*theta(j)*(1-tau)+(1+r)*a(j)+ss(j)-a(j+1);
end

ageA=1:Nj+1;
ageA=ageA+19;

ageC=1:Nj;
ageC=ageC+19;

minJ=20;
maxJ=19+Nj;

norm=c(1);

save('fig_6_1.mat');

figure('Name','Asset')
plot(ageA,a/norm,'LineWidth',3)
%title('Asset : a_{j}')
xlim([minJ maxJ+1])
xlabel('年齢')
grid on 
box on 
set(gca,'Fontsize',16,'FontName','Times New Roman');
saveas(gcf,'fig_a.eps','epsc2')

figure('Name','Consumption')
plot(ageC,c(1:Nj)/norm,'LineWidth',3)
%title('Consumption : c_{j}')
xlim([minJ maxJ])
xlabel('年齢')
grid on 
box on 
set(gca,'Fontsize',16,'FontName','Times New Roman');
saveas(gcf,'fig_c.eps','epsc2')
