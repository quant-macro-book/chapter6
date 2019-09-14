function f=policyA(aX)

global alpha del beta Nj L theta tau ss0 meaJ


a(1)=0; % INITIAL ASSET 
a(Nj+1)=0;

a(2:Nj)=aX;

K=0;
for j=1:Nj
   K=K+a(j)*meaJ(j);
end

r=alpha*(K/L)^(alpha-1)-del; % INTEREST RATE
w=(1-alpha)*(K/L)^alpha; % WAGE

ss=ss0*w;

% EULER EQUATION (u(c)=log(c))
for j=1:Nj-1
   ff=beta*(1+r)*((1+r)*a(j)+w*theta(j)*(1-tau)+ss(j)-a(j+1)); % beta*(1+r)*c, where c=(1+r)a+w*theta*(1-tau)-a' 
   f(j)=(1+r)*a(j+1)+w*theta(j+1)*(1-tau)+ss(j+1)-a(j+2)-ff;   % 0=c'-ff
end


