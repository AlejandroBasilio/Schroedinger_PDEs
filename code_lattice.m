clear all
close all
deltax=0.025; %Space division
deltat=0.0001; %Time division
x=-3.025:deltax:3.025; %Space points vector

%Data rendered by the program above:
E=5.979;
a=1;
bb=0.05;
m=0.1;
nn=1;
V0=-10;
K=2*pi/6*nn;

%Solving the system of equations:
alfa=sqrt(2*m*E);
beta=sqrt(2*m*(V0+E));

e1=exp(1i*(alfa-K)*(a-bb));
e2=exp(-1i*(alfa+K)*(a-bb));
e3=exp(-1i*(beta-K)*bb);
e4=exp(1i*(beta+K)*bb);
AA=[1 1 -1 -1; alfa -alfa -beta beta; e1 e2 -e3 -e4; (alfa-K)*e1 -(alfa+K)*e2 -(beta-K)*e3 (beta+K)*e4];
BB=AA(2:4,2:4);
abs(det(BB))
Br=-[alfa; e1; (alfa-K)*e1];

A=1;
C=linsolve(BB,Br);

u=zeros(243);
for kk=121:160
    if x(kk)>0
        u(kk)=A*exp(1i*(alfa-K)*x(kk))+C(1)*exp(-1i*(alfa+K)*x(kk));
    else
        u(kk)=C(2)*exp(1i*(beta-K)*x(kk))+C(3)*exp(-1i*(beta+K)*x(kk));
    end
end
for kk=2:40
    u(kk)=u(kk+120);
end
for kk=41:80
    u(kk)=u(kk+80);
end
for kk=81:120
    u(kk)=u(kk+40);
end
for kk=161:200
    u(kk)=u(kk-40);
end
for kk=201:242
    u(kk)=u(kk-80);
end


for kk=2:242
    Psi(1,kk)=exp(i*K*x(kk))*u(kk);
end
Psi(1,1)=Psi(1,241);
Psi(1,243)=Psi(1,3);

b=Psi(1,:)*Psi(1,:)'; %Normalize
z=1/sqrt(b*deltax);
Psi(1,:)=z*Psi(1,:);
t=0:deltat:3; %Time points vector

%Potential
for k=1:length(x)
    if x(k)<=-2.975
        V(k)=10;
    elseif x(k)<-2.075
        V(k)=0;
    elseif x(k)<=-1.975
        V(k)=10;
    elseif x(k)<-1.075
        V(k)=0;
    elseif x(k)<=-0.975
        V(k)=10;
    elseif x(k)<-0.075
        V(k)=0;
    elseif x(k)<=0.075
        V(k)=10;
    elseif x(k)<0.975
        V(k)=0;
    elseif x(k)<=1.075
        V(k)=10;
    elseif x(k)<1.975
        V(k)=0;
    elseif x(k)<=2.075
        V(k)=10;
    elseif x(k)<2.975
        V(k)=0;
    else V(k)=10;
    end
end

%The system of equations:
paso=1

r=1/(2*deltax^2);
a1=(2*j/deltat-2*r)*ones(1,length(x)-2)-V(2:end-1);
a2=(2*j/deltat+2*r)*ones(1,length(x)-2)+V(2:end-1);
M=diag(a1)+diag(r.*ones(1,length(x)-3),-1)+diag(r.*ones(1,length(x)-3),1);
N=diag(a2)+diag(-r.*ones(1,length(x)-3),-1)+diag(-r.*ones(1,length(x)-3),1);
n=zeros(1,length(x)-2)';
n(1)=-r*Psi(1,1);
n(end)=-r*Psi(1,end);

for i=2:length(t)
    Psi(i,2:length(x)-1)=M\(N*Psi(i-1,2:length(x)-1).'+n); %+n si es necesario
    Psi(i,1)=Psi(i,length(x)-2);
    Psi(i,end)=Psi(i,3);
    b=Psi(i,:)*Psi(i,:)';  %Normalize
    z=1/sqrt(b*deltax);
    Psi(i,:)=z*Psi(i,:);
    n(1)=-r*Psi(i,1);      %Born-Von Karman boundary conditions
    n(end)=-r*Psi(i,end);
    i
end

paso=2
for i=1:length(t)
    %computation of probability at each space-time point to plot
    if (i-fix(i/100)*100)==0
        i
    end
    for k=1:length(x)
        MPsi2(i,k)=abs(Psi(i,k))^2;
    end
end

N=4000;
paso=2;

figure

yyaxis left

plot(x,MPsi2(1,:))
axis([-3 3 0 1])
hold 
 
%axis([-1 2 0 5])
for i=1+paso:paso:N

yyaxis right

plot(x,V,'k-')

yyaxis left

plot(x,MPsi2(i-paso,:),'w-')
plot(x,MPsi2(i,:),'r-')
c=num2str(i);
title(['Muestra ' c]) %Title
pause(0.01)
end