close all
deltax=0.01; % Space divisions
deltat=0.001; % Time divisions
x=-12:deltax:6; %Space points vector

for k=1:length(x)
    Psi(1,k)=1/sqrt(.1*sqrt(pi))*exp(sqrt(-1)*10*x(k))*exp(-x(k)^2/4/.01); 
end

b=Psi(1,:)*Psi(1,:)'; %Normalization(Psi*Psi_transpuesta)*deltax=1
z=1/sqrt(b*deltax);
Psi(1,:)=z*Psi(1,:);
t=0:deltat:1; %Time points vector

Psi(2:length(t),1)=0; %Boundary conditions
Psi(2:length(t),length(x))=0;

%Potential
for k=1:length(x)
    if x(k)<0.8
        V(k)=0;
    else
        V(k)=150;
    end
end

%Equation
r=1/(2*deltax^2);
a1=(2*j/deltat-2*r)*ones(1,length(x)-2)-V(2:end-1);
a2=(2*j/deltat+2*r)*ones(1,length(x)-2)+V(2:end-1);
M=diag(a1)+diag(r.*ones(1,length(x)-3),-1)+diag(r.*ones(1,length(x)-3),1);
N=diag(a2)+diag(-r.*ones(1,length(x)-3),-1)+diag(-r.*ones(1,length(x)-3),1);
%n=zeros(1,length(x)-2)';
%n(1,1)=-r*Psi(i,1);  % They are equal to zero on the edges
%n(1,end)=-r*Psi(i,end;)
for i=2:length(t)
    Psi(i,2:length(x)-1)=M\(N*Psi(i-1,2:length(x)-1).'); %+n if necessary
    
    b=Psi(i,:)*Psi(i,:)'; %Normalize (Psi*Psi_transpuesta)*deltax=1
    z=1/sqrt(b*deltax);
    Psi(i,:)=z*Psi(i,:);
end

for i=1:length(t) % % Computation of probability at each space-time point to plot
    if (i-fix(i/100)*100)==0
        i
    end
    for k=1:length(x)
        MPsi2(i,k)=abs(Psi(i,k))^2;
    end
end

probl=0;
for k=1:480
    probl=probl+MPsi2(120,k);
end
probr=0;
for k=481:801
    probr=probr+MPsi2(120,k);
end

N=500;  % Until what sample it is shown
paso=2;

for i=1:paso:N
yyaxis right
plot(x,V,'k-')

yyaxis left
plot(x,MPsi2(i,:),'r-')
c=num2str(i);
axis([-1 2 0 5])
title(['Sample ' c]) % Title
pause(0.03)
end