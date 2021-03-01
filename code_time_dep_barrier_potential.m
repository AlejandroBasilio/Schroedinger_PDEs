close all;
deltax=0.01;
deltat=0.001; 
x=-10:deltax:14; %Space points vector

Psi= zeros(1,length(x));
for k=1:length(x)
    Psi(1,k)=1/sqrt(.1*sqrt(pi))*exp(sqrt(-1)*10*x(k))*exp(-x(k)^2/4/.01); 
end

b=Psi(1,:)*Psi(1,:)'; %Normalize
z=1/sqrt(b*deltax);
Psi(1,:)=z*Psi(1,:);
t=0:deltat:3; %Time points vector

Psi(2:length(t),1)=0; %Boundary conditions
Psi(2:length(t),length(x))=0;

%Potential
V = zeros(length(t), length(x));
for i = 1:length(t) 
    for k = 1:length(x)
         if x(k)<1.5
            V(i,k)=0;
        elseif x(k)>2
            V(i,k)=0;
         else
            V(i,k)=150*cos(2*pi*t(i));
         end
    end
end

r=1/(2*deltax^2);
for i=2:length(t)
    a1=(2*1i/deltat-2*r)*ones(1,length(x)-2)-V(i,2:end-1);
    a2=(2*1i/deltat+2*r)*ones(1,length(x)-2)+V(i,2:end-1);
    M=diag(a1)+diag(r.*ones(1,length(x)-3),-1)+diag(r.*ones(1,length(x)-3),1);
    N=diag(a2)+diag(-r.*ones(1,length(x)-3),-1)+diag(-r.*ones(1,length(x)-3),1);
    Psi(i,2:length(x)-1)=M\(N*Psi(i-1,2:length(x)-1).');
    
    b=Psi(i,:)*Psi(i,:)'; %Normalize
    z=1/sqrt(b*deltax);
    Psi(i,:)=z*Psi(i,:);
end

MPsi2 = zeros(length(t), length(x));
for i=1:length(t) %Probability
    if (i-fix(i/100)*100)==0
        i
    end
    for k=1:length(x)
        MPsi2(i,k)=abs(Psi(i,k))^2;
    end
end

N=2000;  
paso=1;

for i=1:paso:N
yyaxis right 
plot(x,V(i,:),'k-')
ylim([-200 200])

yyaxis left
plot(x,MPsi2(i,:),'r-')
c=num2str(i);
axis([-1 4 -4.5 4.5])
title(['V(t) = 150 cos(2\pit), Sample ' c])
xlabel('x')
ylabel('Probability')
pause(0.001)
end