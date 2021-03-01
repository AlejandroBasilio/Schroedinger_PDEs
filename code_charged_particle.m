close all
deltax=0.01;    %space division
deltat=0.001;   %time division

x=-2:deltax:2;  %space vector

%initial conditions (first index is time)
for k=1:length(x)
    Psi(1,k)=1/sqrt(.1*sqrt(pi))*exp(sqrt(-1)*0*x(k))*exp(-x(k)^2/4/.01); 
end
t=0:deltat:10;                  %time vector
 
Psi(2:length(t),1)=0;           %border conditions
Psi(2:length(t),length(x))=0;

for k=1:length(x)
    if x(k)<-0.4-0.1*deltax
        V(k)=1e3/(x(k)+0.4);
    elseif x(k)>-0.4+0.1*deltax
        V(k)=-1e3/(x(k)+0.4);
    else
        V(k)=V(k-1);
    end
end

%system of equations. if the V is depends on t, this part
%goes inside the for
r=1/(2*deltax^2);
a1=(2*j/deltat-2*r)*ones(1,length(x)-2)-V(2:end-1);  % -V(x) if necessary
a2=(2*j/deltat+2*r)*ones(1,length(x)-2)+V(2:end-1);  % +V(x) if necessary
    
M=diag(a1)+diag(r.*ones(1,length(x)-3),-1)+diag(r.*ones(1,length(x)-3),1);
N=diag(a2)+diag(-r.*ones(1,length(x)-3),-1)+diag(-r.*ones(1,length(x)-3),1);
%if the border conditions are not 0, this part must be added
%n=zeros(1,length(x)-2).';
%n(1,1)=-r*Psi(i,1);  
%n(1,end)=-r*Psi(i,end;)
for i=2:length(t)
    Psi(i,2:length(x)-1)=M\(N*Psi(i-1,2:length(x)-1).');    %+n if necessary
    
    b=Psi(i,:)*Psi(i,:)';   %to rescale, sum(Psi*Psi_transpose)*deltax=1
    z=1/sqrt(b*deltax);
    Psi(i,:)=z*Psi(i,:);
end
 
for i=1:length(t) 
    %computation of probability at each space-time point to plot
    for k=1:length(x)
        prob(i,k)=abs(Psi(i,k))^2;
    end
end
%plot
mesh(x,t,prob);

N=4000;  %Until what paso is simulated
paso=1;

figure
for i=1:paso:N
yyaxis right
plot(x,V,'k-')

yyaxis left
plot(x,prob(i,:),'r-')
axis([-2 2 -5 5])
c=num2str(i);
title(['Sample ' c]) %Title
pause(0.1)
end