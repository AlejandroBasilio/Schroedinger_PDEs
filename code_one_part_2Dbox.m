close all
deltax=0.05; % Space divisions
deltat=0.05; % Time divisions
x=-1:deltax:1; % Space points vector
y=-1:deltax:1;
L=length(x);

t=0:deltat:10; % Time points vector
T=length(t);

for k=1:length(x)
    Psi1(1,k)=1/sqrt(.1*sqrt(pi))*exp(sqrt(-1)*0*x(k))*exp(-x(k)^2/4/.01); 
end
for k=1:length(x)
    Psi2(1,k)=1/sqrt(.1*sqrt(pi))*exp(sqrt(-1)*0*x(k))*exp(-x(k)^2/4/.01); 
end

for i = 1:length(x)
    for k = 1:length(x)
           Ps(1,i,k) = Psi1(1,k)*Psi2(1,i);
    end
end

Ps(1:T, 1:L, 1)=0; % Boundary conditions ( = 0)
Ps(1:T, 1:L, L)=0;
Ps(1:T, 1, 1:L)=0;
Ps(1:T, L, 1:L)=0;

b=sum(sum(Ps(1,:,:).*conj(Ps(1,:,:)))); %Normalization (Psi*Psi_transpuesta)*deltax=1
z=1/(sqrt(b)*deltax);
Ps(1,:,:)=z*Ps(1,:,:);

Psi=zeros(1,L^2);
for i=2:L-1
    Psi(1,(i-2)*(L-2)+1:(i-1)*(L-2))=Ps(1,i,2:L-1);
end

% Equation
r=1/(2*deltax^2);
a1=(2*1i/deltat)-4*r;
a2=(2*1i/deltat)+4*r;
M=zeros((L-2)^2, (L-2)^2);
N=zeros((L-2)^2, (L-2)^2);
for k=2:L-1
    M=diag(a1*ones(1,(L-2)^2))+diag(r*ones(1,(L-2)^2-1),-1)+diag(r*ones(1,(L-2)^2-1),1)+diag(r*ones(1,(L-2)^2-L+2),-L+2)+diag(r*ones(1,(L-2)^2-L+2),L-2);
    N=diag(a2*ones(1,(L-2)^2))+diag(-r*ones(1,(L-2)^2-1),-1)+diag(-r*ones(1,(L-2)^2-1),1)+diag(-r*ones(1,(L-2)^2-L+2),-L+2)+diag(-r*ones(1,(L-2)^2-L+2),L-2);
    for p=2:L-2
       M((p-1)*(L-2),(p-1)*(L-2)+1)=0;
       M((p-1)*(L-2)+1,(p-1)*(L-2))=0;
       N((p-1)*(L-2),(p-1)*(L-2)+1)=0;
       N((p-1)*(L-2)+1,(p-1)*(L-2))=0;
    end
end

for i=2:T
    Psi(i,1:(L-2)^2)=M\(N*Psi(i-1,1:(L-2)^2).'); %+n si es necesario
    for k=2:L-1
        Ps(i,k,2:L-1)=Psi(i,(k-2)*(L-2)+1:(k-1)*(L-2));
    end
    
    b=sum(sum(Ps(1,:,:).*conj(Ps(1,:,:)))); %Normalization(Psi*Psi_transpuesta)*deltax=1
    z=1/(sqrt(b)*deltax);
    Ps(1,:,:)=z*Ps(1,:,:);
end
 
prob = zeros(T,L,L);
for i=1:T % Computation of probability at each space-time point to plot
    for k=1:L
        for p=1:L
            prob(i,k,p)=abs(Ps(i,k,p))^2;
        end
    end
end

for i = 1:length(t)
    probt=reshape(prob(i,:,:),[L,L]);
    mesh(x,y,probt);
    axis([-1 1 -1 1 0 8]);
    c=num2str(i);
    title (['Sample ' c]);
    pause(0.1)
end