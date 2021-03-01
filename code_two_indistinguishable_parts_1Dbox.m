close all
deltax=0.05;
deltat=0.05; 

% Quantum state of the particles, limits of the box and symmetry of the
% function is fixed. Symmetric (A = 0); Antisymmetric (A = 1).
na = 3;
nb = 2;
a = 2;
A = 1;

x1=-a/2:deltax:a/2; % Space points vector.
Ea = (na^2*pi^2*1/(2*a^2));
ka = sqrt(2*Ea);
x2 = -a/2:deltax:a/2;
Eb = (nb^2*pi^2*1/(2*a^2));
kb = sqrt(2*Eb);
c = sqrt(2/a);  % Normalization constant
C = sqrt(1/2);  % Normalization constant for the whole function (symmetric or antisymmetric)
L=length(x1);

if mod(na,2) == 1
     Psi1a(1,1:length(x1))=c*cos(ka*x1); %Initial conditions (the first index is the time)
     Psi2a(1,1:length(x2))=c*cos(ka*x2);
else
     Psi1a(1,1:length(x1))=c*sin(ka*x1); 
     Psi2a(1,1:length(x2))=c*sin(ka*x2);
end

if mod(nb,2) == 1
     Psi1b(1,1:length(x1))=c*cos(kb*x1); 
     Psi2b(1,1:length(x2))=c*cos(kb*x2); 
else
     Psi1b(1,1:length(x1))=c*sin(kb*x1);
     Psi2b(1,1:length(x2))=c*sin(kb*x2); 
end

for i = 1:length(x1)
    for k = 1:length(x2)
        if A == 0
             Ps(1,i,k) = C*(Psi1a(1,i)*Psi2b(1,k)+Psi2a(1,k)*Psi1b(1,i));
        else
             Ps(1,i,k) = C*(Psi1a(1,i)*Psi2b(1,k)-Psi2a(1,k)*Psi1b(1,i));
        end
    end
end


t=0:deltat:10; % Time points vector
T=length(t);

Ps(1:T, 1:L, 1)=0; % Boundary conditions
Ps(1:T, 1:L, L)=0;
Ps(1:T, 1, 1:L)=0;
Ps(1:T, L, 1:L)=0;

Psi=zeros(1,L^2);
for i=2:L-1
    Psi(1,(i-2)*(L-2)+1:(i-1)*(L-2))=Ps(1,i,2:L-1);
end

% Equation 
r=1/(2*deltax^2);
a1=(2*j/deltat)-4*r;
a2=(2*j/deltat)+4*r;
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
    Psi(i,1:(L-2)^2)=M\(N*Psi(i-1,1:(L-2)^2).'); %+n if necessary
    for k=2:L-1
        Ps(i,k,2:L-1)=Psi(i,(k-2)*(L-2)+1:(k-1)*(L-2));
    end
    
    b=sum(sum(Ps(1,:,:).*conj(Ps(1,:,:)))); %Normalization (Psi*Psi_transpuesta)*deltax=1
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

for i = 1:length(t) % Plot along time (it is stable, so it does not change along time)
    probt=reshape(prob(i,:,:),[L,L]);
    mesh(x1,x2,probt);
    axis([-a/2 a/2 -a/2 a/2 0 2]);
    n1=num2str(na);
    n2=num2str(nb);
    if A == 0
        title (['Bosons n_a = ' n1 ', n_b = ' n2]);
    else
        title (['Fermions n_a = ' n1 ', n_b = ' n2]);
    end
    xlabel('x_1')
    ylabel('x_2')
    zlabel('Probability(x_1,x_2)')
    pause(0.1)
end