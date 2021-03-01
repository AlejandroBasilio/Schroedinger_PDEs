% The following program generates initial conditions for the lattice problem
clear all
%Initial data:
E1=0.5;
E2=9.1;
a=1;
b=0.05;
m=0.1;
n=1;
V0=-10;
K=2*pi/6*n;

inc=0.1;
E=E1;
E_min=E1;
d_min=1e20;

% Minimization of the determinant:
while d_min>1e-20 && inc>1e-20
    clear d;
    clear E;
    n=0;
for kk=E1:inc:E2
    n=n+1;
    E(n)=kk;
    alfa=sqrt(2*m*E(n));
    beta=sqrt(2*m*(E(n)+V0));
    e1=exp(1i*(alfa-K)*(a-b));
    e2=exp(-1i*(alfa+K)*(a-b));
    e3=exp(-1i*(beta-K)*b);
    e4=exp(1i*(beta+K)*b);
    A=[1 1 -1 -1; alfa -alfa -beta beta; e1 e2 -e3 -e4; (alfa-K)*e1 -(alfa+K)*e2 -(beta-K)*e3 (beta+K)*e4];
    d(n)=abs(det(A));
end

 d_min=min(d);
    for nn=1:length(d)
        if d(nn)==d_min
            E_min=E(nn);
        end
    end
    E1=E_min-inc;
    E2=E_min+inc;
    inc=inc/10;
end

% Optimal energy and value of the determinant:
E_min
d_min