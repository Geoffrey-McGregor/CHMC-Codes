function[qout,pout,Jac,Reject,Iter,alpha,deltaH]=CHMCSolverJ0(q,p,dt,T,energyTol,maxFPI,PGauss)
% CHMC Code
format compact

% Parameters
d = length(q);
PGauss = 4;

% P-Gaussian Distribution
H = @(x,y)sum(x.^PGauss)/PGauss+dot(y,y)*0.5;
detJ = 1;

% Initialize Vectors
qOld = q;
pOld = p;
H0 = H(qOld,pOld);
Q = zeros(d,1);
P = zeros(d,1);
N = ceil(T/dt);
Iter = zeros(N,1);

for i=1:N
    % Begin with an initial guess
    Q = Q+P*(exp(dt)-1);

    % Newton's Method
    j=1;
    while j<=maxFPI %&& abs(H(Q,P)-H0)>energyTol
        j = j+1;
        G = Q-q-dt*p+dt^2/(PGauss*2)*(Q.*Q+q.*q).*(Q+q);
        Gp = 1+dt^2/(PGauss*2)*(3*Q.*Q+2*Q.*q+q.*q);
        Q = Q-G./Gp;
    end
    P = 2/dt*(Q-q)-p;
    Iter(i) = j;
    q = Q;
    p = P;
end

% Compute Acceptance Probability
deltaH = H(qOld,pOld)-H(Q,P);
alpha = min(1,exp(deltaH)*detJ);
if rand(1) > alpha
    Reject = 1;
    qout = qOld;
    pout = pOld;
else
    Reject = 0;
    qout = Q;
    pout = P;
end

Jac = detJ;
Iter = mean(Iter);
end


