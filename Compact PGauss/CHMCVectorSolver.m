function[qout,pout,Jac,Reject,Iter]=CHMCVectorSolver(q,p,dt,T,energyTol,maxFPI,PGauss)
format compact

%In this version of the code we utilize a max FPI, but not an energy tolerance

%Parameters
d=length(q);


%P-Gaussian Distribution
H=@(x,y)sum(x.^PGauss)/PGauss+dot(y,y)*0.5;
detJ=1;

%Initialize Vectors
qOld = q;
pOld = p;
Q=zeros(d,1);
P=zeros(d,1);

 N = ceil(T/dt);
 Iter = zeros(N,1);

for i=1:N
    %Splitting method initial guess
    Q=Q+P*(exp(dt)-1);

    j=1;
    while j<=maxFPI
        j=j+1;
        %Newton
        if PGauss==2
        G=Q-q-dt*p+dt^2/(PGauss*2)*(Q+q);
        Gp=1+dt^2/(PGauss*2)*(1);
        end
        if PGauss==4
        G=Q-q-dt*p+dt^2/(PGauss*2)*(Q.*Q+q.*q).*(Q+q);
        Gp=1+dt^2/(PGauss*2)*(3*Q.*Q+2*Q.*q+q.*q);
        end
        if PGauss==6
        Q2=Q.*Q;
        q2=q.*q;
        QF1=Q2.*Q+q2.*q;
        QF2=Q2+Q.*q+q2;
        G=Q-q-dt*p+dt^2/(PGauss*2)*(QF1).*(QF2);
        Gp=1+dt^2/(PGauss*2)*((3*Q2).*(QF2)+(QF1).*(2*Q+q));
        end
        Q=Q-G./Gp;
    end
        P=2/dt*(Q-q)-p;
    Iter(i)=j;
    
    q=Q;
    p=P;
end

r=min(1,exp(H(qOld,pOld)-H(Q,P))*detJ);
if rand(1)>r
    Reject=1;
    qout=qOld;
    pout=pOld;

else
    qout=Q;
    pout=P;
    Reject=0;
end
Jac=detJ;
Iter=mean(Iter);
end



