function[qout,pout,Reject]=HMCSolver(q,p,dt,T,PGauss)
format compact

%Hamiltonian
H=@(x,y)sum(x.^PGauss)/PGauss+dot(y,y)*0.5;


%Initialize Vectors
Q=q;
P=p;

N = ceil(T/dt);

%LeapFrog
if PGauss==2
    P=P-dt*Q/2;
for i=1:N-1
    Q=Q+dt*P;
    P=P-dt*Q;
end
Q=Q+dt*P;
P=P-dt*Q/2;
end

if PGauss==4
P=P-dt*Q.*Q.*Q/2;
for i=1:N-1
    Q=Q+dt*P;
    P=P-dt*Q.*Q.*Q;
end
Q=Q+dt*P;
P=P-dt*Q.*Q.*Q/2;
end

if PGauss==6
    P=P-dt*Q.*Q.*Q.*Q.*Q/2;
for i=1:N-1
    Q=Q+dt*P;
    P=P-dt*Q.*Q.*Q.*Q.*Q;
end
Q=Q+dt*P;
P=P-dt*Q.*Q.*Q.*Q.*Q/2;
end

p=-p;
%Acceptance/Rejection
r=min(1,exp(H(q,p)-H(Q,P)));
if rand(1)>r
    Reject=1;
    qout=q;
    pout=p;
else
    qout=Q;
    pout=P;
    Reject=0;
end


end



