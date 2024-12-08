function[qout,pout,Reject,alpha,deltaH]=HMCSolver(q,p,dt,T,PGauss)
    format compact

    % Hamiltonian
    H = @(x,y)sum(x.^PGauss)/PGauss+dot(y,y)*0.5;
    
    % Initialize Vectors
    Q = q;
    P = p;
    
    N = ceil(T/dt);
    
    % LeapFrog for pGauss = 2
    if PGauss==2
        P = P-dt*Q/2;
    for i=1:N-1
        Q = Q+dt*P;
        P = P-dt*Q;
    end
    Q = Q+dt*P;
    P = P-dt*Q/2;
    end
    
    % LeapFrog for pGauss = 4
    if PGauss==4
    P = P-dt*Q.*Q.*Q/2;
    for i=1:N-1
        Q = Q+dt*P;
        P = P-dt*Q.*Q.*Q;
    end
    Q = Q+dt*P;
    P = P-dt*Q.*Q.*Q/2;
    end
    
    % LeapFrog for pGauss = 6
    if PGauss==6
        P = P-dt*Q.*Q.*Q.*Q.*Q/2;
    for i=1:N-1
        Q = Q+dt*P;
        P = P-dt*Q.*Q.*Q.*Q.*Q;
    end
    Q = Q+dt*P;
    P = P-dt*Q.*Q.*Q.*Q.*Q/2;
    end
    
    p = -p;

    deltaH = H(q,p)-H(Q,P);
    alpha = min(1,exp(deltaH));
    if rand(1) > alpha
        Reject = 1;
        qout = q;
        pout = p;
    else
        qout = Q;
        pout = P;
        Reject = 0;
    end


end



