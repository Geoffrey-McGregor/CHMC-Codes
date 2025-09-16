function[qJ0Ch,qLFCh,qFJCh,qJ0ChTime,qLFChTime,qFJChTime]=ChiSampler(d,p,N,Chains,T,dt)
format compact
format short

%Standard Choice of parameters
%T=5;
%dt=0.05;
%Chains=10;
%N=5000;

%Variable Initialization
qJ0Ch=zeros(Chains,N);
qLFCh=zeros(Chains,N);
qFJCh=zeros(Chains,N);

%Computation Time Tracking
qJ0ChTime=zeros(Chains,N-1);
qLFChTime=zeros(Chains,N-1);
qFJChTime=zeros(Chains,N-1);

%Creates functions (the PDF, the Potential Energy U, its derivative DU and the
%Hamiltonian H
true_func = @(x) (p^(1-d/p))/gamma(d/p)*abs(x).^(d-1).*exp(-abs(x).^p/p);
U=@(x)-(d-1)*log(abs(x))+abs(x)^p/p;
DU=@(x)-(d-1)*(1/abs(x))+p*abs(x)^(p-1)/p;
H=@(x,y)U(x)+y^2/2;

DDU=@(Q,q) (U(Q)-U(q))/(Q-q);
dDDU=@(Q,q) (DU(Q)*(Q-q)-U(Q)+U(q))/(Q-q)^2;

TPolyVec = [-1/8,1/7,-1/6,1/5,-1/4,1/3,-1/2,1];
DTPolyVec = [-7/8,6/7,-5/6,4/5,-3/4,2/3,-1/2];

for Ch=1:Chains
    
    %Option to set the maximum number of iterations on CHMC solve
    FPIMAX=2;

    %Variable Initialization
    qJ0=zeros(1,N);
    qLF=zeros(1,N);
    qFJ=zeros(1,N);
    pJ0=zeros(1,N);
    pLF=zeros(1,N);
    pFJ=zeros(1,N);
    qLRej=0;
    pLRej=0;

    %Starting point for both HMC-Leapfrog and CHMC
    peak=(d-1)^(1/p);
    qJ0(1)=normrnd(peak,1/2);
    qLF(1)=qJ0(1);
    qFJ(1)=qJ0(1);

    counter=0;
    counterJ=0;

    for i=2:N
        
        %Momentum Draw
        p0=randn(1);
        pJ=p0;
        pL=p0;
        pF=p0;

        pJ0(i-1)=p0;
        pLF(i-1)=p0;
        pFJ(i-1)=p0;
        qJ=qJ0(i-1);
        qL=qLF(i-1);
        qF=qFJ(i-1);

        PL=pL;
        QL=qL;
        PJ=pJ;
        QJ=qJ;
        PF=pF;
        QF=qF;
            
        %Leapfrog Solve
        tic;
        PL=PL-dt/2*DU(qL);
        for j=1:T/dt-1
            QL=QL+dt*PL;
            PL=PL-dt*DU(QL);
        end
        QL=QL+dt*PL;
        PL=PL-dt/2*DU(qL);

        r=min(1,exp(H(qL,pL)-H(QL,PL)));
        
        %Acceptance/Rejection of sample
        if rand(1)>r
            qLF(i)=qL;
            counter=counter+1;
            pLF(i)=pL;
            qLRej(counter)=QL;
            pLRej(counter)=PL;
            EnergyDif(counter)=abs(H(qL,pL)-H(QL,PL));
        else
            qLF(i)=QL;
            pLF(i)=PL;
        end
        CT=toc;
        qLFChTime(Ch,i)=CT+qLFChTime(Ch,i-1);

        %CHMC Solve
        tic;
        for j=1:T/dt
            %Splitting Method Initial Guess
            QJ=QJ+pJ*(dt);

            %Newton's method
            jj=0;
            qVec = [1,qJ,qJ^2,qJ^3,qJ^4,qJ^5];
            dqVec = [5,4*qJ,3*qJ^2,2*qJ^3,qJ^4];
            
            while jj<=FPIMAX
                jj=jj+1;
                %Newton
                deltaU = U(QJ)-U(qJ);
                deltaQ = QJ-qJ;
                G=deltaQ-dt*pJ+0.5*dt^2*deltaU/deltaQ;
                Gp=1+0.5*dt^2*(DU(QJ)*deltaQ-deltaU)/deltaQ^2;

                QJ=QJ-G/Gp;
            end
            PJ=2/dt*(QJ-qJ)-pJ;
        end

        r=min(1,exp(H(qJ,pJ)-H(QJ,PJ)));

        %Acceptance/Rejection
        if rand(1)>r
            qJ0(i)=qJ;
            pJ0(i)=pJ;
            counterJ=counterJ+1;
        else
            qJ0(i)=QJ;
            pJ0(i)=PJ;
        end
       CT=toc;
       qJ0ChTime(Ch,i)=CT+qJ0ChTime(Ch,i-1);

       %CHMC Full-J Solve
       tic;
        for j=1:T/dt
            %Splitting Method Initial Guess
            QF=QF+pF*(dt);

            %Newton's method
            jj=0;
            % qVec = [1,qJ,qJ^2,qJ^3,qJ^4,qJ^5];
            % dqVec = [5,4*qJ,3*qJ^2,2*qJ^3,qJ^4];

            while jj<=FPIMAX
                jj=jj+1;
                %Newton
                deltaU = U(QF)-U(qF);
                deltaQ = QF-qF;
                G=deltaQ-dt*pF+0.5*dt^2*deltaU/deltaQ;
                Gp=1+0.5*dt^2*(DU(QF)*deltaQ-deltaU)/deltaQ^2;

                QF=QF-G/Gp;
            end
            PF=2/dt*(QF-qF)-pF;
        end

        %Compute Jacobian
        deltaQ=(QF-qF);
        deltaU=U(QF)-U(qF);
        DetJ=(1+dt^2/4*(-DU(qF)*deltaQ+deltaU)/deltaQ^2)/(1+dt^2/4*(DU(QF)*deltaQ-deltaU)/deltaQ^2);

        r=min(1,DetJ*exp(H(qF,pF)-H(QF,PF)));

        %Acceptance/Rejection
        if rand(1)>r
            qFJ(i)=qF;
            pFJ(i)=pF;
            counterJ=counterJ+1;
        else
            qFJ(i)=QF;
            pFJ(i)=PF;
        end
        CT=toc;
        qFJChTime(Ch,i)=CT+qFJChTime(Ch,i-1);

    end

    qJ0Ch(Ch,:)=qJ0;
    qLFCh(Ch,:)=qLF;
    qFJCh(Ch,:)=qFJ;
end
end

