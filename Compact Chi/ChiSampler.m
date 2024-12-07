function[qJ0Ch,qLFCh]=ChiSampler(d,p,N,Chains,T,dt)
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

for Ch=1:Chains
    %Creates functions (the PDF, the Potential Energy U, its derivative DU and the
    %Hamiltonian H
    true_func = @(x) (p^(1-d/p))/gamma(d/p)*abs(x).^(d-1).*exp(-abs(x).^p/p);
    U=@(x)-(d-1)*log(abs(x))+abs(x)^p/p;
    DU=@(x)-(d-1)*(1/abs(x))+p*abs(x)^(p-1)/p;
    H=@(x,y)U(x)+y^2/2;
    
    %Option to set the maximum number of iterations on CHMC solve
    FPIMAX=3;

    %Variable Initialization
    qJ0=zeros(1,N);
    qLF=zeros(1,N);
    pJ0=zeros(1,N);
    pLF=zeros(1,N);
    qLRej=0;
    pLRej=0;

    %Starting point for both HMC-Leapfrog and CHMC
    peak=(d-1)^(1/p);
    qJ0(1)=normrnd(peak,1/2);
    qLF(1)=qJ0(1);

    counter=0;
    counterJ=0;

    for i=2:N
        
        %Momentum Draw
        p0=randn(1);
        pJ=p0;
        pL=p0;

        pJ0(i-1)=p0;
        pLF(i-1)=p0;
        qJ=qJ0(i-1);
        qL=qLF(i-1);

        PL=pL;
        QL=qL;
        PJ=pJ;
        QJ=qJ;

        %Leapfrog Solve
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


        %CHMC Solve
        for j=1:T/dt
            %Splitting Method Initial Guess
            QJ=QJ+pJ*(dt);

            %Newton's method
            jj=0;
            while jj<=FPIMAX
                jj=jj+1;
                %Newton
                G=QJ-qJ-dt*pJ+dt^2/2*(U(QJ)-U(qJ))/(QJ-qJ);
                Gp=1+dt^2/(2)*(DU(QJ)*(QJ-qJ)-(U(QJ)-U(qJ)))/(QJ-qJ)^2;
                QJ=QJ-G./Gp;
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

    end

    qJ0Ch(Ch,:)=qJ0;
    qLFCh(Ch,:)=qLF;
end
end

