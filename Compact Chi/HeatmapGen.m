function[J0MeanWErr,LFMeanWErr,J0MeanKSErr,LFMeanKSErr]=HeatmapGen(d,p,N,Chains)
format compact
format short

%Step sizes
DTs=[1/2,1/4,1/8,1/16];
%Integration times
Ts=[1,2,3,4,5];

%Parameter initialization
J0WErrCh=zeros(Chains,1);
LFWErrCh=zeros(Chains,1);
J0MeanWErr=zeros(length(Ts),length(DTs));
LFMeanWErr=zeros(length(Ts),length(DTs));
J0KSErrCh=zeros(Chains,1);
LFKSErrCh=zeros(Chains,1);
J0MeanKSErr=zeros(length(Ts),length(DTs));
LFMeanKSErr=zeros(length(Ts),length(DTs));

%Exact sampling of distribution
ChiSam=AGSamChi(100000,d,p);

for i=1:5
    for j=1:4
        %Select integration time and step size
        T=Ts(i);
        dt=DTs(j);

        %Extract samples using CHMC and HMC-Leapfrog using these parameters
        [qJ0Ch,qLFCh]=ChiSampler(d,p,N,Chains,T,dt);
        for C=1:Chains
            %Computing CDF from discrete samples for CHMC and HMC-Leapfrog
            [yJ0,xJ0]=ecdf(qJ0Ch(C,:));
            [yLF,xLF]=ecdf(qLFCh(C,:));

            %Computing Wasserstein distance for Leapfrog
            a=d/p+1-1/p;
            b=p^(1/p-1)*gamma(d/p)/gamma(d/p+1-1/p);
            if length(xLF)<1000
                %Using discrete version if number of samples is too small for form
                %cdf
                LFWErrCh(C)=ws_distance(qLFCh(C,:),ChiSam);
            else
                LFWErrCh(C)=sum(abs(yLF(2:length(xLF))-gamcdf(xLF(2:length(xLF)).^p/p,d/p,1)).*(xLF(2:length(xLF))-xLF(1:length(xLF)-1)));
            end
            
            %Computing Wasserstein distance for Leapfrog
            if length(xJ0)<1000
                %Using discrete version if number of samples is too small for form
                %cdf
                J0WErrCh(C)=ws_distance(qJ0Ch(C,:),ChiSam);
            else
                J0WErrCh(C)=sum(abs(yJ0(2:length(xJ0))-gamcdf(xJ0(2:length(xJ0)).^p/p,d/p,1)).*(xJ0(2:length(xJ0))-xJ0(1:length(xJ0)-1)));
            end

            %Computing Kolmogorov-Smirnov distances
            J0KSErrCh(C)=max(abs(yJ0-gamcdf(xJ0.^p/p,d/p,1)));
            LFKSErrCh(C)=max(abs(yLF-gamcdf(xLF.^p/p,d/p,1)));

        end

        %Means of KS and W errors across chains
        J0MeanWErr(i,j)=mean(J0WErrCh);
        LFMeanWErr(i,j)=mean(LFWErrCh);

        J0MeanKSErr(i,j)=mean(J0KSErrCh);
        LFMeanKSErr(i,j)=mean(LFKSErrCh);

    end
end
end


