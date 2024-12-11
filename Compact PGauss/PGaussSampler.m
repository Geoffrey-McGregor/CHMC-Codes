function[covMaxErrJ0,covMaxErrLF,KSMaxErrJ0,KSMaxErrLF,W1MaxErrJ0,W1MaxErrLF]=PGaussSampler(dimArray,p,N,T,dt,Chains,numPoints)
    format compact
    format long
    
    % Size of dimArray 
    sizeOfdimArray = length(dimArray);
        
    %Error variable initialization
    covMaxErrJ0 = zeros(sizeOfdimArray,numPoints,Chains);
    covMaxErrLF = zeros(sizeOfdimArray,numPoints,Chains);
    KSMaxErrJ0 = zeros(sizeOfdimArray,numPoints,Chains);
    KSMaxErrLF = zeros(sizeOfdimArray,numPoints,Chains);
    W1MaxErrJ0 = zeros(sizeOfdimArray,numPoints,Chains);
    W1MaxErrLF = zeros(sizeOfdimArray,numPoints,Chains);

    % Selecting p parameter in p-Gaussian
    PGauss=p;
    
    for dimIdx=1:sizeOfdimArray
    
        % Dimension
        d=dimArray(dimIdx);
    
        % CHMC options
        energyTol = 1e-14;
        maxFPI = 2;
    
        %Printing options
        updatePercent = 10; % percent to print progress (set to zero for no updates)
   
        qJ0 = zeros(d,N,Chains);
        pJ0 = zeros(d,N,Chains);
        qLF = zeros(d,N,Chains);
        pLF = zeros(d,N,Chains);
    
        timeLF = zeros(Chains,1);
        timeCHMCJ0 = zeros(Chains,1);

        minLF = zeros(N,Chains);
        deltaHLF = zeros(N,Chains);
        minJ0 = zeros(N,Chains);
        deltaHJ0 = zeros(N,Chains);

        RejectJ0=zeros(N,Chains);
        RejectLF=zeros(N,Chains);
        IterJ0=zeros(N,Chains);

        if updatePercent
                fprintf('\n===========================================\n')
                fprintf(strcat(' Sampling PGauss in d =', 32, ' ', sprintf('%5d',dimArray(dimIdx)),' (', num2str(dimIdx), ' out of', 32, num2str(sizeOfdimArray),')'))
                fprintf('\n===========================================\n\n')
        end
    
        for j=1:Chains
            p=zeros(d,1);
            pLF(:,1,j)=p;
            pJ0(:,1,j)=p;
    
            q=normrnd(0,1,[d,1]);
            qLF(:,1,j)=q;
            qJ0(:,1,j)=q;

            if updatePercent
                fprintf(strcat('Chain #', num2str(j), ' progress at:'))
            end
            for i=2:N
                p=randn(d,1);
    
                % HMC-LF
                tic
                [qLF(:,i,j),pLF(:,i,j),RejectLF(i-1,j),alpha,deltaH] = HMCSolver(qLF(:,i-1,j),p,dt,T,PGauss);
                timeLF(j) = timeLF(j)+toc;
                minLF(i,j) = alpha;
                deltaHLF(i,j) = deltaH;
    
                % CHMC-J0
                tic
                [qJ0(:,i,j),pJ0(:,i,j),JacJ0,RejectJ0(i-1,j),IterJ0(i,j),alpha,deltaH] = CHMCSolverJ0(qJ0(:,i-1,j),p,dt,T,energyTol,maxFPI,PGauss);
                timeCHMCJ0(j) = timeCHMCJ0(j)+toc;
                minJ0(i,j) = alpha;
                deltaHJ0(i,j) = deltaH;
    
                if mod(i,floor(N/100*updatePercent))==0 && updatePercent
                    fprintf(' %2.0f%%',i/N*100)
                end
            end
            if updatePercent
                fprintf("\n")
            end
        end
        if updatePercent
            fprintf("\n")
        end
    
        headingStr = {'Leapfrog','CHMC J0'};
        columnStr = {'Sample Accept %', 'Acceptance Prob. %', 'Mean Energy Error', 'Mean # of FPIs', 'Total time (s)'};
    
        filename = strcat('comparison-info',datestr(now,'_dd-mm-yy_HH-MM-SS'),'.txt');
        fid = fopen(filename,'w');
    
        for j=1:Chains
            % Printing Chain number
            fprintf(strcat('Chain #', num2str(j), ':\n'))
            fprintf(fid, strcat('Chain #', num2str(j), ':\n'));
    
            sampleAccept = 100*[ 1-sum(RejectLF(:,j))./length(RejectLF(:,j)),1-sum(RejectJ0(:,j))./length(RejectJ0(:,j))];
            acceptProb = 100*[ mean(minLF(2:N,j)),  mean(minJ0(2:N,j))];
            energyError = [ mean(abs(deltaHLF(2:N,j))), mean(abs(deltaHJ0(2:N,j)))];
            meanFPIs = [ 1,  mean(IterJ0(2:N,j))];
            timeArray = [ mean(timeLF(j)),  mean(timeCHMCJ0(j))];
    
            % Printing Heading
            fprintf(' %20s ', ' ')
            fprintf(fid,' %20s ', ' ');
            fprintf(' ')
            fprintf(fid,' ');
            fprintf(' %20s |', headingStr{1}, headingStr{2})
            fprintf(fid,' %20s |', headingStr{1}, headingStr{2});
            fprintf('\n')
            fprintf(fid,'\n');
    
            % Printing Column Entries
            fprintf('------------------------------------------------------------------------------------------------------------------------------------------\n')
            fprintf(fid,'------------------------------------------------------------------------------------------------------------------------------------------\n');
            fprintf(' %20s ', columnStr{1})
            fprintf(fid,' %20s ', columnStr{1});
            fprintf('|')
            fprintf(fid,'|');
            fprintf(' %20.3f |', sampleAccept(1), sampleAccept(2))
            fprintf(fid,' %20.3f |', sampleAccept(1), sampleAccept(2));
            fprintf('\n')
            fprintf(fid,'\n');
            fprintf(' %20s ', columnStr{2})
            fprintf(fid,' %20s ', columnStr{2});
            fprintf('|')
            fprintf(fid,'|');
            fprintf(' %20.3f |', acceptProb(1), acceptProb(2))
            fprintf(fid,' %20.3f |', acceptProb(1), acceptProb(2));
            fprintf('\n')
            fprintf(fid,'\n');
            fprintf(' %20s ', columnStr{3})
            fprintf(fid,' %20s ', columnStr{3});
            fprintf('|')
            fprintf(fid,'|');
            fprintf(' %20.3d |', energyError(1), energyError(2))
            fprintf(fid,' %20.3d |', energyError(1), energyError(2));
            fprintf('\n')
            fprintf(fid,'\n');
            fprintf(' %20s ', columnStr{4})
            fprintf(fid,' %20s ', columnStr{4});
            fprintf('|')
            fprintf(fid,'|');
            fprintf(' %20.3f |', meanFPIs(1), meanFPIs(2))
            fprintf(fid,' %20.3f |', meanFPIs(1), meanFPIs(2));
            fprintf('\n')
            fprintf(fid,'\n');
            fprintf(' %20s ', columnStr{5})
            fprintf(fid,' %20s ', columnStr{5});
            fprintf('|')
            fprintf(fid,'|');
            fprintf(' %20.3f |', timeArray(1), timeArray(2))
            fprintf(fid,' %20.3f |', timeArray(1), timeArray(2));
            fprintf('\n')
            fprintf(fid,'\n');
            fprintf('------------------------------------------------------------------------------------------------------------------------------------------\n')
            fprintf(fid,'------------------------------------------------------------------------------------------------------------------------------------------\n');
        end
    
        % Printing Chain number
        fprintf(strcat('Mean values over all chains:\n'))
        fprintf(fid, strcat('Mean values over all chains:\n'));
        sampleAccept = 100*[ 1-sum(sum(RejectLF))./length(RejectLF)/Chains,1-sum(sum(RejectJ0))./length(RejectJ0)/Chains];
        acceptProb = 100*[ mean(mean(minLF(2:N,:))),  mean(mean(minJ0(2:N,:)))];
        energyError = [ mean(mean(abs(deltaHLF(2:N,:)))), mean(mean(abs(deltaHJ0(2:N,:))))];
        meanFPIs = [ 1,  mean(mean(IterJ0(2:N,:)))];
        timeArray = [ mean(mean(timeLF)),  mean(mean(timeCHMCJ0))];
    
        % Printing Heading
        fprintf(' %20s ', ' ')
        fprintf(fid,' %20s ', ' ');
        fprintf(' ')
        fprintf(fid,' ');
        fprintf(' %20s |', headingStr{1}, headingStr{2})
        fprintf(fid,' %20s |', headingStr{1}, headingStr{2});
        fprintf('\n')
        fprintf(fid,'\n');
    
        % Printing Column Entries
        fprintf('------------------------------------------------------------------------------------------------------------------------------------------\n')
        fprintf(fid,'------------------------------------------------------------------------------------------------------------------------------------------\n');
        fprintf(' %20s ', columnStr{1})
        fprintf(fid,' %20s ', columnStr{1});
        fprintf('|')
        fprintf(fid,'|');
        fprintf(' %20.3f |', sampleAccept(1), sampleAccept(2))
        fprintf(fid,' %20.3f |', sampleAccept(1), sampleAccept(2));
        fprintf('\n')
        fprintf(fid,'\n');
        fprintf(' %20s ', columnStr{2})
        fprintf(fid,' %20s ', columnStr{2});
        fprintf('|')
        fprintf(fid,'|');
        fprintf(' %20.3f |', acceptProb(1), acceptProb(2))
        fprintf(fid,' %20.3f |', acceptProb(1), acceptProb(2));
        fprintf('\n')
        fprintf(fid,'\n');
        fprintf(' %20s ', columnStr{3})
        fprintf(fid,' %20s ', columnStr{3});
        fprintf('|')
        fprintf(fid,'|');
        fprintf(' %20.3d |', energyError(1), energyError(2))
        fprintf(fid,' %20.3d |', energyError(1), energyError(2));
        fprintf('\n')
        fprintf(fid,'\n');
        fprintf(' %20s ', columnStr{4})
        fprintf(fid,' %20s ', columnStr{4});
        fprintf('|')
        fprintf(fid,'|');
        fprintf(' %20.3f |', meanFPIs(1), meanFPIs(2))
        fprintf(fid,' %20.3f |', meanFPIs(1), meanFPIs(2));
        fprintf('\n')
        fprintf(fid,'\n');
        fprintf(' %20s ', columnStr{5})
        fprintf(fid,' %20s ', columnStr{5});
        fprintf('|')
        fprintf(fid,'|');
        fprintf(' %20.3f |', timeArray(1), timeArray(2))
        fprintf(fid,' %20.3f |', timeArray(1), timeArray(2));
        fprintf('\n')
        fprintf(fid,'\n');
        fprintf('------------------------------------------------------------------------------------------------------------------------------------------\n')
        fprintf(fid,'------------------------------------------------------------------------------------------------------------------------------------------\n');
    
        % Print parameters
        fprintf(' d = %10d | dt = %10.3f | N = %10d | T = %10.3f | eTOL = %.3d | maxFPI = %d | Chains = %d |\n\n', d, dt, N, T, energyTol, maxFPI, Chains)
        fprintf(fid,' d = %10d | dt = %10.3f | N = %10d | T = %10.3f | eTOL = %.3d | maxFPI = %d | Chains = %d |\n\n', d, dt, N, T, energyTol, maxFPI, Chains);
        fclose(fid);

        % Compute Error Metrics

        %Number of points on convergence plot
        plotStepSize= floor(N/numPoints);
        
        %Exact covariance
        covExact = (PGauss^2)^(1/PGauss)*gamma(3/PGauss)/gamma(1/PGauss);
        if updatePercent
            fprintf(strcat('Computing Covariance Error:\n'))
        end
        
        %Computing covariance error for each chain
        tic
        for j=1:Chains
            fprintf(strcat('Chain #', num2str(j), ':'))
            for i=1:numPoints
                covError=abs(VarComp(qLF(:,1:i*floor(N/numPoints),j))-covExact);
                covMaxErrLF(dimIdx,i,j)=max(covError);
                if updatePercent
                    fprintf(' %2.1f%%',(2*i-1)/2/numPoints*100)
                end
                covError=abs(VarComp(qJ0(:,1:i*floor(N/numPoints),j))-covExact);
                covMaxErrJ0(dimIdx,i,j)=max(covError);
                if updatePercent
                    fprintf(' %2.1f%%',i/numPoints*100)
                end
            end
            if(updatePercent)
                fprintf('\n')
            end
        end
        covCalcTime = toc;
        fprintf("Time to compute Covariance Error: %.3f seconds\n\n",covCalcTime)  
        
        %Computing Wasserstein and Kolmogorov-Smirnov Errors for each chain
        KSVecErrJ0 = zeros(d,1);
        KSVecErrLF = zeros(d,1);
        W1VecErrJ0 = zeros(d,1);
        W1VecErrLF = zeros(d,1);

        if updatePercent
            fprintf(strcat('Computing Wasserstein and Kolmogorov-Smirnov Errors:\n'))
        end
        tic
        for C=1:Chains
            ct=0;
            fprintf(strcat('Chain #', num2str(C), ':'))
            for j=plotStepSize:plotStepSize:N
                ct=ct+1;
                for i=1:d
                    [yJ0cdf,xJ0cdf]=ecdf(qJ0(i,1:j,C));
                    lxJ0 = length(xJ0cdf);
                    KSVecErrJ0(i)=max(abs(yJ0cdf-PGenCDF(xJ0cdf,PGauss)));
                    W1VecErrJ0(i)=sum(abs(yJ0cdf(2:lxJ0)-PGenCDF(xJ0cdf(2:lxJ0),PGauss)).*(xJ0cdf(2:lxJ0)-xJ0cdf(1:lxJ0-1)));
                end
                if updatePercent
                    fprintf(' %2.1f%%',(2*ct-1)/2/numPoints*100)
                end
                for i=1:d
                    [yLFcdf,xLFcdf]=ecdf(qLF(i,1:j,C));
                    lxLF = length(xLFcdf);
                    KSVecErrLF(i)=max(abs(yLFcdf-PGenCDF(xLFcdf,PGauss)));
                    W1VecErrLF(i)=sum(abs(yLFcdf(2:lxLF)-PGenCDF(xLFcdf(2:lxLF),PGauss)).*(xLFcdf(2:lxLF)-xLFcdf(1:lxLF-1)));
                end
                if updatePercent
                    fprintf(' %2.1f%%',ct/numPoints*100)
                end
                KSMaxErrJ0(dimIdx,ct,C)=max(KSVecErrJ0);
                KSMaxErrLF(dimIdx,ct,C)=max(KSVecErrLF);
                W1MaxErrJ0(dimIdx,ct,C)=max(W1VecErrJ0);
                W1MaxErrLF(dimIdx,ct,C)=max(W1VecErrLF);
            end
            if(updatePercent)
                fprintf('\n')
            end
        end
        KSW1CalcTime = toc;
        fprintf("Time to compute KS and W1 Errors: %.3f seconds\n\n",KSW1CalcTime)
    end
end