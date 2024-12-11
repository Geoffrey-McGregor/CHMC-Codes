function[xiJ0d,xiFullJd,xiLFd,JacJ0d,JacFullJd,minJ0d,minFullJd,minLFd,deltaHJ0d,deltaHFullJd,deltaHLFd,RejectionTracker]=PGaussSampler(dimArray,p,N,T,dt)
    format compact
    format long
    
    % Storage of vectors across multiple dimensions
    sizeOfdimArray = length(dimArray);
    JacJ0d=zeros(sizeOfdimArray,N);
    JacFullJd=zeros(sizeOfdimArray,N);
    MaxD=dimArray(sizeOfdimArray);
    Chunk=MaxD/5;
    RejectionTracker=zeros(MaxD/Chunk,3);
    
    xiJ0d = zeros(sizeOfdimArray,N);
    xiFullJd = zeros(sizeOfdimArray,N);
    xiLFd = zeros(sizeOfdimArray,N);
    
    minJ0d = zeros(sizeOfdimArray,N);
    minFullJd = zeros(sizeOfdimArray,N);
    minLFd = zeros(sizeOfdimArray,N);

    deltaHJ0d = zeros(sizeOfdimArray,N);
    deltaHFullJd = zeros(sizeOfdimArray,N);
    deltaHLFd = zeros(sizeOfdimArray,N);
    
    acceptJ0 = cell(5,1);
    acceptFullJ = cell(5,1);
    acceptLF = cell(5,1);
    
    % Selecting p parameter in p-Gaussian
    PGauss=p;
    
    for dimIdx=1:sizeOfdimArray
    
        % Dimension
        d=dimArray(dimIdx);
    
        % Number of chains
        Chains=1;
    
        % CHMC options
        energyTol = 1e-14;
        maxFPI = 2;
    
        %Printing options
        updatePercent = 10; % percent to print progress (set to zero for no updates)
    
        numSteps = ceil(T/dt);
    
        qJ0 = zeros(d,N,Chains);
        pJ0 = zeros(d,N,Chains);
        minJ0 = zeros(N,Chains);
        JacJ0 = zeros(N,Chains);
        deltaHJ0 = zeros(N,Chains);
        IterJ0 = zeros(numSteps,Chains);
        RejectJ0 = zeros(N,Chains);
    
        qFullJ = zeros(d,N,Chains);
        pFullJ = zeros(d,N,Chains);
        minFullJ = zeros(N,Chains);
        JacFullJ = zeros(N,Chains);
        deltaHFullJ = zeros(N,Chains);
        IterFullJ = zeros(numSteps,Chains);
        RejectFullJ = zeros(N,Chains);
    
        qLF = zeros(d,N,Chains);
        pLF = zeros(d,N,Chains);
        minLF = zeros(N,Chains);
        deltaHLF = zeros(N,Chains);
        RejectLF = zeros(N,Chains);
    
        timeLF = zeros(Chains,1);
        timeCHMCJ0 = zeros(Chains,1);
        timeCHMCFullJ = zeros(Chains,1);
    
    
        for j=1:Chains
            p=zeros(d,1);
            pLF(:,1,j)=p;
            pJ0(:,1,j)=p;
            pFullJ(:,1,j)=p;
    
            q=normrnd(0,1,[d,1]);
            qLF(:,1,j)=q;
            qJ0(:,1,j)=q;
            qFullJ(:,1,j)=q;
    
            RejectLF(1,j)=0;
            RejectJ0(1,j)=0;
            RejectFullJ(1,j)=0;
    
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
                [qJ0(:,i,j),pJ0(:,i,j),JacJ0(i-1,j),RejectJ0(i-1,j),IterJ0(i,j),alpha,deltaH] = CHMCSolverJ0(qJ0(:,i-1,j),p,dt,T,energyTol,maxFPI);
                timeCHMCJ0(j) = timeCHMCJ0(j)+toc;
                minJ0(i,j) = alpha;
                deltaHJ0(i,j) = deltaH;
    
                % CHMC-FullJ
                tic
                [qFullJ(:,i,j),pFullJ(:,i,j),JacFullJ(i-1,j),RejectFullJ(i-1,j),IterFullJ(i,j),alpha,deltaH] = CHMCSolverFullJ(qFullJ(:,i-1,j),p,dt,T,energyTol,maxFPI);
                timeCHMCFullJ(j) = timeCHMCFullJ(j)+toc;
                minFullJ(i,j) = alpha;
                deltaHFullJ(i,j) = deltaH;
    
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
    
        headingStr = {'Leapfrog','CHMC J0','CHMC FullJ'};
        columnStr = {'Sample Accept %', 'Acceptance Prob. %', 'Mean Energy Error', 'Mean # of FPIs', 'Total time (s)'};
    
        filename = strcat('comparison-info',datestr(now,'_dd-mm-yy_HH-MM-SS'),'.txt');
        fid = fopen(filename,'w');
    
        for j=1:Chains
            % Printing Chain number
            fprintf(strcat('Chain #', num2str(j), ':\n'))
            fprintf(fid, strcat('Chain #', num2str(j), ':\n'));
    
            sampleAccept = 100*[ 1-sum(RejectLF(:,j))./length(RejectLF(:,j)),1-sum(RejectJ0(:,j))./length(RejectJ0(:,j)),1-sum(RejectFullJ(:,j))./length(RejectFullJ(:,j))];
            acceptProb = 100*[ mean(minLF(2:N,j)),  mean(minJ0(2:N,j)), mean(minFullJ(2:N,j))];
            energyError = [ mean(abs(deltaHLF(2:N,j))), mean(abs(deltaHJ0(2:N,j))), mean(abs(deltaHFullJ(2:N,j)))];
            meanFPIs = [ 1,  mean(IterJ0(2:N,j)), mean(IterFullJ(2:N,j))];
            timeArray = [ mean(timeLF(j)),  mean(timeCHMCJ0(j)), mean(timeCHMCFullJ(j))];
    
            % Printing Heading
            fprintf(' %20s ', ' ')
            fprintf(fid,' %20s ', ' ');
            fprintf(' ')
            fprintf(fid,' ');
            fprintf(' %20s |', headingStr{1}, headingStr{2}, headingStr{3})
            fprintf(fid,' %20s |', headingStr{1}, headingStr{2}, headingStr{3});
            fprintf('\n')
            fprintf(fid,'\n');
    
            % Printing Column Entries
            fprintf('------------------------------------------------------------------------------------------------------------------------------------------\n')
            fprintf(fid,'------------------------------------------------------------------------------------------------------------------------------------------\n');
            fprintf(' %20s ', columnStr{1})
            fprintf(fid,' %20s ', columnStr{1});
            fprintf('|')
            fprintf(fid,'|');
            fprintf(' %20.3f |', sampleAccept(1), sampleAccept(2), sampleAccept(3))
            fprintf(fid,' %20.3f |', sampleAccept(1), sampleAccept(2), sampleAccept(3));
            fprintf('\n')
            fprintf(fid,'\n');
            fprintf(' %20s ', columnStr{2})
            fprintf(fid,' %20s ', columnStr{2});
            fprintf('|')
            fprintf(fid,'|');
            fprintf(' %20.3f |', acceptProb(1), acceptProb(2), acceptProb(3))
            fprintf(fid,' %20.3f |', acceptProb(1), acceptProb(2), acceptProb(3));
            fprintf('\n')
            fprintf(fid,'\n');
            fprintf(' %20s ', columnStr{3})
            fprintf(fid,' %20s ', columnStr{3});
            fprintf('|')
            fprintf(fid,'|');
            fprintf(' %20.3d |', energyError(1), energyError(2), energyError(3))
            fprintf(fid,' %20.3d |', energyError(1), energyError(2), energyError(3));
            fprintf('\n')
            fprintf(fid,'\n');
            fprintf(' %20s ', columnStr{4})
            fprintf(fid,' %20s ', columnStr{4});
            fprintf('|')
            fprintf(fid,'|');
            fprintf(' %20.3f |', meanFPIs(1), meanFPIs(2), meanFPIs(3))
            fprintf(fid,' %20.3f |', meanFPIs(1), meanFPIs(2), meanFPIs(3));
            fprintf('\n')
            fprintf(fid,'\n');
            fprintf(' %20s ', columnStr{5})
            fprintf(fid,' %20s ', columnStr{5});
            fprintf('|')
            fprintf(fid,'|');
            fprintf(' %20.3f |', timeArray(1), timeArray(2), timeArray(3))
            fprintf(fid,' %20.3f |', timeArray(1), timeArray(2), timeArray(3));
            fprintf('\n')
            fprintf(fid,'\n');
            fprintf('------------------------------------------------------------------------------------------------------------------------------------------\n')
            fprintf(fid,'------------------------------------------------------------------------------------------------------------------------------------------\n');
        end
    
        % Printing Chain number
        fprintf(strcat('Mean values over all chains:\n'))
        fprintf(fid, strcat('Mean values over all chains:\n'));
        sampleAccept = 100*[ 1-sum(sum(RejectLF))./length(RejectLF)/Chains,1-sum(sum(RejectJ0))./length(RejectJ0)/Chains,1-sum(sum(RejectFullJ))./length(RejectFullJ)/Chains];
        acceptProb = 100*[ mean(mean(minLF(2:N,:))),  mean(mean(minJ0(2:N,:))),  mean(mean(minFullJ(2:N,:)))];
        energyError = [ mean(mean(abs(deltaHLF(2:N,:)))), mean(mean(abs(deltaHJ0(2:N,:)))), mean(mean(abs(deltaHFullJ(2:N,:))))];
        meanFPIs = [ 1,  mean(mean(IterJ0(2:N,:))), mean(mean(IterFullJ(2:N,:)))];
        timeArray = [ mean(mean(timeLF)),  mean(mean(timeCHMCJ0)), mean(mean(timeCHMCFullJ))];
    
        % Printing Heading
        fprintf(' %20s ', ' ')
        fprintf(fid,' %20s ', ' ');
        fprintf(' ')
        fprintf(fid,' ');
        fprintf(' %20s |', headingStr{1}, headingStr{2}, headingStr{3})
        fprintf(fid,' %20s |', headingStr{1}, headingStr{2}, headingStr{3});
        fprintf('\n')
        fprintf(fid,'\n');
    
        % Printing Column Entries
        fprintf('------------------------------------------------------------------------------------------------------------------------------------------\n')
        fprintf(fid,'------------------------------------------------------------------------------------------------------------------------------------------\n');
        fprintf(' %20s ', columnStr{1})
        fprintf(fid,' %20s ', columnStr{1});
        fprintf('|')
        fprintf(fid,'|');
        fprintf(' %20.3f |', sampleAccept(1), sampleAccept(2), sampleAccept(3))
        fprintf(fid,' %20.3f |', sampleAccept(1), sampleAccept(2), sampleAccept(3));
        fprintf('\n')
        fprintf(fid,'\n');
        fprintf(' %20s ', columnStr{2})
        fprintf(fid,' %20s ', columnStr{2});
        fprintf('|')
        fprintf(fid,'|');
        fprintf(' %20.3f |', acceptProb(1), acceptProb(2), acceptProb(3))
        fprintf(fid,' %20.3f |', acceptProb(1), acceptProb(2), acceptProb(3));
        fprintf('\n')
        fprintf(fid,'\n');
        fprintf(' %20s ', columnStr{3})
        fprintf(fid,' %20s ', columnStr{3});
        fprintf('|')
        fprintf(fid,'|');
        fprintf(' %20.3d |', energyError(1), energyError(2), energyError(3))
        fprintf(fid,' %20.3d |', energyError(1), energyError(2), energyError(3));
        fprintf('\n')
        fprintf(fid,'\n');
        fprintf(' %20s ', columnStr{4})
        fprintf(fid,' %20s ', columnStr{4});
        fprintf('|')
        fprintf(fid,'|');
        fprintf(' %20.3f |', meanFPIs(1), meanFPIs(2), meanFPIs(3))
        fprintf(fid,' %20.3f |', meanFPIs(1), meanFPIs(2), meanFPIs(3));
        fprintf('\n')
        fprintf(fid,'\n');
        fprintf(' %20s ', columnStr{5})
        fprintf(fid,' %20s ', columnStr{5});
        fprintf('|')
        fprintf(fid,'|');
        fprintf(' %20.3f |', timeArray(1), timeArray(2), timeArray(3))
        fprintf(fid,' %20.3f |', timeArray(1), timeArray(2), timeArray(3));
        fprintf('\n')
        fprintf(fid,'\n');
        fprintf('------------------------------------------------------------------------------------------------------------------------------------------\n')
        fprintf(fid,'------------------------------------------------------------------------------------------------------------------------------------------\n');
    
    
        % Print parameters
        fprintf(' d = %10d | dt = %10.3f | N = %10d | T = %10.3f | eTOL = %.3d | maxFPI = %d | Chains = %d |\n\n', d, dt, N, T, energyTol, maxFPI, Chains)
        fprintf(fid,' d = %10d | dt = %10.3f | N = %10d | T = %10.3f | eTOL = %.3d | maxFPI = %d | Chains = %d |\n\n', d, dt, N, T, energyTol, maxFPI, Chains);
        fclose(fid);
    
        acceptJ0{dimIdx} = minJ0;
        acceptFullJ{dimIdx} = minFullJ;
        acceptLF{dimIdx} = minLF;
    
        RejectionTracker(dimIdx,:) = sampleAccept;
        JacJ0d(dimIdx,:) = JacJ0(:,j);
        JacFullJd(dimIdx,:) = JacFullJ(:,j);
        xiJ0d(dimIdx,:) = vecnorm(qJ0,4);
        xiFullJd(dimIdx,:) = vecnorm(qFullJ,4);
        xiLFd(dimIdx,:) = vecnorm(qLF,4);
        minJ0d(dimIdx,:) = minJ0;
        minFullJd(dimIdx,:) = minFullJ;
        minLFd(dimIdx,:) = minLF;
        deltaHLFd(dimIdx,:) = deltaHLF;
        deltaHJ0d(dimIdx,:) = deltaHJ0;
        deltaHFullJd(dimIdx,:) = deltaHFullJ;
    end

end