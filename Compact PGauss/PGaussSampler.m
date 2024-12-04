clearvars
format compact
format long

%Selecting p (Choose 2,4 or 6) 
PGauss=4;

%Number of dimensions (Figure 4 includes 1024, 2560, 5120, 10240, 20480 or 40960)
d=128;

%Number of iterations
N=5000;

%Warmup period
Nw=10;

%Number of chains
Chains=10;

%CHMC options
energyTol = 1e-14;
maxFPI = 2;

%Printing options
updatePercent = 10; % percent to print progress (set to zero for no updates)

%Integration Parameters (If PGauss=6, need to reduce dt for stability)
dt=.1;
T=4;
numSteps = ceil(T/dt);

%Initializing variables
qJ0 = zeros(d,N,Chains);
pJ0 = zeros(d,N,Chains);
minJ0=zeros(N,Chains);
JacJ0=zeros(N,Chains);
IterJ0=zeros(numSteps,Chains);
RejectJ0=zeros(N,Chains);
J0KSErr=zeros(d,Chains);
LFKSErr=zeros(d,Chains);
J0WErr=zeros(d,Chains);
LFWErr=zeros(d,Chains);
qLF = zeros(d,N,Chains);
pLF = zeros(d,N,Chains);
minLF=zeros(N,Chains);
RejectLF=zeros(N,Chains);

timeLF = zeros(Chains,1);
timeCHMCJ0 = zeros(Chains,1);

energyErrLF=zeros(N,Chains);
energyErrJ0=zeros(N,Chains);

%Hamiltonian function
H=@(x,y)sum(x.^PGauss)/PGauss+dot(y,y)/2;

for j=1:Chains
    %Initialization of momentum
    p=zeros(d,1);
    pLF(:,1,j)=p;
    pJ0(:,1,j)=p;

    %Initialization of q-variable
    q=normrnd(0,1,[d,1]);
    qLF(:,1,j)=q;
    qJ0(:,1,j)=q;

    RejectLF(1,j)=0;
    RejectJ0(1,j)=0;

    if updatePercent
        fprintf(strcat('Chain #', num2str(j), ' progress at:'))
    end
    for i=2:N
        % Draw momentum
        p=randn(d,1);

        % Obtain sample using HMC-Leapfrog
        tic
        [qLF(:,i,j),pLF(:,i,j),RejectLF(i-1,j)]=HMCSolver(qLF(:,i-1,j),p,dt,T,PGauss);
        timeLF(j) = timeLF(j)+toc;
        energyErrLF(i,j) = H(qLF(:,i-1,j),p)-H(qLF(:,i,j),pLF(:,i,j));
        minLF(i,j)=min(1,exp(energyErrLF(i,j)));


        % Obtain sample using CHMC
        tic
        [qJ0(:,i,j),pJ0(:,i,j),JacJ0(i-1,j),RejectJ0(i-1,j),IterJ0(i,j)]=CHMCVectorSolver(qJ0(:,i-1,j),p,dt,T,energyTol,maxFPI,PGauss);
        timeCHMCJ0(j) = timeCHMCJ0(j)+toc;

        energyErrJ0(i,j) = H(qJ0(:,i-1,j),p)-H(qJ0(:,i,j),pJ0(:,i,j));
        minJ0(i,j)=min(1,exp(energyErrJ0(i,j)));


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
    energyError = [ mean(abs(energyErrLF(2:N,j))), mean(abs(energyErrJ0(2:N,j)))];
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
energyError = [ mean(mean(abs(energyErrLF(2:N,:)))), mean(mean(abs(energyErrJ0(2:N,:))))];
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


% Setting plot colours for Leapfrog and CHMC
colorLF = [0, 0.4470, 0.7410];
colorJ0 = [0.4660, 0.6740, 0.1880];


alphaLevel = 0.1;

clf
%Number of points on convergence plot
numPoints = 20;

%Error variable initialization
covErrLF=zeros(numPoints,Chains);
covErrJ0=zeros(numPoints,Chains);

%Exact covariance
covExact = (PGauss^2)^(1/PGauss)*gamma(3/PGauss)/gamma(1/PGauss);
if updatePercent
    fprintf(strcat('Computing Covariance Error:\n'))
end

%Computing covariance error for each chain
for j=1:Chains
    fprintf(strcat('Chain #', num2str(j), ':'))
    for i=1:numPoints
        covError=abs(VarComp(qLF(:,1:i*floor(N/numPoints),j))-covExact);
        covErrLF(i,j)=max(covError);
        if updatePercent
            fprintf(' %2.1f%%',(2*i-1)/2/numPoints*100)
        end
        covError=abs(VarComp(qJ0(:,1:i*floor(N/numPoints),j))-covExact);
        covErrJ0(i,j)=max(covError);
        if updatePercent
            fprintf(' %2.1f%%',i/numPoints*100)
        end
    end
    if(updatePercent)
        fprintf('\n')
    end
end
covCalcTime = toc;
if updatePercent
    fprintf("\n")
end
fprintf("Time to compute Covariance: %.3f seconds\n\n",covCalcTime)
meanCELF = mean(covErrLF');
meanCEJ0 = mean(covErrJ0');

%%%%%%%%%%%%%%% Accuracy computation %%%%%%%%%%%%%%%%%%%%
Acc=0.15;
LFIndAcc=find(meanCELF-Acc<0,1)*floor(N/numPoints);
J0IndAcc=find(meanCEJ0-Acc<0,1)*floor(N/numPoints);

%Use exact sampling for Wasserstein computation
YRej=AGSam(100000,PGauss);

DD=N/20;
%Variable initialization
J0KSMaxErr=zeros(length(DD:DD:N),Chains);
LFKSMaxErr=zeros(length(DD:DD:N),Chains);
J0WMaxErr=zeros(length(DD:DD:N),Chains);
LFWMaxErr=zeros(length(DD:DD:N),Chains);

%Computing Wasserstein and Kolmogorov-Smirnov Errors accross all chains
for C=1:Chains
    tic
    ct=0;
    X=['Currently on chain ', num2str(C), ' of the Error Computation'];
    disp(X)
    %C
    for j=DD:DD:N
        ct=ct+1;
        for i=1:d
            [yJ0cdf,xJ0cdf]=ecdf(qJ0(i,1:j,C));
            [yLFcdf,xLFcdf]=ecdf(qLF(i,1:j,C));
                
            J0KSErr(i,C)=max(abs(yJ0cdf-PGenCDF(xJ0cdf,PGauss)));
            LFKSErr(i,C)=max(abs(yLFcdf-PGenCDF(xLFcdf,PGauss)));

            yJ0W(i,C)=sum(abs(yJ0cdf(2:length(xJ0cdf))-PGenCDF(xJ0cdf(2:length(xJ0cdf)),PGauss)).*(xJ0cdf(2:length(xJ0cdf))-xJ0cdf(1:length(xJ0cdf)-1)));
            yLFW(i,C)=sum(abs(yLFcdf(2:length(xLFcdf))-PGenCDF(xLFcdf(2:length(xLFcdf)),PGauss)).*(xLFcdf(2:length(xLFcdf))-xLFcdf(1:length(xLFcdf)-1)));
            
        end
        J0KSMaxErr(ct,C)=max(J0KSErr(:,C));
        LFKSMaxErr(ct,C)=max(LFKSErr(:,C));
        J0WMaxErr(ct,C)=max(yJ0W(:,C));
        LFWMaxErr(ct,C)=max(yLFW(:,C));
    end
    toc
end

%Computing error means accross chains for plotting
meanJ0KS=mean(J0KSMaxErr');
meanLFKS=mean(LFKSMaxErr');
meanJ0W=mean(J0WMaxErr');
meanLFW=mean(LFWMaxErr');

%Plotting
figure(1)
hold on
semilogy(1,meanLFKS(1),'color',colorLF,'linewidth',3);
semilogy(1,meanLFW(1),'color',colorLF,'linewidth',3,'LineStyle','--');
semilogy(1,meanCELF(1),'color',colorLF,'linewidth',3,'LineStyle',':');

semilogy(1,meanJ0KS(1),'color',colorJ0,'linewidth',3);
semilogy(1,meanJ0W(1),'color',colorJ0,'linewidth',3,'LineStyle','--');
semilogy(1,meanCEJ0(1),'color',colorJ0,'linewidth',3,'LineStyle',':');

for j=1:Chains
    ph = semilogy(DD:DD:N,LFKSMaxErr(:,j),'color',colorLF,'linewidth',2);
    ph.Color = [colorLF, alphaLevel];
    ph = semilogy(DD:DD:N,J0KSMaxErr(:,j),'color',colorJ0,'linewidth',2);
    ph.Color = [colorJ0, alphaLevel];
    ph = semilogy(DD:DD:N,LFWMaxErr(:,j),'color',colorLF,'linewidth',2,'LineStyle','--');
    ph.Color = [colorLF, alphaLevel];
    ph = semilogy(DD:DD:N,J0WMaxErr(:,j),'color',colorJ0,'linewidth',2,'LineStyle','--');
    ph.Color = [colorJ0, alphaLevel];
    ph = semilogy(floor(N/numPoints):floor(N/numPoints):N,covErrLF(:,j),'color',colorLF,'linewidth',2,'LineStyle',':');
    ph.Color = [colorLF, alphaLevel];
    ph = semilogy(floor(N/numPoints):floor(N/numPoints):N,covErrJ0(:,j),'color',colorJ0,'linewidth',2,'LineStyle',':');
    ph.Color = [colorJ0, alphaLevel];
end
semilogy(DD:DD:N,meanLFKS,'color',colorLF,'linewidth',3);
semilogy(DD:DD:N,meanLFW,'color',colorLF,'linewidth',3,'LineStyle','--');
semilogy(floor(N/numPoints):floor(N/numPoints):N,meanCELF,'color',colorLF,'linewidth',3,'LineStyle',':');

semilogy(DD:DD:N,meanJ0KS,'color',colorJ0,'linewidth',3);
semilogy(DD:DD:N,meanJ0W,'color',colorJ0,'linewidth',3,'LineStyle','--');
semilogy(floor(N/numPoints):floor(N/numPoints):N,meanCEJ0,'color',colorJ0,'linewidth',3,'LineStyle',':');
xlim([250 5000])
xlabel("HMC Iterations", 'interpreter','latex')
set(gca, 'YScale', 'log')
hold off
legend('HMC-LF (KS Distance)', 'HMC-LF ($W_1$ Distance)', 'HMC-LF (Covariance)', 'CHMC (KS Distance)','CHMC ($W_1$ Distance)','CHMC (Covariance)','interpreter','latex')
convPlotStr = strcat('AllConvPlot-d',num2str(d),'-maxFPI',num2str(maxFPI),datestr(now,'_dd-mm-yy_HH-MM-SS'));
grid on
print('-dpng','-r400',convPlotStr)

