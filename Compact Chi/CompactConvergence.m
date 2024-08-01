format compact
format short

%Parameter choices
d=400;
p=6;
Chains=10;
N=5000;   %Samples
T=5;      %Integration time
dt=0.1;  %Time step


%Exact sampling for comparisons
[qJ0Ch,qLFCh]=ChiSampler(d,p,N,Chains,T,dt);
ChiSam=AGSamChi(100000,d,p);

%Initialization
DD=250;
Num=length(DD:DD:N);

J0KSErrCh=zeros(Num,Chains);
LFKSErrCh=zeros(Num,Chains);
J0KSMeanErr=zeros(Num,1);
LFKSMeanErr=zeros(Num,1);

J0WErrCh=zeros(Num,Chains);
LFWErrCh=zeros(Num,Chains);
J0WMeanErr=zeros(Num,1);
LFWMeanErr=zeros(Num,1);

counter=0;
for j=DD:DD:N
    counter=counter+1;
    for i=1:Chains
        %Computing ecdf
        [yJ0,xJ0]=ecdf(qJ0Ch(i,1:j));
        [yLF,xLF]=ecdf(qLFCh(i,1:j));
        a=d/p+1-1/p;
        b=p^(1/p-1)*gamma(d/p)/gamma(d/p+1-1/p);

        %Computing Wasserstein
        J0WErrCh(counter,i)=ws_distance(qJ0Ch(i,1:j)',ChiSam,1);
        LFWErrCh(counter,i)=ws_distance(qLFCh(i,1:j)',ChiSam,1);

        %Computing KS Norm
        J0KSErrCh(counter,i)=max(abs(yJ0-gamcdf(xJ0.^p/p,d/p,1)));
        LFKSErrCh(counter,i)=max(abs(yLF-gamcdf(xLF.^p/p,d/p,1)));
    end
    %Computing means accros chains for each norm
    J0KSMeanErr(counter)=mean(J0KSErrCh(counter,:));
    LFKSMeanErr(counter)=mean(LFKSErrCh(counter,:));
    J0WMeanErr(counter)=mean(J0WErrCh(counter,:));
    LFWMeanErr(counter)=mean(LFWErrCh(counter,:));
end

%Plotting
colorLF = [0, 0.4470, 0.7410];
colorJ0 = [0.4660, 0.6740, 0.1880];
alphaLevel = 0.1;

figure(4)
clf
hold on
semilogy(1,LFKSMeanErr(1),'color',colorLF,'linewidth',2);
semilogy(1,LFWMeanErr(1),'color',colorLF,'linewidth',2,'LineStyle','--');
semilogy(1,J0KSMeanErr(1),'color',colorJ0,'linewidth',2);
semilogy(1,J0WMeanErr(1),'color',colorJ0,'linewidth',2,'LineStyle','--');
for j=1:Chains
    ph = semilogy(DD:DD:N,LFKSErrCh(:,j),'color',colorLF,'linewidth',2);
    ph.Color = [colorLF, alphaLevel];
    ph = semilogy(DD:DD:N,J0KSErrCh(:,j),'color',colorJ0,'linewidth',2);
    ph.Color = [colorJ0, alphaLevel];
    ph = semilogy(DD:DD:N,LFWErrCh(:,j),'color',colorLF,'linewidth',2,'LineStyle','--');
    ph.Color = [colorLF, alphaLevel];
    ph = semilogy(DD:DD:N,J0WErrCh(:,j),'color',colorJ0,'linewidth',2,'LineStyle','--');
    ph.Color = [colorJ0, alphaLevel];
end
semilogy(DD:DD:N,LFKSMeanErr,'color',colorLF,'linewidth',3);
semilogy(DD:DD:N,LFWMeanErr,'color',colorLF,'linewidth',3,'LineStyle','--');
semilogy(DD:DD:N,J0KSMeanErr,'color',colorJ0,'linewidth',3);
semilogy(DD:DD:N,J0WMeanErr,'color',colorJ0,'linewidth',3,'LineStyle','--');
set(gca, 'YScale', 'log')
hold off
legend('HMC-LF (KS Distance)', 'HMC-LF ($W_1$ Distance)', 'CHMC (KS Distance)','CHMC ($W_1$ Distance)','interpreter','latex')
convPlotStr = strcat('BothConvPlotChi-d',num2str(d),'p',num2str(p),datestr(now,'_dd-mm-yy_HH-MM-SS'));
xlim([250 5000])
ylim([10^(-3.2) 10^(0.37)])
grid on
print('-dpng','-r400',convPlotStr)