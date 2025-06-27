%%%%%% BER performance for CP-OFDM, CP-OTFS, and CP-ODSS using a rectangular pulse
%%%%%% shaping waveform

clear
clc
close all
poolsize = 4;

%% ODSS modulation parameters
BW_sys=10e3;
q_geo=1.001^2;%q>=alphaMax^2
M_subcarr=40;% number of subcarriers (symbols in freq. domain)
fL=10e3;fH=fL+BW_sys;

%% OFDM modulation parameters
delta_f_ofdm=BW_sys/M_subcarr;
f_edge_ofdm=fL+(0:M_subcarr)*delta_f_ofdm;
f_center_ofdm=(f_edge_ofdm(1:end-1)+f_edge_ofdm(2:end))/2;
%% Channel parameters
Np = 5;
tauMax=25e-3;% maximum delay spread of the channel
alphaMax=1.001;
m_outList=(0:M_subcarr-1)';m_inList=(0:M_subcarr-1);
SNRlist =0:3:21; %in dB
SNR = 10.^(SNRlist/10);

%% OFDM modulation matrix
Tcp_ofdm=0;%alphaMax*tauMax; 
Fs_new=BW_sys;
time_ofdm=-Tcp_ofdm:1/Fs_new:(1/delta_f_ofdm-1/Fs_new);
freq_ofdm=fL+(0:1:length(time_ofdm)-1)*Fs_new/length(time_ofdm);
temp_rx=time_ofdm.*ones(length(f_center_ofdm),length(time_ofdm));
Q_rx_ofdm=sqrt(delta_f_ofdm/Fs_new)*exp(-1i*2*pi*(f_center_ofdm.').*time_ofdm);
Q_rx_ofdm(temp_rx<0)=0;Q_rx_ofdm(temp_rx>=1/delta_f_ofdm)=0;
temp_tx=(time_ofdm.').*ones(length(time_ofdm),length(f_center_ofdm));
Q_tx_ofdm=sqrt(delta_f_ofdm/Fs_new)*exp(1i*2*pi*(f_center_ofdm).*time_ofdm.');
Q_tx_ofdm(temp_tx<-Tcp_ofdm)=0;Q_tx_ofdm(temp_tx>=1/delta_f_ofdm)=0;
Q_tx_eff_ofdm=Q_tx_ofdm;
Q2_ofdm=(freq_ofdm(2)-freq_ofdm(1))*exp(1i*2*pi*freq_ofdm.*(time_ofdm.'));
Q_rx_eff_ofdm=Q_rx_ofdm*Q2_ofdm;

displayRate = 30/60; %displays/plots are updated @ displayRate (minutes)
tDispl  = tic; % start the timer for displaying results
saveRate = 10; %results are saved @ saveRate (minutes)
tSave = tic; % start the timer for saving results
tExec = tic; % start timer for timing the loop execution

%% chEst parameters
Ntau=50;Nalpha=5;
tauGrid = linspace(0,tauMax,Ntau);%gridding along tau
res_tau=tauGrid(2)-tauGrid(1);
tauGrid=reshape(repmat(tauGrid,Nalpha,1),Nalpha*Ntau,1)';
betaGrid=-(Nalpha-1)/2:(Nalpha-1)/2;
res_beta=betaGrid(2)-betaGrid(1);
betaGrid=repmat(betaGrid,1,Ntau);
q_alpha=alphaMax^(1/abs(betaGrid(1)));
alphaGrid =q_alpha.^betaGrid ;%gridding along alpha
nmseHFig = 28; figure(nmseHFig), pause(0.01);

NMSE_mmse_OMP=zeros(1,length(SNRlist));
NMSE_mmse_NOMP=zeros(1,length(SNRlist));
NMSE_mmse_VBI=zeros(1,length(SNRlist));
NMSE_mmse_FVB=zeros(1,length(SNRlist));
NMSE_mmse_SVB=zeros(1,length(SNRlist));

for iSNR = 1:length(SNRlist) %loop across SNR

    %% P.3) Monte-Carlo Run Parameters
    if SNRlist(iSNR) <= -3
        Ntrials = 1e1;%2*5000;%10000;%20000;%7000;%25000;% # of Monte Carlo trials to estimate BER
        %--------------------------------------------------------------------------
    elseif SNRlist(iSNR) <=0
        Ntrials=1e1;
    elseif SNRlist(iSNR) <=15
        Ntrials=1e2;
    else
        Ntrials=1e2;
    end
    %% Monte Carlo Trials
    nTrialsPerBlock = poolsize;%~~~~~~~~~~~~~~~~~~~~2*poolsize;
    disp(['Number of trials per trial-block : ' num2str(nTrialsPerBlock) ' ...'])

    Ntrials = ceil(Ntrials/nTrialsPerBlock)*nTrialsPerBlock; %make Ntrials the least multiple of "nTrialsPerBlock" greater or equal to "Ntrials"
    disp(...
        ['Number of trials (adjusted to the least multiple of nTrialsPerBlock >= original Ntrials) : '...
        num2str(Ntrials) ' ...']...
        );

    nTrialBlocks = Ntrials/nTrialsPerBlock; %number of trial blocks to run
    disp(['Number of trial-blocks : ' num2str(nTrialBlocks)])

    tExecSnr = tic; %start timer for timing the Monte Carlo run execution for current SNR value


    err1_nmseH_ofdm_OMP=zeros(1,nTrialBlocks);
    err1_nmseH_ofdm_FVB=zeros(1,nTrialBlocks);
    err1_nmseH_ofdm_SVB=zeros(1,nTrialBlocks);

    for iTrialBlock = 1 : nTrialBlocks %run each trial-block

        tExecTrialBlock = tic; %start timer for timing the Monte Carlo run execution for current SNR value

        disp(['SNR = ' num2str(SNRlist(iSNR)) ' dB,'...
            ' trial-block = ' num2str(iTrialBlock)...
            ' of ' num2str(nTrialBlocks) ' trial-blocks.']);

        err2_nmseH_ofdm_OMP=zeros(1,nTrialsPerBlock);
        err2_nmseH_ofdm_FVB=zeros(1,nTrialsPerBlock);
        err2_nmseH_ofdm_SVB=zeros(1,nTrialsPerBlock);

        parfor trial = 1:nTrialsPerBlock

            %                         x=2*(randi(2,tot_sym_no,1)-1)-1;%BPSK symbols generation
            x=2*(randi(2,M_subcarr,1)-1)-1;%BPSK symbols generation

            %% channel
            hp = (randn(Np,1) + 1i*randn(Np,1))/sqrt(2);
            taup=tauMax*rand(Np,1);%tauGrid(ind);%
            betap=sort(-(Nalpha-1)/2+(Nalpha-1)*rand(Np,1));
            alphap=q_alpha.^betap;
            H_act_channel_ofdm=gen_channel_mat(hp,alphap,taup,freq_ofdm,time_ofdm,Q_tx_eff_ofdm,Q_rx_eff_ofdm);
            [A_act,~,~,~,~] = gen_dictionary_mat_data(Q_tx_eff_ofdm*x,freq_ofdm,time_ofdm,Q_rx_eff_ofdm,taup,betap,q_alpha);
            rx_ofdm=A_act*hp;
            P_rx_ofdm=rx_ofdm'*rx_ofdm;
            noiseVar_ofdm=P_rx_ofdm/(SNR(iSNR)*length(rx_ofdm));
            noise_ofdm=sqrt(noiseVar_ofdm/2)*(randn(length(rx_ofdm),1)+1i*randn(length(rx_ofdm),1));
            y_ofdm=rx_ofdm+noise_ofdm;
     
            %% channel estimation
            [A_data_dic_ofdm,~,~,~,~] = gen_dictionary_mat_data(Q_tx_eff_ofdm*x,freq_ofdm,time_ofdm,Q_rx_eff_ofdm,tauGrid,betaGrid,q_alpha);
            % using OMP
            [hsparse_ofdm_OMP,~,support_est_ofdm_OMP]=omp(A_data_dic_ofdm,y_ofdm,Np);
            support_est_ofdm_OMP=sort(support_est_ofdm_OMP);
            alpha_est_ofdm_OMP=alphaGrid(support_est_ofdm_OMP);
            tau_est_ofdm_OMP=tauGrid(support_est_ofdm_OMP);
            h_est_ofdm_OMP=hsparse_ofdm_OMP(support_est_ofdm_OMP);
            H_data_est_ofdm_OMP=gen_channel_mat(h_est_ofdm_OMP,alpha_est_ofdm_OMP,tau_est_ofdm_OMP,freq_ofdm,time_ofdm,Q_tx_eff_ofdm,Q_rx_eff_ofdm);
            err2_nmseH_ofdm_OMP(trial)=norm(H_act_channel_ofdm-H_data_est_ofdm_OMP,'fro')^2/norm(H_act_channel_ofdm,'fro')^2;
            
            %                     FVB
            [hsparse_ofdm_FVB,tauGrid_ofdm_FVB,betaGrid_ofdm_FVB,~]=FVB(y_ofdm,1e-3,1e2,noiseVar_ofdm,1e-6,1e-6,hsparse_ofdm_OMP,tauGrid,betaGrid,Np,res_tau,Q_tx_eff_ofdm*x,freq_ofdm,time_ofdm,Q_rx_eff_ofdm,q_alpha);
            [~,pos_ofdm_FVB]=maxk(hsparse_ofdm_FVB,Ntau*Nalpha,'ComparisonMethod','abs');
            support_est_ofdm_FVB=sort(pos_ofdm_FVB);
            h_est_ofdm_FVB=hsparse_ofdm_FVB(support_est_ofdm_FVB);
            beta_est_ofdm_FVB=betaGrid_ofdm_FVB(support_est_ofdm_FVB);
            alpha_est_ofdm_FVB=q_alpha.^beta_est_ofdm_FVB;
            tau_est_ofdm_FVB=tauGrid_ofdm_FVB(support_est_ofdm_FVB);
            H_est_ofdm_FVB=gen_channel_mat(h_est_ofdm_FVB,alpha_est_ofdm_FVB,tau_est_ofdm_FVB,freq_ofdm,time_ofdm,Q_tx_eff_ofdm,Q_rx_eff_ofdm);
            err2_nmseH_ofdm_FVB(trial)=norm(H_act_channel_ofdm-H_est_ofdm_FVB,'fro')^2/norm(H_act_channel_ofdm,'fro')^2;
            %NVBI
            [h_est_ofdm_SVB,tau_est_ofdm_SVB,beta_est_ofdm_SVB,~] = SVB(y_ofdm,1e-3,1e2,noiseVar_ofdm,1e-6,1e-6,hsparse_ofdm_OMP,tauGrid,betaGrid,q_alpha,Np,res_tau,res_beta,Q_tx_eff_ofdm*x,freq_ofdm,time_ofdm,Q_rx_eff_ofdm);
            alpha_est_ofdm_SVB=q_alpha.^beta_est_ofdm_SVB;
            H_est_ofdm_NVBI=gen_channel_mat(h_est_ofdm_SVB,alpha_est_ofdm_SVB,tau_est_ofdm_SVB,freq_ofdm,time_ofdm,Q_tx_eff_ofdm,Q_rx_eff_ofdm);
            err2_nmseH_ofdm_SVB(trial)=norm(H_act_channel_ofdm-H_est_ofdm_NVBI,'fro')^2/norm(H_act_channel_ofdm,'fro')^2;

        end

        err1_nmseH_ofdm_OMP(iTrialBlock)=mean(err2_nmseH_ofdm_OMP);
        
        err1_nmseH_ofdm_FVB(iTrialBlock)=mean(err2_nmseH_ofdm_FVB);
        err1_nmseH_ofdm_SVB(iTrialBlock)=mean(err2_nmseH_ofdm_SVB);
        toc(tExecTrialBlock),
    end

    NMSE_mmse_OMP(iSNR)=mean(err1_nmseH_ofdm_OMP);
    NMSE_mmse_FVB(iSNR)=mean(err1_nmseH_ofdm_FVB);
    NMSE_mmse_SVB(iSNR)=mean(err1_nmseH_ofdm_SVB);

    disp(['%%%%%%%%% Run Time = ' num2str(toc(tExecSnr)) 's:'...
        ' SNR = ' num2str(SNRlist(iSNR)) ' dB,']);
    set(0, 'CurrentFigure', nmseHFig) , cla,
    semilogy(SNRlist(1:iSNR),NMSE_mmse_OMP(1:iSNR),'bo-','LineWidth',2,'MarkerSize',8);hold on; grid on;
    semilogy(SNRlist(1:iSNR),NMSE_mmse_FVB(1:iSNR),'r*-','LineWidth',2,'MarkerSize',8);hold on; grid on;
    semilogy(SNRlist(1:iSNR),NMSE_mmse_SVB(1:iSNR),'b*-','LineWidth',2,'MarkerSize',8);hold on; grid on;
    title('NMSE plot');
    xlabel('SNR(dB)');ylabel('NMSE');
    legend('OMP','FVB','SVB');
    pause(0.01); %wait a bit for all updates to finish...

end %end of SNR loop
%----------------------------------------------------------------------
disp(['********* Total run time: ' num2str(toc(tExec)) ' seconds *********']);
