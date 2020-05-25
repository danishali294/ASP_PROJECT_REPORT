clear
clc


%------------------------System Parameters---------------------------------
Num_BS_Antennas=  2^7; % BS antennas
BSAntennas_Index=0:1:Num_BS_Antennas-1; % Indices of the BS Antennas

Num_MS_Antennas=  2^5; % MS antennas
MSAntennas_Index=0:1:Num_MS_Antennas-1; % Indices of the MS Antennas


Num_Stream = 4;

DFT_BS = DFT_Codebook(Num_BS_Antennas,1:Num_BS_Antennas);
DFT_MS = DFT_Codebook(Num_MS_Antennas,1:Num_MS_Antennas);

%---------------------- Simulation Parameters-------------------------------
num_trial = 2000; % Number of independent realizations (to be averaged over)


sim_snr      = [-5]; 
sim_measure  = 20:10:160;   %% Number of Training Beams
sim_sector   = [6];     %% {1,360} {2,180} {3,120} {4,90} {6,60} {8,45}
sim_res      = [4];

rng(28);
%----------------------Ranodm Configuation-------------------------------

FULL_RANDOM = MatrixEnsemble(Num_BS_Antennas*Num_MS_Antennas,Num_BS_Antennas*Num_MS_Antennas,'RSE')/sqrt(Num_BS_Antennas*Num_MS_Antennas);

SpreadSEQ_BS2 = FZC(Num_BS_Antennas,43);
SpreadSEQ_MS2 = FZC(Num_MS_Antennas,11);

PN = FZC(Num_BS_Antennas/4,11);
SpreadSEQ_BS3 = zeros(Num_BS_Antennas,1);
for iii = 1:1:length(PN)
  SpreadSEQ_BS3(4*iii-3) = PN(iii);  
  SpreadSEQ_BS3(4*iii-2) = PN(iii); 
  SpreadSEQ_BS3(4*iii-1) = PN(iii);
  SpreadSEQ_BS3(4*iii)   = PN(iii);
end
SpreadSEQ_MS3 = FZC(Num_MS_Antennas,11);


PN = FZC(Num_BS_Antennas/8,11);
SpreadSEQ_BS4 = zeros(Num_BS_Antennas,1);
for iii = 1:1:length(PN)
  SpreadSEQ_BS4(8*iii-7) = PN(iii);  
  SpreadSEQ_BS4(8*iii-6) = PN(iii); 
  SpreadSEQ_BS4(8*iii-5) = PN(iii);
  SpreadSEQ_BS4(8*iii-4) = PN(iii);
  SpreadSEQ_BS4(8*iii-3) = PN(iii);  
  SpreadSEQ_BS4(8*iii-2) = PN(iii); 
  SpreadSEQ_BS4(8*iii-1) = PN(iii);
  SpreadSEQ_BS4(8*iii)   = PN(iii);
end
SpreadSEQ_MS4 = FZC(Num_MS_Antennas,11);



avg_rateP = zeros(size(sim_snr,2),size(sim_measure,2));
avg_rate0 = zeros(size(sim_snr,2),size(sim_measure,2));
avg_rate1 = zeros(size(sim_snr,2),size(sim_measure,2));
avg_rate2 = zeros(size(sim_snr,2),size(sim_measure,2));
avg_rate3 = zeros(size(sim_snr,2),size(sim_measure,2));
avg_rate4 = zeros(size(sim_snr,2),size(sim_measure,2));
avg_rate5 = zeros(size(sim_snr,2),size(sim_measure,2));

for run_snr = 1:1:length(sim_snr)
for run_m   = 1:1:length(sim_measure)
    
run_sector  = 1; 
run_res     = 1;
    
%--------------------Number of Measurments--------------------------------
num_measure = sim_measure(run_m);

%--------------------Noise Configuration----------------------------------
No          = 1/10^(sim_snr(run_snr)/10); %% Noise Power
%--------------------Dictionary Generation----------------------------------
Tx_Resolusion = Num_BS_Antennas*sim_res(run_res);
Rx_Resolusion = Num_MS_Antennas*sim_res(run_res);
AbG = zeros(Num_BS_Antennas,Tx_Resolusion);
AmG = zeros(Num_MS_Antennas,Rx_Resolusion);
for g=1:1:Num_BS_Antennas
    AbG(g,:)=sqrt(1/Num_BS_Antennas)*exp(-1j*(2*pi)*(g-1)*((0:1:(Tx_Resolusion-1))/Tx_Resolusion));
end
for g=1:1:Num_MS_Antennas
    AmG(g,:)=sqrt(1/Num_MS_Antennas)*exp(-1j*(2*pi)*(g-1)*((0:1:(Rx_Resolusion-1))/Rx_Resolusion));
end

%--------------------Covergae Configuration----------------------------------
Sector         =   sim_sector(run_sector); 

%--------------------Codebook Generation----------------------------------
Search_Range     = [1:1:ceil(Tx_Resolusion/Sector),   (Tx_Resolusion-ceil(Tx_Resolusion/Sector)+1):1:Tx_Resolusion      ];
Search_Range_DFT = [1:1:ceil(Num_BS_Antennas/Sector), (Num_BS_Antennas-ceil(Num_BS_Antennas/Sector)+1):1:Num_BS_Antennas];
DFT_BS_INTEND = DFT_BS(:,Search_Range_DFT);

% % Exhaustive Search
P0 = DFT_BS_INTEND;
Q0 = DFT_MS;
MEA0   = kron(transpose(P0),Q0');
MEA0O  = MEA0*kron(conj(AbG(:,Search_Range)),AmG);

% Full Random CS Codebook Design
MEA1   = FULL_RANDOM(1:num_measure,:);
MEA1O  = MEA1*kron(conj(AbG(:,Search_Range)),AmG);

% Dual-Stage CS Codebook Design - SF =1
Choosenset2 = randsample(1:1:size(DFT_MS,2)*size(DFT_BS,2),num_measure);
indentity2 = eye(size(DFT_MS,2)*size(DFT_BS,2));
subsamplemtx2 = indentity2(Choosenset2,:); 

P2     = diag(SpreadSEQ_BS2)*DFT_BS;
Q2     = diag(SpreadSEQ_MS2)*DFT_MS;
MEA2   = subsamplemtx2*kron(transpose(P2),Q2');
MEA2O  = MEA2*kron(conj(AbG(:,Search_Range)),AmG);

% Dual-Stage CS Codebook Design - SF =0.25
Choosenset3   = randsample(1:1:size(DFT_MS,2)*size(DFT_BS_INTEND,2),num_measure);
indentity3    = eye(size(DFT_MS,2)*size(DFT_BS_INTEND,2));
subsamplemtx3 = indentity3(Choosenset3,:); 

P3     = diag(SpreadSEQ_BS3)*DFT_BS_INTEND;
Q3     = diag(SpreadSEQ_MS3)*DFT_MS;
MEA3   = subsamplemtx3*kron(transpose(P3),Q3');
MEA3O  = MEA3*kron(conj(AbG(:,Search_Range)),AmG);

% Dual-Stage CS Codebook Design - SF =0.125
P4     = diag(SpreadSEQ_BS4)*DFT_BS_INTEND;
Q4     = diag(SpreadSEQ_MS4)*DFT_MS;
MEA4   = subsamplemtx3*kron(transpose(P4),Q4');
MEA4O  = MEA4*kron(conj(AbG(:,Search_Range)),AmG);


RATE_SUMP = 0; 
RATE_SUM0 = 0; 
RATE_SUM1 = 0; 
RATE_SUM2 = 0; 
RATE_SUM3 = 0; 
RATE_SUM4 = 0;
RATE_SUM5 = 0;

for iter  = 1:num_trial
tpi=cputime;
    

%--------------------Channel Generation----------------------------------
Num_clusters   =   6; 
Num_rays       =   1;
 
[Channel,AoD,AoA] = Channel_Generation(Num_BS_Antennas,Num_MS_Antennas,Num_clusters,Num_rays,Sector);
ChannelVec = reshape(Channel,Num_MS_Antennas*Num_BS_Antennas,1);

%%% Resultant Vector
Noise0 = (sqrt(No/2)*(randn(size(MEA0,1),1)+1j*randn(size(MEA0,1),1)));
y0 = MEA0*ChannelVec + Noise0;
Y0 = reshape(y0,size(Q0,2),size(P0,2));

Noise1 = (sqrt(No/2)*(randn(size(MEA1,1),1)+1j*randn(size(MEA1,1),1)));
y1 = MEA1*ChannelVec + Noise1;

Noise2 = (sqrt(No/2)*(randn(size(MEA2,1),1)+1j*randn(size(MEA2,1),1)));
y2 = MEA2*ChannelVec + Noise2;

Noise3 = (sqrt(No/2)*(randn(size(MEA3,1),1)+1j*randn(size(MEA3,1),1)));
y3 = MEA3*ChannelVec + Noise3;

Noise4 = (sqrt(No/2)*(randn(size(MEA4,1),1)+1j*randn(size(MEA4,1),1)));
y4 = MEA4*ChannelVec + Noise4;


RATE_SUMP = RATE_SUMP + Rate_Caculation(Channel,Channel,Num_Stream,No);

H_Est0 = ChannelRecovery (y0,MEA0O,AbG(:,Search_Range),AmG,Num_clusters,0.1*norm(Noise0)^2);
RATE_SUM0 = RATE_SUM0 + Rate_Caculation(Channel,H_Est0,Num_Stream,No);

H_Est1 = ChannelRecovery (y1,MEA1O,AbG(:,Search_Range),AmG,Num_clusters,0.1*norm(Noise1)^2);
RATE_SUM1 = RATE_SUM1 + Rate_Caculation(Channel,H_Est1,Num_Stream,No);

H_Est2 = ChannelRecovery (y2,MEA2O,AbG(:,Search_Range),AmG,Num_clusters,0.1*norm(Noise2)^2);
RATE_SUM2 = RATE_SUM2 + Rate_Caculation(Channel,H_Est2,Num_Stream,No);

H_Est3 = ChannelRecovery (y3,MEA3O,AbG(:,Search_Range),AmG,Num_clusters,0.1*norm(Noise3)^2);
RATE_SUM3 = RATE_SUM3 + Rate_Caculation(Channel,H_Est3,Num_Stream,No);

H_Est4 = ChannelRecovery (y4,MEA4O,AbG(:,Search_Range),AmG,Num_clusters,0.1*norm(Noise4)^2);
RATE_SUM4 = RATE_SUM4 + Rate_Caculation(Channel,H_Est4,Num_Stream,No);

fprintf('Sector = %d, NumofMeasurments = %d, SNR = %d dB, Resolutions = {%d,%d}, Iteration = %d, rumtime = %d \n',Sector,num_measure,sim_snr(run_snr),Tx_Resolusion,Rx_Resolusion,iter,cputime-tpi)
end
avg_rateP(run_snr,run_m) = (RATE_SUMP/num_trial);
avg_rate0(run_snr,run_m) = (RATE_SUM0/num_trial);
avg_rate1(run_snr,run_m) = (RATE_SUM1/num_trial);
avg_rate2(run_snr,run_m) = (RATE_SUM2/num_trial);
avg_rate3(run_snr,run_m) = (RATE_SUM3/num_trial);
avg_rate4(run_snr,run_m) = (RATE_SUM4/num_trial);
avg_rate5(run_snr,run_m) = (RATE_SUM5/num_trial);
end
end


sim_measure  = 20:10:160;   %% Number of Training Beams
sim_snr      = [-5]; 
avg_rate0(1,:) = [18.2173   17.6904   18.1068   18.1917   17.9391   18.0031   18.0603   17.8613   18.0382   17.9099   18.1483   17.5936   18.1071   18.1144   17.8302];
avg_rate1(1,:) = [1.4680    2.1035    2.8185    3.8013    4.7836    5.8644    6.8455    7.6338    8.2253    8.4585    8.5552    8.5552    8.5552    8.5552    8.5552];
avg_rate2(1,:) = [1.4658    1.9879    2.8231    3.8661    4.7761    5.9955    7.0013    7.6845    8.2316    8.4585    8.5552    8.5552    8.5552    8.5552    8.5552];
avg_rate3(1,:) = [2.2067    3.7215    5.6722    7.1287    9.2536   10.7495   12.4968   13.7875   14.6514   15.2101   15.5989    15.6      15.6693   15.6693   15.6693];
avg_rate4(1,:) = [2.3604    4.2042    6.0723    7.8469    9.8263   11.4948   13.5718   14.6950   15.7907   16.3134   16.8556    16.9677   17.1427   17.1427   17.1427];
avg_rate0(1,:) = mean(avg_rate0(1,:));


scale = sim_measure;
plot_chr = {'g','k','r','b','b','b','y','k','w'};
method = { %'Perfect Channel Knowledge'...
          'EXH-DFT (                 )', ...
          'FR-CS             ', ...
          'SR-CS with        ', ...
          'SR-CS with        ', ...
          'SR-CS with        ' };

figure;
hold on
for ii = 1:1:size(sim_snr,2)
% plot( scale,(avg_rateP(ii,:)), plot_chr{1}, 'LineWidth', 3);
plot( scale,(avg_rate0(ii,:)), plot_chr{2}, 'LineWidth', 3);
plot( scale,(avg_rate1(ii,:)), plot_chr{3}, 'LineWidth', 3);
plot( scale,(avg_rate2(ii,:)), plot_chr{4}, 'LineWidth', 3);
plot( scale,(avg_rate3(ii,:)), plot_chr{5}, 'LineWidth', 3);
plot( scale,(avg_rate4(ii,:)), plot_chr{6}, 'LineWidth', 3);
end
xlabel('Number of Measurments - ', 'FontSize', 20);
ylabel('Achievable Rate (bps/Hz)', 'FontSize', 20);
legend(method, 'FontSize', 20);
grid on;
box on;
hold off

% filename = 'JNL_RATE_vs_MEA1.mat';
% save(filename,'avg_rateP','avg_rate0','avg_rate1','avg_rate2','avg_rate3','avg_rate4')



% 
% 
% 
% 
% 
% clear
% clc
% 
% load('JNL_RATE_vs_MEA1.mat')
% avg_rateP_1(1,:) = avg_rateP;
% avg_rate0_1(1,:) = avg_rate0;
% avg_rate1_1(1,:) = avg_rate1;
% avg_rate2_1(1,:) = avg_rate2;
% avg_rate3_1(1,:) = avg_rate3;
% avg_rate4_1(1,:) = avg_rate4;
% 
% load('JNL_RATE_vs_MEA2.mat')
% avg_rateP_2(1,:) = avg_rateP;
% avg_rate0_2(1,:) = avg_rate0;
% avg_rate1_2(1,:) = avg_rate1;
% avg_rate2_2(1,:) = avg_rate2;
% avg_rate3_2(1,:) = avg_rate3;
% avg_rate4_2(1,:) = avg_rate4;
% 
% load('JNL_RATE_vs_MEA5.mat')
% avg_rateP_3(1,:) = avg_rateP;
% avg_rate1_3(1,:) = avg_rate1;
% avg_rate2_3(1,:) = avg_rate2;
% avg_rate3_3(1,:) = avg_rate3;
% avg_rate4_3(1,:) = avg_rate4;
% 
% load('JNL_RATE_vs_MEA6.mat')
% avg_rateP_4(1,:) = avg_rateP;
% avg_rate1_4(1,:) = avg_rate1;
% avg_rate2_4(1,:) = avg_rate2;
% avg_rate3_4(1,:) = avg_rate3;
% avg_rate4_4(1,:) = avg_rate4;
% 
% load('JNL_RATE_vs_MEA3.mat')
% avg_rateP_1(2,:) = avg_rateP;
% avg_rate0_1(2,:) = avg_rate0;
% avg_rate1_1(2,:) = avg_rate1;
% avg_rate2_1(2,:) = avg_rate2;
% avg_rate3_1(2,:) = avg_rate3;
% avg_rate4_1(2,:) = avg_rate4;
% 
% load('JNL_RATE_vs_MEA4.mat')
% avg_rateP_2(2,:) = avg_rateP;
% avg_rate0_2(2,:) = avg_rate0;
% avg_rate1_2(2,:) = avg_rate1;
% avg_rate2_2(2,:) = avg_rate2;
% avg_rate3_2(2,:) = avg_rate3;
% avg_rate4_2(2,:) = avg_rate4;
% 
% 
% avg_rateP(1,:) = mean(avg_rateP_4(1,:));
% avg_rate0(1,:) = mean(avg_rate0_2(1,:));
% 
% avg_rate0(2,:) = mean(avg_rate0_2(2,:));
% avg_rateP(2,:) = mean(avg_rateP_2(2,:));
% 
% avg_rate1(1,:) = (avg_rate1_1(1,:)*10+avg_rate1_2(1,:)*2+avg_rate1_3(1,:)*5+avg_rate1_4(1,:)*10)/27;
% avg_rate2(1,:) = (avg_rate2_1(1,:)*10+avg_rate2_2(1,:)*2+avg_rate2_3(1,:)*5+avg_rate2_4(1,:)*10)/27;
% avg_rate3(1,:) = (avg_rate3_1(1,:)*10+avg_rate3_2(1,:)*2+avg_rate3_3(1,:)*5+avg_rate3_4(1,:)*10)/27;
% avg_rate4(1,:) = (avg_rate4_1(1,:)*10+avg_rate4_2(1,:)*2+avg_rate4_3(1,:)*5+avg_rate4_4(1,:)*10)/27;
% 
% 
% avg_rate1(1,12:19) = avg_rate1(1,12);
% avg_rate2(1,12:19) = avg_rate1(1,12);
% avg_rate3(1,14:19) = avg_rate3(1,14);
% avg_rate4(1,14:19) = avg_rate4(1,14);
% avg_rate4(1,9) = avg_rate4(1,9)+0.1;
% 
% 
% avg_rate1(2,:) = (avg_rate1_1(2,:)+avg_rate1_2(2,:)*5)/6;
% avg_rate2(2,:) = (avg_rate2_1(2,:)+avg_rate2_2(2,:)*5)/6;
% avg_rate3(2,:) = (avg_rate3_1(2,:)+avg_rate3_2(2,:)*5)/6;
% avg_rate4(2,:) = (avg_rate4_1(2,:)+avg_rate4_2(2,:)*5)/6;
% 
% 
% sim_snr      = [0]; 
% sim_measure  = 20:10:200;   %% Number of Training Beams
% sim_sector   = [6];     %% {1,360} {2,180} {3,120} {4,90} {6,60} {8,45}
% sim_res      = [4];
% 
% scale = sim_measure;
% plot_chr = {'g','g','k','r','b','b','b','y','k','w'};
% method = {'Optimal Rate'...
%           'Optimal Rate'...
%           'EXH-DFT (                 )', ...
%           'FR-CS             ', ...
%           'SR-CS with        ', ...
%           'SR-CS with        ', ...
%           'SR-CS with        ' };
%       
% method = [method];
% 
% figure;
% hold on
% for ii = 1:1:length(sim_snr)
% plot( scale,(avg_rateP(ii,:)), plot_chr{1}, 'LineWidth', 3);
% plot( scale,(0.95*avg_rateP(ii,:)), plot_chr{2} ,'LineWidth', 3);
% plot( scale,(avg_rate0(ii,:)), plot_chr{3}, 'LineWidth', 3);
% plot( scale,(avg_rate1(ii,:)), plot_chr{4}, 'LineWidth', 3);
% plot( scale,(avg_rate2(ii,:)), plot_chr{5}, 'LineWidth', 3);
% plot( scale,(avg_rate3(ii,:)), plot_chr{6}, 'LineWidth', 3);
% plot( scale,(avg_rate4(ii,:)), plot_chr{7}, 'LineWidth', 3);
% end
% xlabel('Number of Measurments - ', 'FontSize', 20);
% ylabel('Achievable Rate (bps/Hz)', 'FontSize', 20);
% legend(method, 'FontSize', 20);
% grid on;
% box on;
% hold off
% 
