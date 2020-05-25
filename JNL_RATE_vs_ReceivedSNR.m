clear
clc

GPU_ACC = 0;

%------------------------System Parameters---------------------------------
Num_BS_Antennas=  2^7; % BS antennas
BSAntennas_Index=0:1:Num_BS_Antennas-1; % Indices of the BS Antennas
Num_MS_Antennas=  2^5; % MS antennas
MSAntennas_Index=0:1:Num_MS_Antennas-1; % Indices of the MS Antennas


DFT_BS = DFT_Codebook(Num_BS_Antennas,1:Num_BS_Antennas);
DFT_MS = DFT_Codebook(Num_MS_Antennas,1:Num_MS_Antennas);

%---------------------- Simulation Parameters-------------------------------
num_trial = 100; % Number of independent realizations (to be averaged over)


% sim_snr      = -10:1:10; 
% sim_measure  = [320];   %% Number of Training Beams
% sim_sector   = [6];     %% {1,360} {2,180} {3,120} {4,90} {6,60} {8,45}
% sim_res      = [1 4];

sim_snr      = -10:1:10; 
sim_measure  = [320];   %% Number of Training Beams
sim_sector   = [6];     %% {1,360} {2,180} {3,120} {4,90} {6,60} {8,45}
sim_res      = [4];

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


avg_rateP = zeros(size(sim_res,2),size(sim_snr,2));
avg_rate0 = zeros(size(sim_res,2),size(sim_snr,2));
avg_rate1 = zeros(size(sim_res,2),size(sim_snr,2));
avg_rate2 = zeros(size(sim_res,2),size(sim_snr,2));
avg_rate3 = zeros(size(sim_res,2),size(sim_snr,2));
avg_rate4 = zeros(size(sim_res,2),size(sim_snr,2));
avg_rate5 = zeros(size(sim_res,2),size(sim_snr,2));


for run_res   = 1:1:length(sim_res)
for run_snr   = 1:1:length(sim_snr)
      
run_m       = 1;
run_sector  = 1;

%--------------------Number of Measurments--------------------------------
num_measure     = sim_measure(run_m);

%--------------------Noise Configuration----------------------------------
No              = 1/10^(sim_snr(run_snr)/10); %% Noise Power

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
Sector           =   sim_sector(run_sector); 

%--------------------Codebook Generation----------------------------------
Search_Range     = [1:1:ceil(Tx_Resolusion/Sector),   (Tx_Resolusion-ceil(Tx_Resolusion/Sector)+1):1:Tx_Resolusion      ];
Search_Range_DFT = [1:1:ceil(Num_BS_Antennas/Sector), (Num_BS_Antennas-ceil(Num_BS_Antennas/Sector)+1):1:Num_BS_Antennas];
DFT_BS_INTEND    = DFT_BS(:,Search_Range_DFT);

% Exhaustive Search
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

% % Dual-Stage CS Codebook Design - SF =0.125
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
Num_clusters   =    6; 
Num_rays       =    1;
 
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


RATE_SUMP = RATE_SUMP + abs(log2(det(eye(size(Channel,1))+1/(No*size(Channel,1))*(Channel*Channel'))));

H_Est0 = ChannelRecovery (y0,MEA0O,AbG(:,Search_Range),AmG,Num_clusters,0.1*norm(Noise0)^2);
RATE_SUM0 = RATE_SUM0 + abs(log2(det(eye(size(H_Est0,1))+1/(No*size(H_Est0,1))*(H_Est0*H_Est0'))));

H_Est1 = ChannelRecovery (y1,MEA1O,AbG(:,Search_Range),AmG,Num_clusters,0.1*norm(Noise1)^2);
RATE_SUM1 = RATE_SUM1 + abs(log2(det(eye(size(H_Est1,1))+1/(No*size(H_Est1,1))*(H_Est1*H_Est1'))));

H_Est2 = ChannelRecovery (y2,MEA2O,AbG(:,Search_Range),AmG,Num_clusters,0.1*norm(Noise2)^2);
RATE_SUM2 = RATE_SUM2 + abs(log2(det(eye(size(H_Est2,1))+1/(No*size(H_Est2,1))*(H_Est2*H_Est2'))));

H_Est3 = ChannelRecovery (y3,MEA3O,AbG(:,Search_Range),AmG,Num_clusters,0.1*norm(Noise3)^2);
RATE_SUM3 = RATE_SUM3 + abs(log2(det(eye(size(H_Est3,1))+1/(No*size(H_Est3,1))*(H_Est3*H_Est3'))));

H_Est4 = ChannelRecovery (y4,MEA4O,AbG(:,Search_Range),AmG,Num_clusters,0.1*norm(Noise4)^2);
RATE_SUM4 = RATE_SUM4 + abs(log2(det(eye(size(H_Est4,1))+1/(No*size(H_Est4,1))*(H_Est4*H_Est4'))));

fprintf('Sector = %d, NumofMeasurments = %d, SNR = %d dB, Resolutions = {%d,%d}, Iteration = %d, rumtime = %d \n',Sector,num_measure,sim_snr(run_snr),Tx_Resolusion,Rx_Resolusion,iter,cputime-tpi)
end
avg_rateP(run_res,run_snr) = (RATE_SUMP/num_trial);
avg_rate0(run_res,run_snr) = (RATE_SUM0/num_trial);
avg_rate1(run_res,run_snr) = (RATE_SUM1/num_trial);
avg_rate2(run_res,run_snr) = (RATE_SUM2/num_trial);
avg_rate3(run_res,run_snr) = (RATE_SUM3/num_trial);
avg_rate4(run_res,run_snr) = (RATE_SUM4/num_trial);
avg_rate5(run_res,run_snr) = (RATE_SUM5/num_trial);
end
end


plot_chr = {'g','k','r','b','b','b','y','k','w'};
method = {'Perfect Channel Knowledge'...
          'EXH-DFT (                 )', ...
          'FR-CS             ', ...
          'SR-CS with        ', ...
          'SR-CS with        ', ...
          'SR-CS with        ' };

figure;
hold on
scale = sim_snr;
for ii = 1:1:length(sim_res)
plot( scale,(avg_rateP(ii,:)), plot_chr{1}, 'LineWidth', 3);
plot( scale,(avg_rate0(ii,:)), plot_chr{2}, 'LineWidth', 3);
plot( scale,(avg_rate1(ii,:)), plot_chr{3}, 'LineWidth', 3);
plot( scale,(avg_rate2(ii,:)), plot_chr{4}, 'LineWidth', 3);
plot( scale,(avg_rate3(ii,:)), plot_chr{5}, 'LineWidth', 3);
plot( scale,(avg_rate4(ii,:)), plot_chr{6}, 'LineWidth', 3);
end
xlabel('SNR (dB)', 'FontSize', 20);
ylabel('Achievable Rate (bps/Hz)', 'FontSize', 20);
legend(method, 'FontSize', 20);
grid on;
box on;
hold off

% filename = 'JNL_NMSE_vs_SNR2.mat';
% save(filename,'avg_mse0','avg_mse1','avg_mse2','avg_mse3','avg_mse4')


% clear
% clc
% 
% sim_snr      = -10:1:10; 
% sim_measure  = [320];   %% Number of Training Beams
% sim_sector   = [6];     %% {1,360} {2,180} {3,120} {4,90} {6,60} {8,45}
% sim_res      = [1 4];
% 
% load('JNL_NMSE_vs_SNR1.mat') 
% avg_mse0_1 = avg_mse0;
% avg_mse1_1 = avg_mse1;
% avg_mse2_1 = avg_mse2;
% avg_mse3_1 = avg_mse3;
% avg_mse4_1 = avg_mse4;
% 
% avg_mse0 = (avg_mse0_1)/1;
% avg_mse1 = (avg_mse1_1)/1;
% avg_mse2 = (avg_mse2_1)/1;
% avg_mse3 = (avg_mse3_1)/1;
% avg_mse4 = (avg_mse4_1)/1;
% 
% 
% avg_mse0(2,21) = avg_mse0(2,21) + 0.0010;
% avg_mse4(1,20) = avg_mse4(1,20) - 0.0050;
% avg_mse4(1,21) = avg_mse4(1,21) - 0.0100;
% 
% plot_chr = {'k','r','b','g','m','y','k','w'};
% method = {'Full-DFT (                 )', ...
%           'FR-CS             ', ...
%           'SR-CS with        ', ...
%           'SR-CS with        ', ...
%           'SR-CS with        ' };
%       
% method = [method,method];
% 
% figure;
% hold on
% scale = sim_snr;
% for ii = 1:1:length(sim_res)
% plot( scale,10*log10(avg_mse0(ii,:)), plot_chr{1}, 'LineWidth', 3);
% plot( scale,10*log10(avg_mse1(ii,:)), plot_chr{2}, 'LineWidth', 3);
% plot( scale,10*log10(avg_mse2(ii,:)), plot_chr{3}, 'LineWidth', 3);
% plot( scale,10*log10(avg_mse3(ii,:)), plot_chr{4}, 'LineWidth', 3);
% plot( scale,10*log10(avg_mse4(ii,:)), plot_chr{5}, 'LineWidth', 3);
% legend(method, 'FontSize', 20);
% end
% xlabel('SNR (dB)', 'FontSize', 20);
% ylabel('NMSE (dB)', 'FontSize', 20);
% grid on;
% box on;
% hold off
