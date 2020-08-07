%%% Written by Van-Dinh Nguyen
%%% Authors: Van-Dinh Nguyen, Hoang Duong Tuan, Trung Q. Duong, Oh-Soon Shin, and H. Vincent Poor
%%% Corresponding to the paper: "Joint Fractional Time Allocation and Beamforming for Downlink Multiuser MISO Systems," 
% % % IEEE Communications Letters, vol. 21, no. 12, pp. 2650-2653, Dec. 2017
%%%  


close all
clear
clc

addpath(genpath('./yalmip'));
addpath(genpath('./SeDuMi_1_3'));
addpath(genpath('./SDPT3-4.0'));

tic

%% Initialization

%load('test.mat')
load('C:\Users\WCL-SON\Dropbox\Matlab CODE\FractionalTime\Sum Rate\LayoutParameters\D_Hk.mat')
load('C:\Users\WCL-SON\Dropbox\Matlab CODE\FractionalTime\Sum Rate\LayoutParameters\D_Gm.mat')
load('C:\Users\WCL-SON\Dropbox\Matlab CODE\FractionalTime\Sum Rate\LayoutParameters\positionDownlinkUsers.mat')


NumberOfRunning = 1;

Conven = 0;
Fixed_Timegroup = 0;

MaxIteration = 40;

N_t = 8; % number of antennas at BS

G = 2; % 2 regions - inner and outer zone

K = 4; % number of users in inner zone
M = 4; % number of users in outer zone

Rate_Threshold = 0.5;

Pbs_dB = 42; % dowlink power [dBm]

 Pbs = 10^(Pbs_dB/10); % downlink power


Noise = -174; % dBm/Hz
BW = 10; % MHz

sigma_k = 10^(Noise/10)*BW*10^6;

sigma_m = 10^(Noise/10)*BW*10^6;

RadiusOfCell = 500;
InnerZoneRadius = 200;
RadiusOfNearestUser = 10;
StandardDeviation = 8;
ploss = 3;

Parameters = [RadiusOfCell InnerZoneRadius RadiusOfNearestUser StandardDeviation ploss];

   %[D_Hk, D_Gm, positionDownlinkUsers] = CreateLargeScaleFading( K(end), M(end), Parameters);

Layout = Plot_Layout( RadiusOfCell, InnerZoneRadius, positionDownlinkUsers, K(end) );

if (length(K)>1)
    XDir = K;
elseif (length(Rate_Threshold)>1)
    XDir = Rate_Threshold;
elseif (length(Pbs)>1)
    XDir = Pbs;
else
    XDir = 1;
    MaxIteration = 2*MaxIteration;
end

Alg1_Rate_all_zone1 = zeros(NumberOfRunning, length(XDir));
Alg1_Rate_all_zone2 = zeros(NumberOfRunning, length(XDir));

Alg1_Rate_all = zeros(NumberOfRunning, length(XDir));
Alg2_Rate_all = zeros(NumberOfRunning, length(XDir));
Alg3_Rate_all = zeros(NumberOfRunning, length(XDir));
Alg4_Rate_all = zeros(NumberOfRunning, length(XDir));
Alg5_Rate_all = zeros(NumberOfRunning, length(XDir));

DownlinkRate_PerGroupPerUser = zeros(K, G, NumberOfRunning);


%h = waitbar(0, 'Please wait ...','Name','Percent of completion');

% Clus = parcluster('local');
% Clus.NumWorkers = 4;

% poolobj = parpool(Clus, Clus.NumWorkers);

i_NumOfSim = 0;

while (i_NumOfSim<NumberOfRunning)
    
    i_NumOfSim = i_NumOfSim + 1
    
    % Channels
    if (length(K)==1)
        Hk_init = CreateChannel(1, 1, N_t, K)*D_Hk;
        Gm_init = CreateChannel(1, 1, N_t, M)*D_Gm;
        if (Conven)
            Hk = [Hk_init Gm_init];
            Gm = [];
        else
            Hk = Hk_init;
            Gm = Gm_init;
        end
    end
    
    
%     percent = (i_NumOfSim-1) / (NumberOfRunning);
%     waitbar(percent, h,[ num2str(floor(percent*100)) ' % Completed'])


for iXDir = XDir
% parfor i = 1:1:length(rho_dB)

    iK = K(1);
    iRateThreshold = Rate_Threshold(1);
    iPbs = Pbs(1);
    
    if (length(K)>1)
        iK = iXDir;
         Hk_init = CreateChannel(1, 1, N_t, iK)*D_Hk(1:iK, 1:iK);
         Gm_init = CreateChannel(1, 1, N_t, iK)*D_Gm(1:iK, 1:iK);
        if (Conven)
            Hk = [Hk Gm];
            Gm = [];
        else
            Hk = Hk_init;
            Gm = Gm_init;
        end
    elseif (length(Rate_Threshold)>1)
        iRateThreshold = iXDir;
    elseif (length(Pbs)>1)
        iPbs = iXDir;
    end
    
    size(Hk)
    size(Gm)
 
    
    %% Fractional Time


    
      [isBreak1,  global_OptValue1,global_OptValue2, Alg1_OptValueChain, DownlinkRate_PerGroupPerUser1, DownlinkRate_PerGroupPerUser2 ]  = ProposedAlg( iPbs, Hk, Gm, iRateThreshold, MaxIteration, Conven, Fixed_Timegroup);

%     [isBreak1, Alg1_OptValue, Alg1_OptValueChain, DownlinkRate_PerGroupPerUser1, DownlinkRate_PerGroupPerUser2 ]  = ProposedAlg( iPbs, Hk, Gm, iRateThreshold, MaxIteration, Conven, Fixed_Timegroup);

%     [isBreak1, Alg2_OptValue, Alg2_OptValueChain, DownlinkRate_PerGroupPerUser1, DownlinkRate_PerGroupPerUser2 ]  = ProposedAlgNOMA1( iPbs, Hk, Gm, iRateThreshold, MaxIteration, Conven, Fixed_Timegroup); 
% 
%      [isBreak1, Alg3_OptValue, Alg3_OptValueChain, DownlinkRate_PerGroupPerUser1, DownlinkRate_PerGroupPerUser2 ]  = ProposedAlgNOMAperzone( iPbs, Hk, Gm, iRateThreshold, MaxIteration, Conven, Fixed_Timegroup);
%     
%       [isBreak1, Alg4_OptValue, Alg4_OptValueChain, DownlinkRate_PerGroupPerUser1, DownlinkRate_PerGroupPerUser2 ]  = ProposedAlgNOMAzone2( iPbs, Hk, Gm, iRateThreshold, MaxIteration, Conven, Fixed_Timegroup);
% 
%      [isBreak1, Alg5_OptValue, Alg5_OptValueChain, DownlinkRate_PerGroupPerUser1, DownlinkRate_PerGroupPerUser2 ]  = ProposedAlgNOMAzone1( iPbs, Hk, Gm, iRateThreshold, MaxIteration, Conven, Fixed_Timegroup);
  
%      Alg1_Rate_all(i_NumOfSim, iXDir) = Alg1_OptValue;
%      Alg2_Rate_all(i_NumOfSim, iXDir) = Alg2_OptValue;
%      Alg3_Rate_all(i_NumOfSim, iXDir) = Alg3_OptValue;
%      Alg4_Rate_all(i_NumOfSim, iXDir) = Alg4_OptValue; 
%      Alg5_Rate_all(i_NumOfSim, iXDir) = Alg5_OptValue;
    
       Alg1_Rate_all_zone1(i_NumOfSim, iXDir) = global_OptValue1;
       Alg1_Rate_all_zone2(i_NumOfSim, iXDir) = global_OptValue2;
    
end
%     waitbar(i_NumOfSim / NumberOfRunning, h,[ str_Group ':  ' num2str(floor(i_NumOfSim*100 / NumberOfRunning)) ' % Completed'])
    
    



end



if (NumberOfRunning == 1)
   Alg1_Ratezone1 = Alg1_Rate_all_zone1;
   Alg1_Ratezone2 = Alg1_Rate_all_zone2;
%     Alg1_Rate = Alg1_Rate_all
%     Alg2_Rate = Alg2_Rate_all
%     Alg3_Rate = Alg3_Rate_all
%     Alg4_Rate = Alg4_Rate_all
%     Alg5_Rate = Alg5_Rate_all
else
    Alg1_Ratezone1 = mean(Alg1_Rate_all_zone1);
    Alg1_Ratezone2 = mean(Alg1_Rate_all_zone2);
%     Alg1_Rate = mean(Alg1_Rate_all)
%     Alg2_Rate = mean(Alg2_Rate_all)
%     Alg3_Rate = mean(Alg3_Rate_all)
%     Alg4_Rate = mean(Alg4_Rate_all)
%     Alg5_Rate = mean(Alg5_Rate_all)
end


%close(h);

figure;

hold on

if (XDir==1)
%     plot([1:1:length(Alg1_OptValueChain)], Alg1_OptValueChain, 'b-', 'linewidth', 2, 'markersize',9);hold on
%     plot([1:1:length(Alg2_OptValueChain)], Alg2_OptValueChain, 'k-', 'linewidth', 2, 'markersize',9);hold on
%     plot([1:1:length(Alg3_OptValueChain)], Alg3_OptValueChain, 'r-', 'linewidth', 2, 'markersize',9);hold on
%     plot([1:1:length(Alg4_OptValueChain)], Alg4_OptValueChain, 'g-', 'linewidth', 2, 'markersize',9);hold on
%     plot([1:1:length(Alg5_OptValueChain)], Alg5_OptValueChain, 'c-', 'linewidth', 2, 'markersize',9);
else
    plot(XDir, Alg1_Rate, 'bs--', 'linewidth', 2, 'markersize',9);
end
hold on


time = toc
hour = 0; min=0;
if (time>3600)
    hour = floor(time/3600);
    time = mod(time,3600);
end
if (time>60)
    min = floor(time/60);
    time = mod(time,60);
end
disp(['Running Time = ' num2str(hour) ' h ' num2str(min) ' m ' num2str(time) ' s.' ]);