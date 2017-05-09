clear;
clc;
close all;

rng(0);
set(0,'defaultlinelinewidth',1.5)
Sparsity_Indx = 10;

%% Simulation Paramters %%
TxScenario                          =   3;     % 1- Contigious Single CC
% 2- Multi-Cluster Single CC
% 3- Intra-Band CA Two CC
N_symbols                           =  800;    % No of OFDM Symbols for Simulation
TxSignalType                        =   2;     % 1- OFDM for DL Transmission
% 2- SCFDMA for UL Transmission
ModulationType                      =   3;     % 1- QPSK
% 2- 16-QAM
% 3- 64-QAM
POWER_PLOT_1MHZ                     =   0;     % Use 1 MHz power plot

% RF and PA paramters
Pout_dBm                            =   23;    % Desired Tx Output Power after the PA in dBm

MemoryLessPA                        =   1;     % 1- MemoryLess PA model
% 0- PH PA model with memory
MemoryLessDPD                       =   1;     % 1- MemoryLess DPD
% 0- MemoryDPD

IM3_right   = 1;      %Supress RHS
IM3_left    = 0;      %Supress LHS NOT FUNCTIONAL YET

IM5_right   = 0; %NOT FUNCTIONAL YET
IM5_left    = 0; %NOT FUNCTIONAL YET

IM3_order_3 = 0;
IM3_order_5 = 0;
IM3_order_7 = 0;
IM3_order_9 = 1;

USE_WARP = 1;                %Use a WARP Board
DO_Training = 1;             %Train coefficients or run code with predefined coefficients.
DPD_LearningRate = 1;         % Decorrelating DPD Learning Rate
DPD_LearningBlockSize  = 528;    % Decorrelating DPD Learning Block size
DPD_FilteringBlockSize = 1024;    % Decorrelating DPD Filtering Block size
NumLearningSamples = 0.5*200000;


% Gain 3,45
%  IM3_ThirdOrder_Coeffs   =    6.4367 - 1.6282i;%39.2965 +10.5928i;
%  IM3_FifthOrder_Coeffs   =    7.5791 - 0.6598i;%1.5774 + 3.5535i;
%  IM3_SeventhOrder_Coeffs =        -4.4727 + 0.6793i;%  -7.2355 - 2.0117i;
%  IM3_NinthOrder_Coeffs   =  0;%3.8273 + 0.1256i;

%Gain 3,55
% IM3_ThirdOrder_Coeffs   =    40.6080 + 1.2534i;%39.2965 +10.5928i;
% IM3_FifthOrder_Coeffs   =    -4.7023 + 3.2298i;%1.5774 + 3.5535i;
% IM3_SeventhOrder_Coeffs =    -7.1796 - 0.9978i;%  -7.2355 - 2.0117i;
% IM3_NinthOrder_Coeffs   =  0;%3.8273 + 0.1256i;

%Interpolate for 3,50
 IM3_ThirdOrder_Coeffs   =  23.5223 - 0.1874i;%39.2965 +10.5928i;
 IM3_FifthOrder_Coeffs   =  1.4384 + 1.2850i;%1.5774 + 3.5535i;
 IM3_SeventhOrder_Coeffs =   -5.8262 - 0.1593i;%  -7.2355 - 2.0117i;
 IM3_NinthOrder_Coeffs   =  0;%3.8273 + 0.1256i;


IM3_ThirdOrder_Coeffs_start   = IM3_ThirdOrder_Coeffs;
IM3_FifthOrder_Coeffs_start   = IM3_FifthOrder_Coeffs;
IM3_SeventhOrder_Coeffs_start = IM3_SeventhOrder_Coeffs;
IM3_NinthOrder_Coeffs_start   = IM3_NinthOrder_Coeffs;

WARP_ID = 00031;

if(USE_WARP)
   Max_input_scale_factor = 0.8;   %Scale the input between 0 and 1
   Dont_trust_cyclosync = 1;
   
   % RX variables
   Channel_RX = 6;
   RxGainRF   = 2;  % Rx RF Gain in [1:3] (ignored if USE_AGC is true)
   RxGainBB   = 14; % Rx Baseband Gain in [0:31] (ignored if USE_AGC is true)
   RX_LPF     = 2;  % [0,1,2,3] for approx ![7.5,9.5,14,18]MHz corner
   
   % TX variables
   Channel_TX = 6;
   TxGainBB   = 3;  % [0,1,2,3] for approx ![-5, -3, -1.5, 0]dB baseband gain
   TXGainRF   = 50; % [0:63] for approx [0:31]dB RF gain
   TX_LPF     = 3;  % [1,2,3] for approx [12,18,24]MHz corner frequencies ([24,36,48]MHz bandwidths)
   
   global Spacing phase USE_WARP;
   Spacing = (Channel_RX - Channel_TX) * (5*10^6); %5 MHz per channel
end

%% Baseband Tx signal paramters
LTE_Bandwidth = 1.4; % Carrier BW of the LTE Tx Signal
if LTE_Bandwidth == 1.4
   NRB1 = 6;            % Number of RB allocated to first CC
   NRB2 = 6;            % Number of RB allocated to second CC
end
if LTE_Bandwidth == 3
   NRB1 = 15;            % Number of RB allocated to first CC
   NRB2 = 15;            % Number of RB allocated to second CC
end
if LTE_Bandwidth == 5
   NRB1 = 25;            % Number of RB allocated to first CC
   NRB2 = 25;            % Number of RB allocated to second CC
end
if LTE_Bandwidth == 10
   NRB1 = 50;            % Number of RB allocated to first CC
   NRB2 = 50;            % Number of RB allocated to second CC
end

CarrierSpacing = 6;%20/3; % Spacing in MHz between the 2 CC

IM3_Freq = 3*(CarrierSpacing/2);
IM5_Freq = 5*(CarrierSpacing/2);
IM7_Freq = 7*(CarrierSpacing/2);
IM9_Freq = 9*(CarrierSpacing/2);
Signal_Bandwidth = NRB1*0.18;


%% Power Amplifier Model
% Estimated PA parameters at 26 dBm output power, 1.88 GHz, using a 10 MHz signal with 120 MHz Fs
PA_Power_Measured = 26;
MemorylessPA_Paramters = [0.9512 - 0.0946i;
   0.0239 + 0.1632i;
   0.0082 - 0.0727i;
   -0.0016 + 0.0147i;
   -0.0001 - 0.0011i];
% Memory depth = 4 per non-lineaity order
Memory_depth = 4;
MemoryPA_Paramters =  [0.5823 - 0.0608i;
   1.1417 - 0.1198i;
   -1.1964 + 0.1273i;
   0.4264 - 0.0407i;
   0.0074 + 0.1609i;
   0.0292 + 0.0037i;
   -0.0185 + 0.0002i;
   0.0032 - 0.0016i;
   0.0096 - 0.0727i;
   0.0034 + 0.0010i;
   -0.0060 - 0.0026i;
   0.0028 + 0.0012i;
   -0.0017 + 0.0149i;
   -0.0011 - 0.0008i;
   0.0014 + 0.0012i;
   -0.0007 - 0.0005i;
   -0.0001 - 0.0012i;
   0.0001 + 0.0001i;
   -0.0001 - 0.0001i;
   0.0001 + 0.0001i];



%% Baseband Equivilent LTE Signal Transmitter
[LTE_Signal, CC1, CC2, SystemFs, UpsamplingFactor] = LTE_Transmitter(LTE_Bandwidth,CarrierSpacing,NRB1,NRB2,...
   N_symbols,TxSignalType,ModulationType,TxScenario);

% Scale the Baseband generated signal to have a unit RMS power
if(USE_WARP)
   MAX_real = max(abs(real(LTE_Signal)));
   MAX_imag = max(abs(imag(LTE_Signal)));
   ScalingForPA = Max_input_scale_factor/max(MAX_imag,MAX_real);
else
   TX_PowerScale = sqrt(10^((PA_Power_Measured-Pout_dBm)/10));
   RMS_PAin = sqrt(mean(abs(LTE_Signal).^2));
   ScalingForPA = 1/(RMS_PAin*TX_PowerScale);
end

PA_InputSignal = LTE_Signal*ScalingForPA;
CC1 = CC1*ScalingForPA;
CC2 = CC2*ScalingForPA;

%% SETUP WARP
if(USE_WARP)
   global nodes node_tx ifc_ids node_rx eth_trig ts_tx DC SystemFs delay;
   nodes = wl_initNodes(1);
   node_tx = nodes(1);
   node_rx = nodes(1);
   eth_trig = wl_trigger_eth_udp_broadcast;
   wl_triggerManagerCmd(nodes, 'add_ethernet_trigger', [eth_trig]);
   trig_in_ids  = wl_getTriggerInputIDs(nodes(1));
   trig_out_ids = wl_getTriggerOutputIDs(nodes(1));
   wl_triggerManagerCmd(nodes, 'output_config_input_selection', [trig_out_ids.BASEBAND], [trig_in_ids.ETH_A]);
   ifc_ids = wl_getInterfaceIDs(nodes(1));
   wl_interfaceCmd(nodes, ifc_ids.RF_A, 'channel', 2.4, Channel_TX);
   wl_interfaceCmd(nodes, ifc_ids.RF_B, 'channel', 2.4, Channel_RX);
   wl_interfaceCmd(nodes, ifc_ids.RF_ALL, 'rx_gain_mode', 'manual');
   wl_interfaceCmd(nodes, ifc_ids.RF_ALL, 'rx_gains', RxGainRF, RxGainBB);
   wl_interfaceCmd(nodes, ifc_ids.RF_ALL, 'tx_gains', TxGainBB, TXGainRF);
   wl_interfaceCmd(nodes, ifc_ids.RF_ALL, 'tx_lpf_corn_freq',TX_LPF);
   wl_interfaceCmd(nodes, ifc_ids.RF_ALL, 'rx_lpf_corn_freq',RX_LPF );
   
   ts_tx   = 1 / (wl_basebandCmd(nodes(1), 'tx_buff_clk_freq'));
   ts_rx   = 1 / (wl_basebandCmd(nodes(1), 'rx_buff_clk_freq'));
   ts_rssi = 1 / (wl_basebandCmd(nodes(1), 'rx_rssi_clk_freq'));
   maximum_buffer_len = wl_basebandCmd(node_tx, ifc_ids.RF_A, 'tx_buff_max_num_samples');
   
   %    wl_basebandCmd(node_tx, 'continuous_tx', true);
   
   if(WARP_ID == 31)
      DC_real = 0.0210;
      DC_imag =  - 0.0110;
   end
   
   DC = DC_real + DC_imag*i;
   
end


%% Power Amplifier
if(USE_WARP)
   MAX_real = max(abs(real(PA_InputSignal)));
   MAX_imag = max(abs(imag(PA_InputSignal)));
   MAX = max(MAX_real,MAX_imag);
   display('MAX PAInput (just LTE) =');
   display(MAX);
   
   PA_OutputSignal = WARP_broadcast(PA_InputSignal);
   
   % Synchronization
   synchronization(Channel_RX, Channel_TX,PA_OutputSignal,PA_InputSignal);
   
   PA_OutputSignal = PA_OutputSignal(delay:delay+length(PA_InputSignal)-1);
   MAX_real = max(abs(real(PA_OutputSignal)));
   MAX_imag = max(abs(imag(PA_OutputSignal)));
   MAX = max(MAX_real,MAX_imag);
   display('MAX abs(PA_Output) =');
   display(MAX);
else
   if MemoryLessPA
      PA_OutputSignal = MemoryLess_PA(PA_InputSignal,MemorylessPA_Paramters);
   else
      PA_OutputSignal = MemoryPH_PA(PA_InputSignal,MemoryPA_Paramters,Memory_depth);
   end
end

global IM3_BasisThirdOrder IM3_BasisFifthOrder IM3_BasisSeventhOrder IM3_BasisNinthOrder;

%% First Generate the righthand 3rd basis function and do training
if(DO_Training)
   GenerateBasis(CC1,CC2,1,3,3);
   GenerateBasis(CC1,CC2,1,3,5);
   GenerateBasis(CC1,CC2,1,3,7);
   %GenerateBasis(CC1,CC2,1,3,9);
   
   [LoopDelay, FeedBackFilter] = DPD_LoopDelayEst(PA_InputSignal, IM3_BasisThirdOrder, MemoryLessPA, ...
      SystemFs,Signal_Bandwidth,IM3_Freq,MemorylessPA_Paramters,MemoryPA_Paramters,Memory_depth,0);
   
   %Apply initial guess
   %All_Coeffs = [IM3_ThirdOrder_Coeffs; IM3_FifthOrder_Coeffs; IM3_SeventhOrder_Coeffs; IM3_NinthOrder_Coeffs];
   AdaptiveFilterDelay  = length(IM3_ThirdOrder_Coeffs) - 1;
   IM3_DPD_Signal = conv(IM3_BasisThirdOrder,IM3_ThirdOrder_Coeffs) + ...
      conv(IM3_BasisFifthOrder,IM3_FifthOrder_Coeffs) +...
      conv(IM3_BasisSeventhOrder,IM3_SeventhOrder_Coeffs);
   %conv(IM3_BasisNinthOrder,IM3_NinthOrder_Coeffs);
   IM3_DPD_Signal = IM3_DPD_Signal(1:end-AdaptiveFilterDelay);
   n = (1:length(IM3_DPD_Signal)).';
   IM3_DPD_Signal = IM3_DPD_Signal.*exp(2*pi*1i*n*IM3_Freq*1e6/SystemFs);
   PAin_IM3_ThirdOrderDPD   = PA_InputSignal + IM3_DPD_Signal;
   
   %Check How good theses still are
   if(USE_WARP)
      PA_Output_IM3_ThirdOrderDPD = WARP_broadcast(PAin_IM3_ThirdOrderDPD);
   else
      PA_Output_IM3_ThirdOrderDPD = MemoryLess_PA(PAin_IM3_ThirdOrderDPD,MemorylessPA_Paramters);
   end
   PA_Output_IM3_ThirdOrderDPD = PA_Output_IM3_ThirdOrderDPD(100:end-100);
   
   
   figure(100);
   plot_freqdomain(PA_OutputSignal,SystemFs,'','r',UpsamplingFactor,POWER_PLOT_1MHZ,Pout_dBm);     %Plot without DPD
   hold on;
   plot_freqdomain(PA_Output_IM3_ThirdOrderDPD,SystemFs,'','r',UpsamplingFactor,POWER_PLOT_1MHZ,Pout_dBm);     %Plot without DPD
   
   %Train with that
   DPD_LearningRate = 3;
   [IM3_ThirdOrder_Coeffs, Coeff_3rd] = BlockDecorrDPD_MEM(PAin_IM3_ThirdOrderDPD,IM3_BasisThirdOrder.',MemoryLessPA,MemoryLessDPD, ...
      SystemFs,IM3_Freq,MemorylessPA_Paramters,MemoryPA_Paramters,Memory_depth,LoopDelay,....
      0,DPD_LearningBlockSize,DPD_FilteringBlockSize,FeedBackFilter,...
      DPD_LearningRate,NumLearningSamples,1,Sparsity_Indx,IM3_ThirdOrder_Coeffs);
   disp('DPD Filtering BlockSize in msec = ');
   disp(1e3*DPD_FilteringBlockSize/SystemFs);
   
   %Check to see if I need more suppression
   
   % Applying Third Order IM3+ Decorr DPD
   AdaptiveFilterDelay  = length(IM3_ThirdOrder_Coeffs) - 1;
   IM3_DPD_Signal = conv(IM3_BasisThirdOrder,IM3_ThirdOrder_Coeffs);
   IM3_DPD_Signal = IM3_DPD_Signal(1:end-AdaptiveFilterDelay);
   n = (1:length(IM3_DPD_Signal)).';
   IM3_DPD_Signal = IM3_DPD_Signal.*exp(2*pi*1i*n*IM3_Freq*1e6/SystemFs);
   PAin_IM3_ThirdOrderDPD   = PAin_IM3_ThirdOrderDPD + IM3_DPD_Signal;
   
   % Check to see that I'm not saturating anything
   MAX_real = max(abs(real(PAin_IM3_ThirdOrderDPD)));
   MAX_imag = max(abs(imag(PAin_IM3_ThirdOrderDPD)));
   MAX = max(MAX_real,MAX_imag);
   display('MAX PAin_IM3_ThirdOrderDPD =');
   display(MAX);
   if(USE_WARP)
      PA_Output_IM3_ThirdOrderDPD = WARP_broadcast(PAin_IM3_ThirdOrderDPD);
   else
      PA_Output_IM3_ThirdOrderDPD = MemoryLess_PA(PAin_IM3_ThirdOrderDPD,MemorylessPA_Paramters);
   end
   PA_Output_IM3_ThirdOrderDPD = PA_Output_IM3_ThirdOrderDPD(100:end-100);
   
   figure(100);
   plot_freqdomain(PA_Output_IM3_ThirdOrderDPD,SystemFs,'','b',UpsamplingFactor,POWER_PLOT_1MHZ,Pout_dBm);
   grid on;
   hold on;
   
   %% Do additional training as needed.
   
   %GenerateBasis(CC1,CC2,1,3,5);
   DPD_LearningRate = 1;
   [IM3_FifthOrder_Coeffs, Coeff_5th] = BlockDecorrDPD_MEM(PAin_IM3_ThirdOrderDPD,IM3_BasisFifthOrder.',MemoryLessPA,MemoryLessDPD, ...
      SystemFs,IM3_Freq,MemorylessPA_Paramters,MemoryPA_Paramters,Memory_depth,LoopDelay,....
      0,DPD_LearningBlockSize,DPD_FilteringBlockSize,FeedBackFilter,...
      DPD_LearningRate,NumLearningSamples,1,Sparsity_Indx,IM3_FifthOrder_Coeffs);
   disp('DPD Filtering BlockSize in msec = ');
   disp(1e3*DPD_FilteringBlockSize/SystemFs);
   
   % Applying 5th order IM3+ Decorr DPD
   IM3_DPD_Signal = conv(IM3_BasisFifthOrder,IM3_FifthOrder_Coeffs);
   IM3_DPD_Signal = IM3_DPD_Signal(1:end-AdaptiveFilterDelay);
   n = (1:length(IM3_DPD_Signal)).';
   IM3_DPD_Signal = IM3_DPD_Signal.*exp(2*pi*1i*n*IM3_Freq*1e6/SystemFs);
   PAin_IM3_FifthOrderDPD   = PAin_IM3_ThirdOrderDPD + IM3_DPD_Signal;
   
   % Check to see that I'm not saturating anything
   MAX_real = max(abs(real(PAin_IM3_FifthOrderDPD)));
   MAX_imag = max(abs(imag(PAin_IM3_FifthOrderDPD)));
   MAX = max(MAX_real,MAX_imag);
   display('MAX PAin_IM3_FifthOrderDPD =');
   display(MAX);
   if(USE_WARP)
      PA_Output_IM3_FifthOrderDPD = WARP_broadcast(PAin_IM3_FifthOrderDPD);
   else
      PA_Output_IM3_FifthOrderDPD = MemoryLess_PA(PAin_IM3_FifthOrderDPD,MemorylessPA_Paramters);
   end
   PA_Output_IM3_FifthOrderDPD = PA_Output_IM3_FifthOrderDPD(100:end-100);
   
   figure(100);
   plot_freqdomain(PA_Output_IM3_FifthOrderDPD,SystemFs,'','k',UpsamplingFactor,POWER_PLOT_1MHZ,Pout_dBm);
   grid on;
   
   %% Do additional training as needed. 7th
   
   %GenerateBasis(CC1,CC2,1,3,7);
   DPD_LearningRate = 1;
   [IM3_SeventhOrder_Coeffs, Coeff_7th] = BlockDecorrDPD_MEM(PAin_IM3_FifthOrderDPD,IM3_BasisSeventhOrder.',MemoryLessPA,MemoryLessDPD, ...
      SystemFs,IM3_Freq,MemorylessPA_Paramters,MemoryPA_Paramters,Memory_depth,LoopDelay,....
      0,DPD_LearningBlockSize,DPD_FilteringBlockSize,FeedBackFilter,...
      DPD_LearningRate,NumLearningSamples,1,Sparsity_Indx,IM3_SeventhOrder_Coeffs);
   disp('DPD Filtering BlockSize in msec = ');
   disp(1e3*DPD_FilteringBlockSize/SystemFs);
   
   % Applying 7th order IM3+ Decorr DPD
   IM3_DPD_Signal = conv(IM3_BasisSeventhOrder,IM3_SeventhOrder_Coeffs);
   IM3_DPD_Signal = IM3_DPD_Signal(1:end-AdaptiveFilterDelay);
   n = (1:length(IM3_DPD_Signal)).';
   IM3_DPD_Signal = IM3_DPD_Signal.*exp(2*pi*1i*n*IM3_Freq*1e6/SystemFs);
   PAin_IM3_SeventhOrderDPD   = PAin_IM3_FifthOrderDPD + IM3_DPD_Signal;
   
   % Check to see that I'm not saturating anything
   MAX_real = max(abs(real(PAin_IM3_SeventhOrderDPD)));
   MAX_imag = max(abs(imag(PAin_IM3_SeventhOrderDPD)));
   MAX = max(MAX_real,MAX_imag);
   display('MAX PAin_IM3_SeventhOrderDPD =');
   display(MAX);
   if(USE_WARP)
      PA_Output_IM3_SeventhOrderDPD = WARP_broadcast(PAin_IM3_SeventhOrderDPD);
   else
      PA_Output_IM3_SeventhOrderDPD = MemoryLess_PA(PA_Output_IM3_SeventhOrderDPD,MemorylessPA_Paramters);
   end
   PA_Output_IM3_SeventhOrderDPD = PA_Output_IM3_SeventhOrderDPD(100:end-100);
   
   figure(100);
   plot_freqdomain(PA_Output_IM3_SeventhOrderDPD,SystemFs,'','g',UpsamplingFactor,POWER_PLOT_1MHZ,Pout_dBm);
   grid on;
   
   %% Do additional training as needed. 9th
   
   %     GenerateBasis(CC1,CC2,1,3,9);
   %     DPD_LearningRate = 1;
   %     [IM3_NinthOrder_Coeffs, Coeff_9th] = BlockDecorrDPD_MEM(PAin_IM3_SeventhOrderDPD,IM3_BasisNinthOrder.',MemoryLessPA,MemoryLessDPD, ...
   %         SystemFs,IM3_Freq,MemorylessPA_Paramters,MemoryPA_Paramters,Memory_depth,LoopDelay,....
   %         0,DPD_LearningBlockSize,DPD_FilteringBlockSize,FeedBackFilter,...
   %         DPD_LearningRate,NumLearningSamples,1,Sparsity_Indx,IM3_NinthOrder_Coeffs);
   %     disp('DPD Filtering BlockSize in msec = ');
   %     disp(1e3*DPD_FilteringBlockSize/SystemFs);
   %
   %     % Applying 9th order IM3+ Decorr DPD
   %     IM3_DPD_Signal = conv(IM3_BasisNinthOrder,IM3_NinthOrder_Coeffs);
   %     IM3_DPD_Signal = IM3_DPD_Signal(1:end-AdaptiveFilterDelay);
   %     n = (1:length(IM3_DPD_Signal)).';
   %     IM3_DPD_Signal = IM3_DPD_Signal.*exp(2*pi*1i*n*IM3_Freq*1e6/SystemFs);
   %     PAin_IM3_NinthOrderDPD   = PAin_IM3_SeventhOrderDPD + IM3_DPD_Signal;
   %
   %     % Check to see that I'm not saturating anything
   %     MAX_real = max(abs(real(PAin_IM3_NinthOrderDPD)));
   %     MAX_imag = max(abs(imag(PAin_IM3_NinthOrderDPD)));
   %     MAX = max(MAX_real,MAX_imag);
   %     display('MAX PAin_IM3_NinthOrderDPD =');
   %     display(MAX);
   %     PA_Output_IM3_NinthOrderDPD = WARP_broadcast(PAin_IM3_NinthOrderDPD);
   %     PA_Output_IM3_NinthOrderDPD = PA_Output_IM3_NinthOrderDPD(100:end-100);
   %
   %     figure(100);
   %     plot_freqdomain(PA_Output_IM3_NinthOrderDPD,SystemFs,'','g',UpsamplingFactor,POWER_PLOT_1MHZ,Pout_dBm);
   %     grid on;
   
end


length_coeffs = length(Coeff_3rd);
Coeff_3rd = [Coeff_3rd; (zeros(2*length_coeffs,1)+IM3_ThirdOrder_Coeffs+IM3_ThirdOrder_Coeffs_start)];
Coeff_5th = [zeros(1*length_coeffs,1)+IM3_FifthOrder_Coeffs_start; Coeff_5th; (zeros(1*length_coeffs,1)+IM3_FifthOrder_Coeffs+IM3_FifthOrder_Coeffs_start)];
Coeff_7th = [zeros(2*length_coeffs,1)+IM3_SeventhOrder_Coeffs_start; Coeff_7th];
%Coeff_9th = [zeros(3*162,1)+IM3_NinthOrder_Coeffs_start; Coeff_9th];

Samples = 1:length(Coeff_3rd);
TimeAxis = (Samples*DPD_FilteringBlockSize/SystemFs)*1e3;

figure();
plot(TimeAxis,abs(Coeff_3rd));
hold on;
plot(TimeAxis,abs(Coeff_5th));
plot(TimeAxis,abs(Coeff_7th));
grid on;
%plot(TimeAxis,abs(Coeff_9th));
