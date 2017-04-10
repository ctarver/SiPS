function [LTE_signal CC1 CC2 SystemFs UpsamplingFactor] = LTE_Transmitter(LTE_Bandwidth,CarrierSpacing,NRB1,NRB2,...
                                      N_symbols,TxSignalType,ModulationType,TxScenario)
%I CHANGED THE THING TO DO 2 DIFFERENT SIZES

EXACT_40MHz_Fs = 1;
if EXACT_40MHz_Fs
    %     FFT_Size = 200;
    %     sampling_frequency = 2e6;
    %     UpsamplingFactor = 20;
    %     SystemFs = UpsamplingFactor*sampling_frequency;
    %     NumUsedSubcarriers1 = 100;
    %     NumUsedSubcarriers2 = 100;
    
    %Chance Modification
    SystemFs = 40000000;
    
    
    FFT_Size =   128;
    
    subcarrier_spacing  = 15e3;
    sampling_frequency  = subcarrier_spacing * FFT_Size;
    UpsamplingFactor = round(SystemFs / sampling_frequency);
    RB_length           = 12; % 12 Sub-carriers per RB
    NumUsedSubcarriers1 = RB_length * NRB1;
    NumUsedSubcarriers2 = RB_length * NRB2;
    
    % QAM Modulator Paramters Per Sub-Carrier
    alphabet                                        =   QAM_Alphabet(ModulationType);
    
    % RRC pusle shaping roll off factor
    RRC_RollOffFactor                               =   0.25;  % Roll off factor of Root Raised Cosine Pulse Shape
    
    
    LTE_signal_1     = LTE_SignalGenerator(TxSignalType,N_symbols,alphabet,NumUsedSubcarriers1,...
        FFT_Size,RRC_RollOffFactor,UpsamplingFactor);
    %LTE_signal_1 = LTE_signal_1*1.7;
    
    FFT_Size =   128;
    
    subcarrier_spacing  = 15e3;
    sampling_frequency  = subcarrier_spacing * FFT_Size;
    UpsamplingFactor = round(SystemFs / sampling_frequency);
    RB_length           = 12; % 12 Sub-carriers per RB
    NumUsedSubcarriers1 = RB_length * NRB1;
    NumUsedSubcarriers2 = RB_length * NRB2;
    
    LTE_signal_2     = LTE_SignalGenerator(TxSignalType,N_symbols,alphabet,NumUsedSubcarriers2,...
        FFT_Size,RRC_RollOffFactor,UpsamplingFactor);  %
    
else
    switch LTE_Bandwidth
        case 1.4
            FFT_Size =   128;
        case 3
            FFT_Size =   256;
        case 5
            FFT_Size =   512;
        case 10
            FFT_Size =   1024;
        case 15
            FFT_Size =   1536;
        case 20
            FFT_Size =   2048;
    end
    % Sampling Frequency Calculations
    subcarrier_spacing  = 15e3;
    sampling_frequency  = subcarrier_spacing * FFT_Size;
    RB_length           = 12; % 12 Sub-carriers per RB
    NumUsedSubcarriers1 = RB_length * NRB1;
    NumUsedSubcarriers2 = RB_length * NRB2;
    UpsamplingFactor    = ceil(5*(CarrierSpacing + LTE_Bandwidth)*1e6/sampling_frequency); % For pulse shaping and possible Non-linearity modeling
    SystemFs            = UpsamplingFactor*sampling_frequency;
end


LTE_signal_1 = LTE_signal_1(1:length(LTE_signal_2));


IF_Freq          = (CarrierSpacing/2)*1e6;
n                                   =   1:max(length(LTE_signal_1),length(LTE_signal_2));
Component_Carrier_1                 =   LTE_signal_1.*exp(2*pi*1j*n*IF_Freq/SystemFs).';
Component_Carrier_2                 =   LTE_signal_2.*exp(-2*pi*1j*n*IF_Freq/SystemFs).';

LTE_signal                          =   Component_Carrier_1 + Component_Carrier_2;
% Discard some samples from the end to make the signal cyclic                 
LTE_signal                          =   LTE_signal(1:end - UpsamplingFactor*FFT_Size);              
CC1                                 =   LTE_signal_1(1:end - UpsamplingFactor*FFT_Size);
CC2                                 =   LTE_signal_2(1:end - UpsamplingFactor*FFT_Size);
