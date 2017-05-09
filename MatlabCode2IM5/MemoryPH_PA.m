function PA_Out = MemoryPH_PA(PAin,MemoryPA_Paramters,MemoryDepth)

PH_f1 = MemoryPA_Paramters(1:MemoryDepth);
PH_f3 = MemoryPA_Paramters(MemoryDepth + 1:2*MemoryDepth);
PH_f5 = MemoryPA_Paramters(2*MemoryDepth + 1:3*MemoryDepth);
PH_f7 = MemoryPA_Paramters(3*MemoryDepth + 1:4*MemoryDepth);
PH_f9 = MemoryPA_Paramters(4*MemoryDepth + 1:end);

Epsi_1 = PAin;
PH_Branch_1 = filter(PH_f1,1,Epsi_1);
Epsi_3 = PAin.*abs(PAin).^2;
PH_Branch_3 = filter(PH_f3,1,Epsi_3);
Epsi_5 = PAin.*abs(PAin).^4;
PH_Branch_5 = filter(PH_f5,1,Epsi_5);
Epsi_7 = PAin.*abs(PAin).^6;
PH_Branch_7 = filter(PH_f7,1,Epsi_7);
Epsi_9 = PAin.*abs(PAin).^8;
PH_Branch_9 = filter(PH_f9,1,Epsi_9);

PA_Out = PH_Branch_1 + PH_Branch_3 + PH_Branch_5 + PH_Branch_7 + PH_Branch_9;

% Adding noise to the PA model
if 0
    SNR = 70; % dB
    SNR_ratio = 10^(SNR/10);
    Signal_Pwr = mean(abs(PA_Out).^2);
    Noise_pwr = Signal_Pwr/SNR_ratio;
    Noise = sqrt(Noise_pwr)*randn(size(PA_Out));
    PA_Out = PA_Out + Noise;
end
