function PA_Out = MemoryLess_PA(PAin,MemorylessPA_Paramters)

Epsi_1 = PAin;
Epsi_3 = PAin.*abs(PAin).^2;
Epsi_5 = PAin.*abs(PAin).^4;
Epsi_7 = PAin.*abs(PAin).^6;
Epsi_9 = PAin.*abs(PAin).^8;

PA_Out = MemorylessPA_Paramters(1)*Epsi_1 + ...
         MemorylessPA_Paramters(2)*Epsi_3 + ...
         MemorylessPA_Paramters(3)*Epsi_5 + ...
         MemorylessPA_Paramters(4)*Epsi_7 + ...
         MemorylessPA_Paramters(5)*Epsi_9;

% Adding noise to the PA model
if 0
    SNR = 70; % dB
    SNR_ratio = 10^(SNR/10);
    Signal_Pwr = mean(abs(PA_Out).^2);
    Noise_pwr = Signal_Pwr/SNR_ratio;
    Noise = sqrt(Noise_pwr)*randn(size(PA_Out));
    PA_Out = PA_Out + Noise;
end


