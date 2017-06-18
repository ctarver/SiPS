function f = plot_freqdomain(Rx_Signal,fs,TITLE,colour,UpsamplingFactor,POWER_PLOT_1MHZ,PA_Power_Measured)
hold off;
if POWER_PLOT_1MHZ
    
    L = length(Rx_Signal);
    NFFT = 2^nextpow2(L); % Next power of 2 from length of Rx_Signal
   
    Rx_Signal_FD = fftshift(fft(Rx_Signal,NFFT)/L);
    Rx_Power_FD = (abs(Rx_Signal_FD).^2).';
    Bins_1MHz = round(UpsamplingFactor*1000/15);    
    Rx_1MHz_Power = Bins_1MHz*smooth(Rx_Power_FD,Bins_1MHz,'moving');
    
    f = (fs/2)*linspace(-1,1,NFFT);
    % Plot double-sided amplitude spectrum.
    plot(f/1e6,10*log10(Rx_1MHz_Power)+PA_Power_Measured,colour) 
    title(TITLE)
    xlabel('Frequency (MHz)')
    ylabel('Power in 1 MHz (dBm)')
    axis([-fs/(2*1e6) fs/(2*1e6) -Inf Inf])
    %ylim(axes1,[-80 0]);
else
    Rx_Signal = detrend(Rx_Signal);
    Nfft    = 2048;
    Window  = kaiser(2000,9); 
    Signal_PSD = 10*log10(fftshift(pwelch(Rx_Signal,Window)));
    PA_Max = max(Signal_PSD);
    var = (-1:2/Nfft:1-2/Nfft)*(fs/(2e6))+2437;
    plot([(-1:2/Nfft:1-2/Nfft)*(fs/(2e6))+2437],[-35*ones(length(var),1)],'Color','k')
    hold on;
    plot((-1:2/Nfft:1-2/Nfft)*(fs/(2e6))+2437,Signal_PSD-PA_Max,'LineWidth',1,'Color','k'); 
    title(TITLE)
    xlabel('RF Frequency (MHz)')
    ylabel('Normalized Power (dB)')
    %axis([-fs/(2*1e6)+1 fs/(2*1e6)-1 -Inf Inf])
end