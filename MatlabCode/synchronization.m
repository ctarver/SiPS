function synchronization(Channel_RX, Channel_TX, PA_OutputSignal, PA_InputSignal);
global Spacing SystemFs ts_tx phase delay;

tx_length    = 2^19;

%% Define the preamble
sts_f = zeros(1,64);
sts_f(1:27) = [0 0 0 0 -1-1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0 1+1i 0 0 0 1+1i 0 0 0 1+1i 0 0];
sts_f(39:64) = [0 0 1+1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0 -1-1i 0 0 0 -1-1i 0 0 0 1+1i 0 0 0];
sts_t = ifft(sqrt(13/6).*sts_f, 64);
sts_t = sts_t(1:16);

% LTS for CFO and channel estimation
lts_f = [0 1 -1 -1 1 1 -1 1 -1 1 -1 -1 -1 -1 -1 1 1 -1 -1 1 -1 1 -1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1 1 1 -1 -1 1 1 -1 1 -1 1 1 1 1];
lts_t = ifft(lts_f, 64);

% Use 30 copies of the 16-sample STS for extra AGC settling margin
preamble = [repmat(sts_t, 1, 30)  lts_t(33:64) lts_t lts_t];
preamble = preamble / max(abs(preamble));

%Shift the preamble so it shows up properly at center for RX channel
preamble = preamble.';
preamble = preamble.*exp(2*pi*1i*(0:length(preamble)-1).'*Spacing/SystemFs);


%% GET FREQ OFFSET
payload_length = tx_length - length(preamble);

t = [0:ts_tx:((payload_length - 1) * ts_tx)].'; % Create time vector (Sample Frequency is ts_tx (Hz))
sinusoid_1    = 0.9 * exp(j*2*pi * (20/6) * 1e6 * t);
payload = sinusoid_1;
payload_shifted = payload.*exp(-2*pi*1i*(0:payload_length-1).'*Spacing/SystemFs); %Shift for calculating phase difference later

tx_data  = vertcat(preamble, payload);

rx_iq = WARP_broadcast(tx_data);

% Process Preamble
LTS_CORR_THRESH = 0.8;
payload_ind_array = 0;
FFT_OFFSET = 4;

% Complex cross correlation of Rx waveform with time-domain LTS
lts_corr = abs(conv(conj(fliplr(lts_t)), sign(rx_iq)));

% Skip early and late samples
lts_corr = lts_corr(32:end-32);

% Find all correlation peaks
lts_peaks = find(lts_corr > LTS_CORR_THRESH*max(lts_corr));

% Select best candidate correlation peak as LTS-payload boundary
[LTS1, LTS2] = meshgrid(lts_peaks,lts_peaks);
[lts_second_peak_index,y] = find(LTS2-LTS1 == length(lts_t));

% Stop if no valid correlation peak was found
if(isempty(lts_second_peak_index))
    fprintf('No LTS Correlation Peaks Found!\n');
    return;
end

% Set the sample indices of the payload symbols and preamble
payload_ind = lts_peaks(max(lts_second_peak_index))+32;
payload_ind_array = [payload_ind_array, payload_ind];
lts_ind = payload_ind-160;

if(1) %Apply CFO
    %Extract LTS (not yet CFO corrected)
    rx_lts = rx_iq(lts_ind : lts_ind+159);
    rx_lts1 = rx_lts(-64+-FFT_OFFSET + [97:160]);
    rx_lts2 = rx_lts(-FFT_OFFSET + [97:160]);
    
    %Calculate coarse CFO est
    rx_cfo_est_lts = mean(unwrap(angle(rx_lts1 .* conj(rx_lts2))));
    rx_cfo_est_lts = rx_cfo_est_lts/(2*pi*64);
else
    rx_cfo_est_lts = 0;
end

% Apply CFO correction to raw Rx waveform
rx_cfo_corr_t = exp(1i*2*pi*rx_cfo_est_lts*[0:length(rx_iq)-1]);
rx_dec_cfo_corr = rx_iq .* rx_cfo_corr_t.';

% Re-extract LTS for channel estimate
rx_lts = rx_dec_cfo_corr(lts_ind : lts_ind+159);
rx_lts1 = rx_lts(-64+-FFT_OFFSET + [97:160]);
rx_lts2 = rx_lts(-FFT_OFFSET + [97:160]);

rx_lts1_f = fft(rx_lts1);
rx_lts2_f = fft(rx_lts2);

% Calculate channel estimate
rx_H_est = lts_f .* (rx_lts1_f.' + rx_lts2_f.')/2;

% Extract the payload samples
payload_received = rx_iq(payload_ind : payload_ind+payload_length-1);
delay = payload_ind - length(preamble);

angles = payload_received ./ payload_shifted;
angles2 = payload_shifted ./ payload_received;

old_phase = mean(angles(50:end-50));

% Estimate the frequency offset
Fs = 40e6; % sample rate
angles_dec = resample(angles,1,1000,16); % decimate (and filter out high-frequency content)
angles_dec = angles_dec(11:end-10); % remove transients
L = length(angles_dec);
N = L/2;
cfo = f_cfr_fitz(angles_dec,L,N,Fs/1000) % frequency offset in Hz

Spacing = (Channel_RX - Channel_TX) * (5*10^6) - cfo ; %5 MHz per channel


%% Estimate Phase 
payload_shifted = payload.*exp(-2*pi*1i*(0:payload_length-1).'*Spacing/SystemFs); %New 
angles = payload_received ./ payload_shifted;
angles2 = payload_shifted ./ payload_received;

%throw away first 100 and last 100 and save as phase difference.
phase = mean(angles(100:end-100));
phase = phase / abs(phase);     


 %More Synchronization
 
    PA_OutputSignal = PA_OutputSignal(delay:delay+length(PA_InputSignal)-1);    %Apply delay calculated earlier
    
    payload_shifted = PA_InputSignal.*exp(-2*pi*1i*(0:length(PA_InputSignal)-1).'*Spacing/SystemFs); %Shift for calculating phase difference
    rx_shifted = PA_OutputSignal.*exp(2*pi*1i*(0:length(PA_OutputSignal)-1).'*Spacing/SystemFs); %Shift for calculating phase difference
    
    angles2 = PA_OutputSignal ./ payload_shifted;
    angles3 = rx_shifted ./ PA_InputSignal;
    
    [delay4,phase4,PA_OutputSignal_Synchronized] = cyclosync(payload_shifted, PA_OutputSignal,'Y TO X');
    
    [delay5,phase5,PA_OutputSignal2_Synchronized] = cyclosync(PA_InputSignal, rx_shifted,'Y TO X');
    
    disp('Subsample Delay estimation =')
    disp(delay4)
    
    disp('Cyclosync Phase estimation on LTE signal (shifted input; we are using this one) =')
    disp(angle(phase4))    
    
    disp('Cyclosync Phase estimation on LTE signal (shifted output) =')
    disp(angle(phase5))    
    
    disp('Sinusoid Phase estimation (@+IM3 Freq) =')
    disp(-angle(phase))      
    
    phase = phase4; %Trust cyclosynce
    if(0)
        [delay5,phase5,PA_OutputSignal2_Synchronized] = cyclosync(PA_InputSignal, rx_shifted,'Y TO X');
        
        angle(phase) %From Sinusoid
        angle(phase4) %From cyclosync1
        angle(phase5) %From cyclosync2
        
        FFT_Out = fft(rx_shifted,2^11);
        FFT_Out(1) = 0;
        PAin2 = PA_InputSignal*2;
        FFT_In = fft(PA_InputSignal,2^11);
        FFT_In(1) = 0;
        
        FFT_In(abs(FFT_In) < 5) = 0;
        FFT_Out(abs(FFT_Out) < 5) = 0;
        
        subplot(2,1,1)
        plot(abs(FFT_Out))
        hold on;
        plot(abs(FFT_In))
        
        subplot(2,1,2)
        plot(angle(FFT_Out))
        hold on
        plot(angle(FFT_In))
        
        Diff = angle(FFT_Out)-angle(FFT_In);
        
        figure()
        plot(Diff)
        hold on;
        Line = phase4*ones(size(Diff));
        plot(-angle(Line));
        
        x = input('1 for up, 0 for down')
        if(x == 1)
            for(n = 1:length(Diff))
                if(Diff(n) < 0)
                    Diff(n) = Diff(n)+2*pi;
                end
            end
        else
            for(n = 1:length(Diff))
                if(Diff(n) > 0)
                    Diff(n) = Diff(n)-2*pi;
                end
            end
        end
        
        figure()
        plot(Diff)
        hold on;
        Line = phase4*ones(size(Diff));
        plot(-angle(Line));
        
        
        x = input('phase')
        x = -x;
        [X,Y] = pol2cart(x,1);
        phase = X+Y*i
    end

end