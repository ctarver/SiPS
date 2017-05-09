function LTE_signal = LTE_SignalGenerator(TxSignalType,N_symbols,alphabet,users_bandwidth_subcarriers, ...
                                      FFT_Size,RRC_RollOffFactor,UpsamplingFactor)

% Signal Generation
switch TxSignalType
    
    case 1 % OFDM signal
        
        % symbol_RX                             =   alphabet(randi([1 length(alphabet)], users_bandwidth_subcarriers, N_symbols));
        symbol_RX                               =   alphabet(ceil(length(alphabet)*rand(users_bandwidth_subcarriers,N_symbols)));        
        data_matrix_RX                          =   zeros(FFT_Size,N_symbols);
        data_matrix_RX(2:2+users_bandwidth_subcarriers/2-1,:) = symbol_RX(1:1+users_bandwidth_subcarriers/2-1,:);
        data_matrix_RX(end-users_bandwidth_subcarriers/2+1:end,:) = symbol_RX(users_bandwidth_subcarriers/2+1:end,:);
        OFDM_data                               =   ifft(data_matrix_RX)*sqrt(users_bandwidth_subcarriers);
        signal_RX                               =   OFDM_data(:);   
        pulse_RX                                =   rcosine(1,UpsamplingFactor,'sqrt',RRC_RollOffFactor,50);
        s_over_RX                               =   zeros(UpsamplingFactor*length(signal_RX),1);
        s_over_RX(1:UpsamplingFactor:end)       =   signal_RX;
        S_IQ_RX                                 =   conv(pulse_RX,s_over_RX);
        LTE_signal                              =   S_IQ_RX((length(pulse_RX)+1)/2:end-(length(pulse_RX)+1)/2);
        
    case 2 % SC-FDMA signal
        
        users_subchannel_locations              =   [-(users_bandwidth_subcarriers/2) (users_bandwidth_subcarriers/2)];
        users_subchannel_index                  =   transpose(users_subchannel_locations(1):users_subchannel_locations(2));
        users_N_subcarriers                     =   length(users_subchannel_index(:));
        % create QAM modulated symbols and put them in data_matrix
        % symbol_RX                             =   alphabet(randi([1 length(alphabet)], users_N_subcarriers, N_symbols));
        symbol_RX                               =   alphabet(ceil(length(alphabet)*rand(users_N_subcarriers,N_symbols)));        
        symbol_dft                              =   fft(symbol_RX)/sqrt(length(symbol_RX));
        data_matrix                             =   symbol_dft;
        % Subcarrier Mapping
        index_vector                            =   mod(users_subchannel_index, FFT_Size) + 1;
        index_vector_blocks                     =   repmat(index_vector, 1, N_symbols) + repmat(FFT_Size*[0:N_symbols-1],  users_N_subcarriers, 1);
        symbol_mapped                           =   cell(FFT_Size,N_symbols);
        load_zeros                              =   cellfun('isempty',symbol_mapped);
        symbol_mapped(load_zeros)               =   {0};
        symbol_mapped_double                    =   cell2mat(symbol_mapped);
        symbol_mapped_double(index_vector_blocks)=  data_matrix;
        % M-point IDFT
        s_ifft                                  =   sqrt(FFT_Size).*ifft(symbol_mapped_double);
        % P/S conversion
        signal                                  =   s_ifft(:);
        % Pulse shaping
        pulse                                   =   rcosine(1,UpsamplingFactor,'sqrt',RRC_RollOffFactor,50);
        s_over                                  =   zeros(UpsamplingFactor*length(signal),1);
        s_over(1:UpsamplingFactor:end)          =   signal;
        S_IQ                                    =   conv(pulse,s_over);
        LTE_signal                              =   S_IQ((length(pulse)+1)/2:end-(length(pulse)+1)/2);
             
end

