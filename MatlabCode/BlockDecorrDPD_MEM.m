function [LMSFilterTaps, FinalCoeff] = BlockDecorrDPD_MEM(PA_InputSignal,IM3_Basis_Orthogonal,MemoryLessPA,MemoryLessDPD, ...
     SystemFs,IM3_Freq,MemorylessPA_Paramters,MemoryPA_Paramters,MemoryDepth,LoopDelay,AdditionalDelay,...
     DPD_LearningBlockSize,DPD_FilteringBlockSize,IM3Filter,FilterLearningRate,NumLearningSamples,...
     plot_coeff,Sparsity_Indx,start_coeff)
global USE_WARP delay Spacing;
coeff = [];
delays = [];
NumberOfBasisFunctions = length(IM3_Basis_Orthogonal(:,1));

norm_basis = norm(IM3_Basis_Orthogonal);

% Adaptive Decorrelating DPD using MBF
if MemoryLessDPD 
    NumSamples = NumLearningSamples; % Total number of samples used for learning    
    NumBlocks = floor(NumSamples/DPD_FilteringBlockSize);
    Mu = FilterLearningRate;
    AdaptiveFilterDelay = 1;
    LMSInput = IM3_Basis_Orthogonal;
    W = zeros(NumberOfBasisFunctions,AdaptiveFilterDelay);
    DPD_Coeff = zeros(NumBlocks,NumberOfBasisFunctions);
    Sparsity_Indx = 1;
else
    NumSamples = NumLearningSamples;
    Mu = FilterLearningRate; 
    NumBlocks = floor(NumSamples/DPD_FilteringBlockSize);
    AdaptiveFilterDelay = 2; % 2
    LMSInput = [IM3_Basis_Orthogonal, zeros(NumberOfBasisFunctions,AdaptiveFilterDelay-1)];
    W = zeros(NumberOfBasisFunctions*AdaptiveFilterDelay,1);
    DPD_Coeff = zeros(NumberOfBasisFunctions*AdaptiveFilterDelay,NumBlocks);
    Alpha = zeros(NumberOfBasisFunctions*AdaptiveFilterDelay,1);
    % Sparsity_Indx = 1;
end  

% Reset the delay line
DelayLineReset = 1;
DelayLine(0,0,DelayLineReset);
DelayLineReset = 0;

DPD_BlockIndx = 0;
RMS_data = [];
Correlate = [];
MSE_vector = [];
MU_vector = [];
confidence = 0;

for Sample = 1:DPD_FilteringBlockSize:NumSamples

    PA_In_StartIndx = Sample;
    PA_In_EndIndx   = Sample + DPD_FilteringBlockSize - 1;
    DPD_BlockIndx   = DPD_BlockIndx + 1;
    
    if Sample > (LoopDelay + DPD_FilteringBlockSize) + Sparsity_Indx*(AdaptiveFilterDelay - 1)   
        
        if MemoryLessDPD
            
            % Prepare block of data to be used for learning
            LMS_In_StartIndx = Sample - (LoopDelay + DPD_FilteringBlockSize) + 1;
            LMS_In_EndIndx = LMS_In_StartIndx + DPD_LearningBlockSize - 1;
            LMS_InputBlock = LMSInput(:,LMS_In_StartIndx:LMS_In_EndIndx);
            
            % Correlation between filter input and error block
            MeanCorrelation = (LMS_InputBlock*conj(ErrorBlock))/norm(LMS_InputBlock); %Normalized   
            if (abs(MeanCorrelation) < 0.05)
                confidence = confidence + 1;
            end
            
            if confidence > 7
                Mu=0.7/norm_basis;
            else
                Mu=4/norm_basis;
            end
            MU_vector = [MU_vector ; Mu];
            Correlate = [Correlate; MeanCorrelation];             
            % Block LMS filter update
            W = W - (Mu).*MeanCorrelation;
            
            % Store DPD filter coeff.
            DPD_Coeff(DPD_BlockIndx,:) = W';
                        
            % Apply LMS block filtering
            Alpha = W';
            LMS_FilteringBlock = LMSInput(:,PA_In_StartIndx:PA_In_EndIndx);
            IM3_Block = zeros(size(LMS_FilteringBlock(1,:))).';
            IM3_Block = IM3_Block + (Alpha*LMS_FilteringBlock).';
            
            % Add the LMS filter output to the PA input
            PA_InBlock = PA_InputSignal(PA_In_StartIndx:PA_In_EndIndx) ...
                + IM3_Block.*exp(2*pi*1i*(PA_In_StartIndx:PA_In_EndIndx).'*IM3_Freq*1e6/SystemFs);
           
        else
            
            % Prepare block of data to be used for learning
            LMS_In_StartIndx = Sample - (LoopDelay + DPD_FilteringBlockSize) + 1;
            LMS_In_EndIndx = LMS_In_StartIndx + DPD_LearningBlockSize - 1;
                           
            LMS_InputBlock = [LMSInput(:,LMS_In_StartIndx:LMS_In_EndIndx); ...
                              LMSInput(:,LMS_In_StartIndx - 1*Sparsity_Indx:LMS_In_EndIndx - 1*Sparsity_Indx)];
            
            % Correlation between filter input and error block
            MeanCorrelation = LMS_InputBlock*conj(ErrorBlock)/norm(LMS_InputBlock);;
            
            LMS_FilteringBlock = [LMSInput(:,PA_In_StartIndx:PA_In_EndIndx); ...
                                  LMSInput(:,PA_In_StartIndx - 1*Sparsity_Indx:PA_In_EndIndx - 1*Sparsity_Indx)];

            W = W - (Mu).*MeanCorrelation;
            Alpha = W';
            DPD_Coeff(:,DPD_BlockIndx) = W';
            IM3_Block = Alpha*LMS_FilteringBlock;
            
            % Add the LMS filter output to the PA input
            PA_InBlock = PA_InputSignal(PA_In_StartIndx:PA_In_EndIndx) ...
                + (IM3_Block.*exp(2*pi*1i*(PA_In_StartIndx:PA_In_EndIndx)*IM3_Freq*1e6/SystemFs)).';
            
        end
    else
        PA_InBlock = PA_InputSignal(PA_In_StartIndx:PA_In_EndIndx);             
    end
    
    % PA moldel with PA input block with DPD included
    if(USE_WARP)
        PA_OutBlock = WARP_broadcast(PA_InBlock);
        PA_OutBlock = PA_OutBlock(delay:delay+length(PA_InBlock)-1);
        
        PA_OutBlock_shifted = PA_OutBlock.*exp(2*pi*1i*(PA_In_StartIndx:PA_In_EndIndx).'*Spacing/SystemFs); %Shift for calculating phase difference
        
        [delay4,phase4,PA_OutBlock_cyclosync] = cyclosync(PA_InBlock, PA_OutBlock_shifted,'Y TO X');
        PA_OutBlock_sync = PA_OutBlock * phase4; %PA_angle;  
        
        %coeff = [coeff phase4];
        %delays = [delays delay4];
        
        %Apply Synchronization
        PA_OutBlock = PA_OutBlock_sync;
    else
        if MemoryLessPA
            PA_OutBlock = MemoryLess_PA(PA_InBlock,MemorylessPA_Paramters);
        else
            PA_OutBlock = MemoryPH_SampleBased_PA(PA_InBlock,MemoryPA_Paramters,MemoryDepth);
        end
    end
    
    % Shift the PA output such that the IM3 frequency is at baseband
    PA_OutBlockShifted = PA_OutBlock.*exp(-2*pi*1i*(PA_In_StartIndx+(delay4/2):PA_In_EndIndx+(delay4/2)).'*((IM3_Freq*1e6)-Spacing)/SystemFs);
    
    % IM3 Selection Filter to pick up the IM3 signal used for DPD learning
    if Sample == 1
        [IM3FilteredBlock, Current_State] = filter(IM3Filter,1,PA_OutBlockShifted);
    else
        [IM3FilteredBlock, Current_State] = ... 
            filter(IM3Filter,1,PA_OutBlockShifted,Previous_State);
    end
    Previous_State = Current_State;
    
    % Implement any possible additional delay due to hardware latency
    if AdditionalDelay ~= 0
        IM3FilteredBlockDelayed = zeros(size(IM3FilteredBlock));
        for BlockSample = 1:length(IM3FilteredBlock)
            IM3FilteredSample = IM3FilteredBlock(BlockSample);
            IM3FilteredSampleDelayed = DelayLine(IM3FilteredSample,AdditionalDelay,DelayLineReset);
            IM3FilteredBlockDelayed(BlockSample) = IM3FilteredSampleDelayed;
        end
    else
        IM3FilteredBlockDelayed = IM3FilteredBlock;
    end

    % Extract the IM3 block from the PA output in the feedback receiver
    ErrorBlock = IM3FilteredBlockDelayed(1:DPD_LearningBlockSize);
    MSE_vector = [MSE_vector; ErrorBlock'*ErrorBlock];
   
    RMS_data = [RMS_data; rms(ErrorBlock)];
    
    
    
%     if(rms(ErrorBlock) >= 0.005)
%         Mu = 10;
%     else
%         Mu=1;
%     end
    
end
b = [0.9 0.1];
RMS_data_filtered = filter(b,1,RMS_data);
FinalCoeff = DPD_Coeff(1:NumBlocks,:)+start_coeff;

figure();
plot(abs(FinalCoeff)/max(abs(FinalCoeff)))
hold on
plot(MU_vector/max(MU_vector))
plot(abs(Correlate))

% Plot DPD filter taps convergence
if MemoryLessDPD
    if plot_coeff
        FinalCoeff = DPD_Coeff(1:NumBlocks,:)+start_coeff;
        Samples = 1:length(FinalCoeff(:,1));
        TimeAxis = (Samples*DPD_FilteringBlockSize/SystemFs)*1e3;
        figure;plot(TimeAxis,real(FinalCoeff(:,1)),'k');hold on;plot(TimeAxis,imag(FinalCoeff(:,1)),'--k')
        if NumberOfBasisFunctions >= 2
            plot(TimeAxis,real(FinalCoeff(:,2)),'r');hold on;;plot(TimeAxis,imag(FinalCoeff(:,2)),'--r')
        end
        if NumberOfBasisFunctions >= 3
            plot(TimeAxis,real(FinalCoeff(:,3)),'b');hold on;plot(TimeAxis,imag(FinalCoeff(:,3)),'--b')
        end
        if NumberOfBasisFunctions >= 4
            plot(TimeAxis,real(FinalCoeff(:,4)),'g');hold on;plot(TimeAxis,imag(FinalCoeff(:,4)),'--g')
            %L2= legend('$abs(\alpha_3(n))$','$abs(\alpha_5(n))$','$abs(\alpha_7(n))$','$abs(\alpha_9(n))$');
            %set(L2,'Interpreter','Latex');
        end
        xlabel('Time in msecs'); ylabel('DPD Filter Taps')
        grid on; 
        %setStyle(gcf, 'ieee_column', false)
    end
    LMSFilterTaps = Alpha;
        
else
    if plot_coeff
        Blocks = 1:length(DPD_Coeff(1,:));
        TimeAxis = (Blocks*DPD_FilteringBlockSize/SystemFs)*1e3;
        figure;plot(TimeAxis,abs(DPD_Coeff(1,:)),'r');hold on;
        plot(TimeAxis,abs(DPD_Coeff(2,:)),'g');hold on;
        if NumberOfBasisFunctions >= 2
            plot(TimeAxis,abs(DPD_Coeff(3,:)),'b');hold on;
            plot(TimeAxis,abs(DPD_Coeff(4,:)),'k');hold on;
        end
        if NumberOfBasisFunctions >= 3
            plot(TimeAxis,abs(DPD_Coeff(5,:)),'r--');hold on;
            plot(TimeAxis,abs(DPD_Coeff(6,:)),'g--');hold on;
        end
        if NumberOfBasisFunctions >= 4
            plot(TimeAxis,abs(DPD_Coeff(7,:)),'b--');hold on;
            plot(TimeAxis,abs(DPD_Coeff(8,:)),'k--');hold on;
            L2= legend('$abs(\alpha_{3,0}(n))$','$abs(\alpha_{5,0}(n))$','$abs(\alpha_{7,0}(n))$',...
                           '$abs(\alpha_{9,0}(n))$','$abs(\alpha_{3,1}(n))$','$abs(\alpha_{5,1}(n))$',...
                           '$abs(\alpha_{7,1}(n))$','$abs(\alpha_{9,1}(n))$');
            set(L2,'Interpreter','Latex');
        end
        xlabel('Time in msecs'); ylabel('DPD Filter Coefficients')
        grid on; 
        setStyle(gcf, 'ieee_column', false)
    end
    
    if 0
        Blocks = 1:length(DPD_Coeff(:,1,1));
        TimeAxis = (Blocks*DPD_FilteringBlockSize/SystemFs)*1e3;
        FinalCoeff = DPD_Coeff(1:NumBlocks,:);
        
        figure;plot(TimeAxis,abs(FinalCoeff(:,1)),'k');hold on;
        plot(TimeAxis,abs(FinalCoeff(:,1)),'k');
        if NumberOfBasisFunctions >= 2
            hold on;
            plot(TimeAxis,abs(FinalCoeff(:,2,1)),'r');hold on;
            plot(TimeAxis,abs(FinalCoeff(:,2,2)),'r');hold on;
        end
        if NumberOfBasisFunctions >= 3
            hold on;
            plot(TimeAxis,abs(FinalCoeff(:,3,1)),'b');hold on;
            plot(TimeAxis,abs(FinalCoeff(:,3,2)),'b');hold on;
        end
        if NumberOfBasisFunctions >= 4
            hold on;
            plot(TimeAxis,abs(FinalCoeff(:,4,1)),'g');hold on;
            plot(TimeAxis,abs(FinalCoeff(:,4,2)),'g');
            L2= legend('$abs(\alpha_{3,0}(n))$','$abs(\alpha_{3,1}(n))$','$abs(\alpha_{3,2}(n))$',...
                       '$abs(\alpha_{5,0}(n))$','$abs(\alpha_{5,1}(n))$','$abs(\alpha_{5,2}(n))$',...
                       '$abs(\alpha_{7,0}(n))$','$abs(\alpha_{7,1}(n))$','$abs(\alpha_{7,2}(n))$',...
                       '$abs(\alpha_{9,0}(n))$','$abs(\alpha_{9,1}(n))$','$abs(\alpha_{9,2}(n))$');
            set(L2,'Interpreter','Latex');
        end
        xlabel('Time in msecs'); ylabel('Sub-band DPD Coefficients')
        grid on; 
        setStyle(gcf, 'ieee_column', false)
    end
    
    LMSFilterTaps = DPD_Coeff(:,end);
end


