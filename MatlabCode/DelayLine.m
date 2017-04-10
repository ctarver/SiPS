function DelayedSample = DelayLine(InputSample,Delay,Reset)

persistent DelayBuffer;

if Reset
    DelayBuffer = [];
    DelayedSample = [];
    return;
end
if isempty(DelayBuffer)
    DelayBuffer = zeros(1,Delay);
end

% Delay Buffer Implementation
DelayedSample = DelayBuffer(1);
DelayBuffer_Temp = zeros(1,Delay);
DelayBuffer_Temp(1:end-1) = DelayBuffer(2:end);
DelayBuffer_Temp(end) = InputSample;
DelayBuffer = DelayBuffer_Temp;

% for index = 1:Delay-1
%     DelayBuffer(index) = DelayBuffer(index+1);
% end
% DelayBuffer(Delay) = InputSample;


