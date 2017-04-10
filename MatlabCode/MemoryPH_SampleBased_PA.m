function PA_OutSample = MemoryPH_SampleBased_PA(PAinSample,MemoryPA_Paramters,MemoryDepth)                                         

persistent PreviousState1
persistent PreviousState3
persistent PreviousState5
persistent PreviousState7
persistent PreviousState9

PH_f1 = MemoryPA_Paramters(1:MemoryDepth);
PH_f3 = MemoryPA_Paramters(MemoryDepth + 1:2*MemoryDepth);
PH_f5 = MemoryPA_Paramters(2*MemoryDepth + 1:3*MemoryDepth);
PH_f7 = MemoryPA_Paramters(3*MemoryDepth + 1:4*MemoryDepth);
PH_f9 = MemoryPA_Paramters(4*MemoryDepth + 1:end);

Epsi_1 = PAinSample;
Epsi_3 = PAinSample.*abs(PAinSample).^2;
Epsi_5 = PAinSample.*abs(PAinSample).^4;
Epsi_7 = PAinSample.*abs(PAinSample).^6;
Epsi_9 = PAinSample.*abs(PAinSample).^8;

% keyboard;

if isempty(PreviousState1)
    [PH_SampleBranch1 CurrentState1] = filter(PH_f1,1,Epsi_1);
else
    [PH_SampleBranch1 CurrentState1] = filter(PH_f1,1,Epsi_1,PreviousState1);
end
PreviousState1 = CurrentState1;

if isempty(PreviousState3)
    [PH_SampleBranch3 CurrentState3] = filter(PH_f3,1,Epsi_3);
else
    [PH_SampleBranch3 CurrentState3] = filter(PH_f3,1,Epsi_3,PreviousState3);
end
PreviousState3 = CurrentState3;

if isempty(PreviousState5)
    [PH_SampleBranch5 CurrentState5] = filter(PH_f5,1,Epsi_5);
else
    [PH_SampleBranch5 CurrentState5] = filter(PH_f5,1,Epsi_5,PreviousState5);
end
PreviousState5 = CurrentState5;

if isempty(PreviousState7)
    [PH_SampleBranch7 CurrentState7] = filter(PH_f7,1,Epsi_7);
else
    [PH_SampleBranch7 CurrentState7] = filter(PH_f7,1,Epsi_7,PreviousState7);
end
PreviousState7 = CurrentState7;

if isempty(PreviousState9)
    [PH_SampleBranch9 CurrentState9] = filter(PH_f9,1,Epsi_9);
else
    [PH_SampleBranch9 CurrentState9] = filter(PH_f9,1,Epsi_9,PreviousState9);
end
PreviousState9 = CurrentState9;

PA_OutSample = PH_SampleBranch1 + PH_SampleBranch3 + PH_SampleBranch5 +...
               PH_SampleBranch7 + PH_SampleBranch9;

