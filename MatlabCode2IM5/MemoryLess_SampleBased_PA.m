function PA_OutSample = MemoryLess_SampleBased_PA(PAinSample,Beta_1,Beta_3,Beta_5)

Epsi_1 = PAinSample;
Epsi_3 = PAinSample.*abs(PAinSample).^2;
Epsi_5 = PAinSample.*abs(PAinSample).^4;

PA_OutSample = Beta_1*Epsi_1 + Beta_3*Epsi_3 + Beta_5*Epsi_5;