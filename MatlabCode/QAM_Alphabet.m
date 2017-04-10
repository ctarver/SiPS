function alphabet = QAM_Alphabet(ModulationType)

switch ModulationType
    case 1
        MQAM                        =   4;
    case 2
        MQAM                        =   16;
    case 3
        MQAM                        =   64;
end
alphaMqam                           =   -(sqrt(MQAM)-1):2:(sqrt(MQAM)-1);
A                                   =   repmat(alphaMqam,sqrt(MQAM),1);
B                                   =   flipud(A');
const_qam                           =   A+1j*B;
const_qam                           =   const_qam(:);
alphabet                            =   const_qam;

