function GenerateBasis(CC1,CC2,right,spur,order)
global IM3_BasisThirdOrder IM3_BasisFifthOrder IM3_BasisSeventhOrder IM3_BasisNinthOrder


if ~right       %Swap the contents of CC1 and CC2
    x = CC1;    %Temporary variable 
    CC1 = CC2;
    CC2 = x;
end 

if(spur == 3)
    if(order == 3)
       IM3_BasisThirdOrder = (conj(CC2).*(CC1.^2)); %All basis depend on this.
       IM3_BasisThirdOrder = IM3_BasisThirdOrder ./ norm(IM3_BasisThirdOrder);
    end
    if(order == 5)
       IM3_BasisFifthOrder = (conj(CC2).*(CC1.^2)).*(2*(abs(CC1).^2) + 3*(abs(CC2).^2));  
       projection_5_3 = dot(IM3_BasisThirdOrder,IM3_BasisFifthOrder) / (norm(IM3_BasisThirdOrder)^2);
       IM3_BasisFifthOrder = IM3_BasisFifthOrder - projection_5_3 * IM3_BasisThirdOrder;  %Orthogonalize 
       IM3_BasisFifthOrder = IM3_BasisFifthOrder ./ norm(IM3_BasisFifthOrder);                                        %Normalize
    end   
    if(order == 7)
      IM3_BasisSeventhOrder = (conj(CC2).*(CC1.^2)).*(3*(abs(CC1).^4) + 6*(abs(CC2).^4) + ...
                12*(abs(CC1).^2).*(abs(CC2).^2));
      projection_7_5 = dot(IM3_BasisFifthOrder,IM3_BasisSeventhOrder) / (norm(IM3_BasisFifthOrder)^2);
      projection_7_3 = dot(IM3_BasisThirdOrder,IM3_BasisSeventhOrder) / (norm(IM3_BasisThirdOrder)^2);
      IM3_BasisSeventhOrder = IM3_BasisSeventhOrder - projection_7_5 * IM3_BasisFifthOrder;
      IM3_BasisSeventhOrder = IM3_BasisSeventhOrder - projection_7_3 * IM3_BasisThirdOrder;
      IM3_BasisSeventhOrder = IM3_BasisSeventhOrder ./ norm(IM3_BasisSeventhOrder);            
    end    
    if(order == 9)
      IM3_BasisNinthOrder   = (conj(CC2).*(CC1.^2)).*(4*(abs(CC1).^6) + 10*(abs(CC2).^6) + ...
          30*(abs(CC1).^4).*(abs(CC2).^2) + ...
          40*(abs(CC1).^2).*(abs(CC2).^4));
      IM3_BasisNinthOrder = IM3_BasisNinthOrder - dot(IM3_BasisNinthOrder,IM3_BasisSeventhOrder)*IM3_BasisSeventhOrder;
      IM3_BasisNinthOrder = IM3_BasisNinthOrder - dot(IM3_BasisNinthOrder,IM3_BasisFifthOrder)*IM3_BasisFifthOrder;
      IM3_BasisNinthOrder = IM3_BasisNinthOrder - dot(IM3_BasisNinthOrder,IM3_BasisThirdOrder)*IM3_BasisThirdOrder;
      IM3_BasisNinthOrder = IM3_BasisNinthOrder./norm(IM3_BasisNinthOrder);      
    end 
end
%%QR OPTION 
% IM3_Basis_7th = [IM3_BasisThirdOrder IM3_BasisFifthOrder IM3_BasisSeventhOrder];
% [Q,~] = qr(IM3_Basis_7th,0);
% IM3_Basis_Orth_7th = Q.';
% 
% IM3_BasisThirdOrder   = IM3_Basis_Orth_7th(1,:).';
% IM3_BasisFifthOrder   = IM3_Basis_Orth_7th(2,:).';
% IM3_BasisSeventhOrder = IM3_Basis_Orth_7th(3,:).';

if(spur == 5)
    %% Need to Update this later
    IM5_Basis_5thOrder  = (conj(CC2).^2).*(CC1.^3);
    IM5_Basis_7thOrder  = IM5_Basis_5thOrder.*(4*(abs(CC2).^2) + 3*(abs(CC1).^2));
    IM5_Basis_9thOrder  = IM5_Basis_5thOrder.*(10*(abs(CC2).^4) + 6*(abs(CC1).^4) + ...
        20*(abs(CC2).^2).*(abs(CC1).^2));
    IM5_Basis_9th = [IM5_Basis_5thOrder IM5_Basis_7thOrder IM5_Basis_9thOrder];
    [Q,~] = qr(IM5_Basis_9th,0);
    IM5_Basis = Q.';
    
    IM3_BasisThirdOrder = IM5_Basis_5thOrder; %Hack for loopdelay
end    