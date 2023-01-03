function [pi_mmse] = update_weight(para, hu, Hr, hd, p)
%Update the weights in the WMMSE algorithm
%  [pi_mmse] = update_weight(para, hu, Hr, hd, p)
%Inputs:
%   para: structure of the initial parameters
%   hu: the communication channels from the users to the BS
%   Hr: the round-trip sensing channels
%   hd: the communication channels from the BS to the CS
%   p: the transmit beamformer at the BS
%Outputs:
%   pi_mmse: the optimal WMMSE weights
%Date: 28/02/2021
%Author: Zhaolin Wang


%% WMMSE weight for user-EAP offloading
pi_mmse = zeros(para.K+2,1);
for k = 1:para.K
    hk = hu(:,k);
    Rk = effective_noise(para,hu,Hr,p,k);
    pi_mmse(k) = real((1 + para.Pt * hk'*inv(Rk)*hk)/log(2));
end

%% WMMSE weight for sensing
G = sum(Hr(:,:,2:end), 3);
H0 = Hr(:,:,1);
if para.mode == 1
    B = G*(p*p')*G' + eye(para.N);
elseif para.mode == 0
    B = G*(p*p')*G' + para.Pt * (hu*hu') +eye(para.N);
end
pi_mmse(para.K+1) = real((1 + p'*H0'*inv(B)*H0 * p)/log(2));

%% WMMSE weight for EAP-FAP offloading
pi_mmse(para.K+2) = real((1 + hd'*(p*p')*hd)/log(2));

end



