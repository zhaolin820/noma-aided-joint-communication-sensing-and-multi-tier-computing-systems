function [SINR_user_BS, SINR_sensing, SINR_BS_CS] = SINR(para, hu, Hr, hd, p)
%Calculate the SINR for communication and sensing
%  [SINR_user_BS, SINR_sensing, SINR_BS_CS] = SINR(para, hu, Hr, hd, p)
%Inputs:
%   para: structure of the initial parameters
%   hu: the communication channels from the users to the BS
%   Hr: the round-trip sensing channels
%   hd: the communication channels from the BS to the CS
%   p: the transmit beamformer at the BS
%Outputs:
%   SINR_user_BS: SINR for decoding the offloading signals at the BS
%   SINR_sensing: SINR for sensing
%   SINR_BS_CS: SINR for decoding the offloading signals at the CS
%Date: 28/02/2021
%Author: Zhaolin Wang

G = sum(Hr(:,:,2:end), 3);
H0 = Hr(:,:,1);

% NOMA or SDMA?
if para.mode == 1
    B = G*(p*p')*G' + eye(para.N);
elseif para.mode == 0
    B = G*(p*p')*G' + para.Pt * (hu*hu') +eye(para.N);
end

SINR_sensing = p'*H0'*inv(B)*H0*p;
SINR_user_BS = zeros(para.K,1);
for k = 1:para.K
    hk = hu(:,k);
    Rk = effective_noise(para,hu,Hr,p,k);
    SINR_user_BS(k) = para.Pt*hk'*inv(Rk)*hk;
end

SINR_BS_CS = abs(hd'*p)^2;

end

