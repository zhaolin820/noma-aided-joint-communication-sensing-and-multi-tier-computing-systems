function [wu_mmse,wr_mmse,wd_mmse] = update_receiver(para, hu, Hr, hd, p)
%Update the receivers in the WMMSE algorithm
%  [wu_mmse,wr_mmse,wd_mmse] = update_receiver(para, hu, Hr, hd, p)
%Inputs:
%   para: structure of the initial parameters
%   hu: the communication channels from the users to the BS
%   Hr: the round-trip sensing channels
%   hd: the communication channels from the BS to the CS
%   p: the transmit beamformer at the BS
%Outputs:
%   wu_mmse: the optimal MMSE receiver for user-BS offloading
%   wr_mmse: the optimal MMSE receiver for sensing
%   wd_mmse: the optimal MMSE receiver for BS-CS offloading
%Date: 28/02/2021
%Author: Zhaolin Wang



%% MMSE receiver for user-BS offloading
wu_mmse = zeros(para.N, para.K);
for k = 1:para.K
    hk = hu(:,k);
    Rk = effective_noise(para,hu,Hr,p,k);
    wu_mmse(:,k) = inv(para.Pt*(hk*hk') + Rk)*sqrt(para.Pt)*hk;
end

%% MMSE receiver for sensing
G = sum(Hr(:,:,2:end), 3);
H0 = Hr(:,:,1);
if para.mode == 1
    B = G*(p*p')*G' + eye(para.N);
elseif para.mode == 0
    B = G*(p*p')*G' + para.Pt * (hu*hu') +eye(para.N);
end
wr_mmse = inv(H0*(p*p')*H0' + B)*H0*p;

%% MMSE receiver for BS-CS offloading
wd_mmse = hd'*p / ( hd'*(p*p')*hd + 1 );

end

