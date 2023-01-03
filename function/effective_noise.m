function [Rk] = effective_noise(para, hu, Hr, p, k)
%Calculate the covariance matrix of the effective noise for NOMA and SDMA
%  [Rk] = effective_noise(para, hu, Hr, p, k)
%Inputs:
%   para: structure of the initial parameters
%   hu: the communication channels from the users to the BS
%   Hr: the round-trip sensing channels
%   hd: the communication channels from the BS to the CS
%   k: the index of the current user
%Outputs:
%   Rk: the covariance matrix of the effective noise for decoding the signal from user k
%Date: 28/02/2021
%Author: Zhaolin Wang

Hr_sum = sum(Hr, 3);
if para.mode == 1 % NOMA
    hu_sic = hu(:,k+1:para.K);
    Rk = para.Pt * (hu_sic*hu_sic') + Hr_sum*(p*p')*Hr_sum' + eye(para.N);
elseif para.mode == 0 % SDMA
    hu_i = hu;
    hu_i(:,k) = [];
    Rk = para.Pt * (hu_i*hu_i') + Hr_sum*(p*p')*Hr_sum' + eye(para.N);
elseif para.mode == 2
    Rk = Hr_sum*(p*p')*Hr_sum' + eye(para.N);
else
    disp("wrong mode!");
end

end

