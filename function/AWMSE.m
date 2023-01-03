function [eu, er, ed] = AWMSE(para, hu, Hr, hd, p, wu_mmse,wr_mmse,wd_mmse,pi_mmse)
%Calculate the augmented weighted MSEs for the WMMSE algorithm
%  [eu, er, ed] = AWMSE(para, hu, Hr, hd, p, wu_mmse,wr_mmse,wd_mmse,pi_mmse)
%Inputs:
%   para: structure of the initial parameters
%   hu: the communication channels from the users to the BS
%   Hr: the round-trip sensing channels
%   hd: the communication channels from the BS to the CS
%   p: the transmit beamformer at the BS
%   wu_mmse: the optimal MMSE receiver for user-BS offloading
%   wr_mmse: the optimal MMSE receiver for sensing
%   wd_mmse: the optimal MMSE receiver for BS-CS offloading
%   pi_mmse: the optimal WMMSE weights
%Outputs:
%   eu: the AWMSE for user-BS offloading
%   er: the AWMSE for sensing
%   ed: the AWMSE for BS-CS offloading
%Date: 28/02/2021
%Author: Zhaolin Wang

%% AWMSE for user-BS offloading
Hr_sum = sum(Hr, 3);
eu = [];
for k = 1:para.K
    hk = hu(:,k);
    if para.mode == 1 % NOMA
       hc_i = hu(:,k+1:para.K);
    elseif para.mode == 0 % SDMA
        hc_i = hu;
        hc_i(:,k) = [];
    end


    w_k = wu_mmse(:,k);
    A = matrix_decompisition(Hr_sum'*(w_k*w_k')*Hr_sum);
    B = para.Pt*(hk*hk') + para.Pt*(hc_i*hc_i') + eye(para.N);
    

    eu_k = p'*(A*A')*p + real(w_k'*B*w_k) - 2*real(sqrt(para.Pt)*w_k'*hk) + 1;
    eu_k = pi_mmse(k)*eu_k - log2(pi_mmse(k));
    eu = [eu, eu_k];
end

%% AWMSE for sensing
G = sum(Hr(:,:,2:end), 3);
H0 = Hr(:,:,1);
A = matrix_decompisition(H0'*(wr_mmse*wr_mmse')*H0 + G'*(wr_mmse*wr_mmse')*G);

if para.mode == 1
    er = p'*(A*A')*p + wr_mmse'*wr_mmse- 2*real(wr_mmse'*H0*p) + 1;
elseif para.mode == 0
    er = p'*(A*A')*p + real(wr_mmse'*(para.Pt * (hu*hu') +eye(para.N))*wr_mmse)- 2*real(wr_mmse'*H0*p) + 1;
end

er = pi_mmse(para.K+1)*er - log2(pi_mmse(para.K+1));

%% AWMSE for BS-CS offloading
ed = abs(wd_mmse)^2*(p'*(hd*hd')*p + 1) - 2*real(conj(wd_mmse)*hd'*p) + 1;
ed = pi_mmse(para.K+2)*ed - log2(pi_mmse(para.K+2));
end

function [U_de] = matrix_decompisition(U)
[P,C] = eig(U) ;
U_de = P * (C.^(1/2)) * P';
end