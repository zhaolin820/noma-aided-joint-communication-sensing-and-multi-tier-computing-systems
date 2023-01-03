function [WSR, re, rc, p, WSR_convergence] = algorithm_WMMSE(para, hu, Hr, hd)
%WMMSE algorithm for the partial offloading
%  [WSR, re, rc, p, WSR_convergence] = algorithm_WMMSE(para, hu, Hr, hd)
%Inputs:
%   para: structure of the initial parameters
%   hu: the communication channels from the users to the BS
%   Hr: the round-trip sensing channels
%   hd: the communication channels from the BS to the CS
%Outputs:
%   WSR: optimal weighted sum rate
%   re: optimal computation rate at the BS
%   rc: optimal computation rate at the CS
%   p: optimal transmit beamformer at the BS
%   WSR_convergence: convergence behaviour of the WSR
%Date: 28/02/2021
%Author: Zhaolin Wang

cvx_solver mosek; 
WSR_diff = 100; WSR = 1;
e = 0.3e-3;

WSR_convergence = [];
%% initialization
p = ULA_func(para.target_theta(1)*pi/180, para.N);
p = sqrt(para.Pmax) * p / norm(p);

step = 1;
%% AO algorithm
while WSR_diff > e
    [pi_mmse] = update_weight(para, hu, Hr, hd, p);
    [wu_mmse,wr_mmse,wd_mmse] = update_receiver(para, hu, Hr, hd, p);
    [re, rc, p] = update_rate_beamformer(para, hu, Hr, hd, wu_mmse,wr_mmse,wd_mmse,pi_mmse);
    WSR_current = sum(para.weight .* (re + rc));
    WSR_diff = abs(WSR_current - WSR)/WSR;
    WSR = WSR_current;
    WSR_convergence = [WSR_convergence, WSR];
    [~, SINR_sensing, ~] = SINR(para,hu,Hr,hd,p);
    gamma = real(SINR_sensing);
    disp(['Iteration ' num2str(step) ', rate - ' num2str(WSR)...
        ', sensing SINR - ' num2str(10*log10(gamma)) ' dB']);
    step = step + 1;
end
end


function [re, rc, p] = update_rate_beamformer(para, hu, Hr, hd, wu_mmse,wr_mmse,wd_mmse,pi_mmse)

c = 1/log(2) + log2(log(2));
cvx_begin quiet
    variable re(para.K,1) 
    variable rc(para.K,1)
    variable p(para.N,1) complex
    [eu, er, ed] = AWMSE(para, hu, Hr, hd, p, wu_mmse,wr_mmse,wd_mmse,pi_mmse);

    % constraints
    log2(1 + para.gamma_min) <= c - er;
    for k = 1:para.K
        re(k) + rc(k) <= para.B * (c-eu(k));
        re(k) >= 0;
        rc(k) >= 0;
    end
    sum(rc) <= para.B * (c-ed);
    para.xi*para.phi^3*sum(pow_p(re, 3)) + p'*p <= para.Pmax;
      
    maximize(sum(para.weight .* (re + rc)));
cvx_end

end




