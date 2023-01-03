function [WSR, m, ze, zc, p, WSR_convergence, WSR_convergence_real] = algorithm_ADMM(para, hu, Hr, hd, initial_point)
%ADMM algorithm for the binary offloading
%  [WSR,m,ze,zc,p,WSR_convergence, WSR_convergence_real] = algorithm_ADMM(para, hu, Hr, hd, initial_point)
%Inputs:
%   para: structure of the initial parameters
%   hu: the communication channels from the users to the BS
%   Hr: the round-trip sensing channels
%   hd: the communication channels from the BS to the CS
%   initial_point: the initial point
%Outputs:
%   WSR: optimal weighted sum rate
%   m: binary decision factor in binary offloading
%   ze: optimal computation rate at the BS
%   zc: optimal computation rate at the CS
%   p: optimal transmit beamformer at the BS
%   WSR_convergence: convergence behaviour of the objective value
%   WSR_convergence_real: convergence behaviour of the effective WSR
%Date: 28/02/2021
%Author: Zhaolin Wang


cvx_solver mosek; 
WSR = 0;
h = 100;
e_inner = 1e-2;
e_outer = 1e-3;
%% initialization

WSR_convergence = [];
WSR_convergence_real = [];

re = initial_point.re;
rc = initial_point.rc;
p = initial_point.p;

ze = re; zc = rc;
m = re ./ (re + rc); m_bar = re ./ (re + rc);
lambda = zeros(para.K,1); lambda_bar = zeros(para.K,1);
mu = zeros(para.K,1); mu_bar = zeros(para.K,1);
rho = 1;

step = 0;
%% ADMM AO algorithm
while (h > e_outer || WSR_diff > e_inner) && step <= 500
    [pi_mmse] = update_weight(para, hu, Hr, hd, p);
    [wu_mmse,wr_mmse,wd_mmse] = update_receiver(para, hu, Hr, hd, p);
    [re, rc, m_bar] = update_r_m_bar(para, ze, zc, m, rho, lambda, lambda_bar, mu, mu_bar);
    [ze, zc, m, p] = update_z_m_p(para, hu, Hr, hd, wu_mmse,wr_mmse,wd_mmse,pi_mmse, re, rc, m_bar, rho, lambda, lambda_bar, mu, mu_bar);
    [lambda, lambda_bar, mu, mu_bar] = update_dual_variable(para, lambda, lambda_bar, mu, mu_bar, rho, ze, zc, re, rc, m, m_bar);  

    WSR_current = sum(para.weight .* (ze + zc));
    WSR_diff = abs(WSR_current - WSR);
    WSR = WSR_current;
    WSR_convergence = [WSR_convergence, WSR]; % objective value
    h = vioalation_func(para, m, m_bar, ze, zc, re, rc);

    m_i = round(m);
    ze_bar = ze; zc_bar = zc;
    ze_bar(m_i==0) = 0;
    zc_bar(m_i==1) = 0;
    WSR_real = sum(para.weight .* (ze_bar + zc_bar)); % effective rate
    WSR_convergence_real = [WSR_convergence_real, WSR_real];
    step = step + 1;

    disp(['Iteration ' num2str(step) ', objective value - ' num2str(WSR)...
    ', effective rate - ' num2str(WSR_real) ', constraint violation - ' num2str(h)]);

end
m_i = round(m);
ze(m_i==0) = 0;
zc(m_i==1) = 0;
WSR = sum(para.weight .* (ze + zc));
end


function [ze, zc, m, p] = update_z_m_p(para, hu, Hr, hd, wu_mmse,wr_mmse,wd_mmse,pi_mmse, re, rc, m_bar, rho, lambda, lambda_bar, mu, mu_bar)
c = 1/log(2) + log2(log(2));
cvx_begin quiet
    variable ze(para.K,1) 
    variable zc(para.K,1)
    variable m(para.K,1)
    variable p(para.N,1) complex
    [eu, er, ed] = AWMSE(para, hu, Hr, hd, p, wu_mmse,wr_mmse,wd_mmse,pi_mmse);

    % constraints
    log2(1 + para.gamma_min) <= c - er;
    for k = 1:para.K

        ze(k) + zc(k) <= para.B * (c-eu(k));
        ze(k) >= 0;
        zc(k) >= 0;
        ze(k) <= m(k)*500;
        zc(k) <= (1-m(k))*500;

        m(k) <= 1;
        m(k) >= 0;
    end
    sum(zc) <= para.B * (c-ed);
    para.xi*para.phi^3*sum(pow_p(ze, 3)) + p'*p <= para.Pmax;

    % objective function
    obj = sum(para.weight .* (ze + zc));
    obj = obj - 1/(2*rho) * sum_square(ze - m.*re + rho*lambda) - 1/(2*rho) * sum_square(zc - (1-m).*rc + rho*lambda_bar);
    obj = obj - 10/(2*rho) * sum_square(m - m_bar + rho*mu) - 10/(2*rho) * sum_square(m.*(1-m_bar) + rho*mu_bar);
    maximize(obj);
cvx_end
end

function [re, rc, m_bar] = update_r_m_bar(para, ze, zc, m, rho, lambda, lambda_bar, mu, mu_bar)
re = zeros(para.K,1); rc = zeros(para.K,1); m_bar = zeros(para.K,1);

for k = 1:para.K
    re(k) = (ze(k) + rho*lambda(k)) / m(k);
    rc(k) = (zc(k) + rho*lambda_bar(k)) / (1-m(k));

    if rc(k) < 0 rc(k) = 0; end
    if re(k) < 0 re(k) = 0; end

    m_bar(k) = (m(k)^2 + (rho/10*mu_bar(k) + 1)*m(k) + rho/10*mu(k)) / (m(k)^2 + 1);
    if m_bar(k) > 1 
        m_bar(k) = 1;
    elseif m_bar(k) < 0
        m_bar(k) = 0;
    end
end

end

function [lambda, lambda_bar, mu, mu_bar] = update_dual_variable(para, lambda, lambda_bar, mu, mu_bar, rho, ze, zc, re, rc, m, m_bar)

for k = 1:para.K
    lambda(k) = lambda(k) + 1/rho * (ze(k) - m(k)*re(k));
    lambda_bar(k) = lambda_bar(k) + 1/rho * (zc(k) - (1-m(k))*rc(k));
    mu(k) = mu(k) + 10/rho * (m(k) - m_bar(k));
    mu_bar(k) = mu_bar(k) + 10/rho * m(k)*(1-m_bar(k));
end

end

function [h, h_sum] = vioalation_func(para, m, m_bar, ze, zc, re, rc)

h_sum = 0;
h = zeros(para.K, 1);
for k = 1:para.K
   h(k) = max( [ abs(m(k)*(1-m_bar(k))), abs(m(k)-m_bar(k)), abs(ze(k) - m(k)*re(k)), abs(zc(k) - (1-m(k))*rc(k)) ] );
end
h = max(h);
end


