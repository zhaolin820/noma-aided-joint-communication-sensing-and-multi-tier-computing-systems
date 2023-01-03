function [pattern] = beampattern(para,theta, hu, Hr, p)
%The beampattern of ULA in a DFRC BS with seperate depolyment
%  [pattern] = beampattern(P,Rx,Mc,Mr,theta_degree)
%Inputs:
%   theta: range of direction angle
%   hu: the communication channels from the users to the BS
%   Hr: the round-trip sensing channels
%   p: the transmit beamformer at the BS
%Outputs:
%   pattern: beampattern
%Date: 28/02/2021
%Author: Zhaolin Wang

theta = theta*pi/180;
pattern = zeros(length(theta),1);
G = sum(Hr(:,:,2:end), 3);
H0 = Hr(:,:,1);
if para.mode == 1
    B = G*(p*p')*G' + eye(para.N);
elseif para.mode == 0
    B = G*(p*p')*G' + para.Pt * (hu*hu') +eye(para.N);
end

w = inv(B)*H0*p;


for i = 1:length(theta)
    t = theta(i);
    a = ULA_func(t,para.N);
    A = a*a.';
    if para.mode == 1
        pattern(i) = norm(w'*A*p)^2;
    elseif para.mode == 0
        pattern(i) = norm(w'*(A*p + sum(hu,2)))^2;
    end
end

pattern = 10*log10(pattern ./ max(pattern));
end

