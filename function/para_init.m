function [para] = para_init()
%Construct a struct of the initial values for all the parameters 
%  [values] = para_init()
%Inputs:
%   None
%Outputs:
%   values: a struct
%Date: 08/02/2022
%Author: Zhaolin Wang

para.noise = -80; % noise powe in dBm

para.N = 8; % overall antennas

para.Pmax = 10^(30/10); % overall transmit power at BS
para.Pt = 10^(15/10);
para.n = 1; % equivalent noise power
para.K = 3; % user number
para.target_theta = [0, -60, -30, 30, 60]'; % the first one is the desired directions, the others are the directions of interference sources
para.L = length(para.target_theta);
para.weight = ones(para.K,1); % weight of weighted sum rate 
% para.weight = [0.8,1,1.2]';
para.pathloss_direct =  @(d) 30 + 30*log10(d); % path loss with d in m
para.gamma_min = 10^(30/10);


para.phi = 300;
para.xi = 10^-11; % energy needs for each cycle of CPU (mJ/cycle)
para.B = 30; % bandwidth (MHz)

para.BS_loc = [0,0]; % location of the BS
para.FAN_loc = [100,0]; %location of the CS
para.target_range = [50,50]; % range of target location from BS in m
para.user_range = [60, 60]; % range of user location from user center in m
para.user_center = [0,0];

para.mode = 1; % 1 for NOMA; 0 for SDMA;

end




