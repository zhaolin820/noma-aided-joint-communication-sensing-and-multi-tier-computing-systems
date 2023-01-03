function [user_loc, target_loc, d_BU, d_BT, d_BF, para] = generate_location(para)
%Generate the user locations randomly 
%  [user_loc, d_BU] = generate_user_location(para)
%Inputs:
%   para: structure of the initial parameters
%Outputs:
%   user_loc: locations of users in meters
%   d_BU: distance between MEN and user
%Date: 30/05/2021
%Edit: 27/09/2021
%Author: Zhaolin Wang

% 
% d_CU = rand(para.K,1) * (para.user_range(2) - para.user_range(1))...
%     + para.user_range(1); % distance from user center to user

user_gap = (para.user_range(2) - para.user_range(1)) / (para.K-1);
d_CU = (0:para.K-1)' * user_gap + para.user_range(1);

user_angle = rand(para.K,1) * (pi) + pi/2; % angle of directions from MEN
user_loc = d_CU.*exp(1i*user_angle);
user_loc = [real(user_loc), imag(user_loc)] + para.user_center;
relative_user_loc = user_loc - para.BS_loc;
d_BU = sqrt(relative_user_loc(:,1).^2 + relative_user_loc(:,2).^2);

d_BT = rand(para.L,1) * (para.target_range(2) - para.target_range(1))...
    + para.target_range(1); % distance from MEN to target
% d_BT = sort(d_BT);

target_loc = d_BT.*exp(1i*(para.target_theta+180)*pi/180);
target_loc = [real(target_loc), imag(target_loc)] + para.BS_loc;

d_BF = para.BS_loc - para.FAN_loc;
d_BF = sqrt(d_BF(1)^2 + d_BF(2)^2);

end

