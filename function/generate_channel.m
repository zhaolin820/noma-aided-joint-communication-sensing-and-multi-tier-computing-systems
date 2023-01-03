function [hu, Hr, hd] = generate_channel(para, path_loss_user, path_loss_target, path_loss_FAP)
%Generate the BS-user channels 
%  [h] = generate_channel(para, angle, path_loss)
%Inputs:
%   para: structure of the initial parameters
%   user_angle: struture of the angles
%   path_loss: structure of the path loss
%Outputs:
%   d: BS-user channels
%Date: 14/07/2021
%Edit: 27/09/2021
%Author: Zhaolin Wang

%% BS-user channel
hu = 1/sqrt(2) .* ( randn(para.N,para.K) + 1i*randn(para.N,para.K) );
hu = path_loss_user .* hu;

%% BS-target channel
Hr = zeros(para.N, para.N, para.L);
for l = 1:para.L
    hr = ULA_func(para.target_theta(l)*pi/180, para.N);
    Hr(:,:,l) = path_loss_target(l) * (hr*hr.');
end

%% BS-CS channel
hd = path_loss_FAP * 1/sqrt(2) .* ( randn(para.N,1) + 1i*randn(para.N,1) );

end

