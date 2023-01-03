clc
clear all
close all

addpath('./function/');
%% Parameters
para = para_init();

% %% Generate user location
[user_loc, target_loc, d_BU, d_BT, d_BF] = generate_location(para);
% load("locations.mat");
plot_location(para,user_loc,target_loc);

%% Path loss
path_loss_user = para.pathloss_direct(d_BU)';
path_loss_user = sqrt(10.^((-para.noise - path_loss_user)/10));

path_loss_target = para.pathloss_direct(2*d_BT)';
path_loss_target = sqrt(10.^((-para.noise - path_loss_target)/10));

path_loss_FAP = para.pathloss_direct(d_BF);
path_loss_FAP = sqrt(10.^((-para.noise - path_loss_FAP)/10));

%% Generate channel
[hu, Hr, hd] = generate_channel(para, path_loss_user, path_loss_target, path_loss_FAP);
p = ULA_func(para.target_theta(1)*pi/180, para.N);
p = sqrt(para.Pmax) * p / norm(p);
[~, SINR_sensing, ~] = SINR(para, hu, Hr, hd, p);
SINR_sensing = real(10*log10(SINR_sensing));

%% Partial offloading
[WSR_partial,re,rc,p_partial,WSR_convergence] = algorithm_WMMSE(para, hu, Hr, hd);

figure;
subplot(2,1,1);
plot(WSR_convergence,'-or', 'LineWidth', 1);
xlabel('Number of Iterations');
ylabel("Computation Rate (Mbit/s)");
title("Partial offloading");

subplot(2,1,2);
theta = -90:0.5:90;
[pattern] = beampattern(para,theta, hu, Hr, p_partial);
plot(theta, pattern, '-b', 'LineWidth', 1);
ylim([-60,0]); xlim([-90,90]);
xlabel('Angles (degree)');
ylabel("Beampattern (dB)");



%% Binary offloading

initial_point.re = re;
initial_point.rc = rc;
initial_point.p = p_partial;

[WSR,m,ze,zc,p_binary,WSR_convergence, WSR_convergence_real] = algorithm_ADMM(para, hu, Hr, hd, initial_point);

figure; 
subplot(2,1,1); hold on;
plot(WSR_convergence,'-b','LineWidth', 1);
plot(WSR_convergence_real,'-r','LineWidth', 1);
legend("Objective value", "Effective rate");
ylabel("Computation Rate (Mbit/s)");
xlabel('Number of Iterations');
title("Binary offloading");

subplot(2,1,2);
theta = -90:0.5:90;
[pattern] = beampattern(para,theta, hu, Hr, p_binary);
plot(theta, pattern,  '-b', 'LineWidth', 1);
ylim([-60,0]); xlim([-90,90]);
xlabel('Angles (degree)');
ylabel("Beampattern (dB)");



