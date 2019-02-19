clear all;
close all;
clc;


%% without dimer formation
apaf30 = [69701, 16106, 38149, 8531, 999999, 59912, 90010, 4320, 22996, 58126];
% apaf50 = [9717, 10864, 9233, 3485, 5737, 9612, 20965, 11815, 9526, 9037];
% apaf50 = [12298, 19225, 57402, 4743, 66688, 3173, 20142, 5803, 27519, 17792];
apaf50 = [7688, 6738, 5509, 11410, 26994, 5826, 2988, 27930, 20846, 923];
apaf70 = [9226, 18559, 5840, 9923, 3008, 14884, 18185, 9468, 6274, 6616];
apaf90 = [6747, 6668, 4180, 3293, 12079, 9145, 7431, 3879, 6168, 7064];

apaf30 = apaf30.*10^-4;
apaf50 = apaf50.*10^-4;
apaf70 = apaf70.*10^-4;
apaf90 = apaf90.*10^-4;

%% with dimer formation
dimer = [208217, 250456, 41934, 185861, 797145, 417918, 110671, 274836, 198655, 226183];
dimer = dimer.*10^-4;
%% plotting the result
x = [30 50 70 90];
y1 = [mean(apaf30) mean(apaf50) mean(apaf70) mean(apaf90)];    % mean
y2 = [std(apaf30) std(apaf50) std(apaf70) std(apaf90)];    % variance

figure;
subplot(121), plot(x, y1);
xlabel('Initial Apaf count');
ylabel('Average time-to-apoptosome (sec)');
title('Average time-to-apotosome vs apaf levels');

subplot(122), plot(x, y2);
xlabel('Initial Apaf count');
ylabel('Standard deviation of time-to-apoptosome (sec)');
title('Standard Deviation of time-to-apotosome vs apaf levels');
