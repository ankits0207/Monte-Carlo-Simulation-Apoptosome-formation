clc;
clear all;
close all;

fileID = fopen('apoptosomes','r');
formatSpec = '%f';
y = fscanf(fileID,formatSpec);

timestep = 10^-4;
x = 0:timestep:timestep*length(y);

subplot(211),plot(x(1:end-1),y);
title('Total Apoptosome Count');
xlabel('Simulation Time (sec)');
ylabel('Apoptosome Count');

fileID = fopen('dimers','r');
formatSpec = '%f';
y = fscanf(fileID,formatSpec);

subplot(212),plot(x(1:end-1),y);
title('Total Dimer Count');
xlabel('Simulation Time (sec)');
ylabel('Dimer Count');
