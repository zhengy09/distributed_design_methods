
%% Figure -- convergence results

clc;clear;close all
load ConData.mat

figure;
histogram(ConIterR,10);
xlabel('Iteration'); ylabel('Frequency')