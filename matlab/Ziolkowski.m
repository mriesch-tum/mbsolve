clear all;
close all;

load ~/mbsolve-build/Ziolkowski-Basic.mat

t = 2;

x = 0:GridPointSize:XDim;

figure;
plot(x, dm11(:, t) - dm22(:, t));


figure;
plot(x, e(:, t));
