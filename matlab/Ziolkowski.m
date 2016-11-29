clear all;
close all;

load ~/mbsolve-build/Ziolkowski-Basic.mat

t = 8;

x = 0:GridPointSize:XDim;

figure;
plot(x, dm11(:, t) + dm22(:, t));
xlim([0, XDim]);
ylabel('Trace as sanity check');

figure;
plot(x, dm11(:, t) - dm22(:, t));
ylabel('Population  inversion');
xlim([0, XDim]);

figure;
plot(x, dm11(:, t));
ylabel('dm11');
xlim([0, XDim]);

figure;
plot(x, dm22(:, t));
ylabel('dm22');
xlim([0, XDim]);

figure;
plot(x, e(:, t));
ylabel('E-Field');
xlim([0, XDim]);
