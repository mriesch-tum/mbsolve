clear all;
close all;

[filename, folder] = uigetfile('../../*.mat', 'Select result data');
load(fullfile(folder, filename));

[filename, folder] = uigetfile('../../*.csv', 'Select compare trace');
use_reference = 0;
if (filename ~= 0)
    comp_data = csvread(fullfile(folder, filename));
    use_reference = 1;
end

% explicitely non-complex quantities to real numbers
e = real(e);
rho11 = real(inv12);

% spatial plot
%x = 0:GridPointSize:XDim;
%figure;
%plot(x, e(:, 76));
%ylabel('E-Field');
%xlim([0, XDim]);



%figure;
%plot(x, inv12(:, t));
%%xlim([0, XDim]);
%if (use_reference == 1)
%    hold on;
%    plot(comp_data(:, 1) * 1e-6, comp_data(:, 2));
%end
%ylabel('Population  inversion');

% temporal plot
t = 0:TimeStepSize:SimEndTime;
figure;
plot(t, e);
ylabel('E-Field');
xlim([0, SimEndTime]);

figure;
plot(t, rho11);
if (use_reference == 1)
    hold on;
    plot(comp_data(:, 1) * 1e-15, comp_data(:, 2));
end
ylabel('\rho_{11}');
xlim([0, SimEndTime]);
