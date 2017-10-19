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
rho11 = real(d11);
rho22 = real(d22);
rho33 = real(d33);

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
hold on;
plot(t, rho22);
plot(t, rho33);
if (use_reference == 1)
    plot(comp_data(:, 1) * 1e-15, comp_data(:, 2));
end
ylabel('populations');
xlim([0, SimEndTime]);

trace = rho11 + rho22 + rho33;
trace_err = max(trace) - min(trace);
display(trace_err)

rho11_min = min(rho11)
rho22_min = min(rho22)
rho33_min = min(rho33)
rho11_max = max(rho11)
rho22_max = max(rho22)
rho33_max = max(rho33)
