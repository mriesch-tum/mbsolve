% marskar2011  Displays result of marskar2011 setup.
%              Reference data available in tools/reference-data.

% mbsolve: An open-source solver tool for the Maxwell-Bloch equations.
%
% Copyright (c) 2016, Computational Photonics Group, Technical University
% of Munich.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software Foundation,
% Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA

clear all;
close all;

%  choose hdf5 file
[filename, folder] = uigetfile('../../*.hdf', 'Select result data');
f = fullfile(folder, filename);
if (filename == 0)
   return;
end

% try to determine level count from filename
parts = split(filename, { '-', '_' });
if length(parts) > 1
    N = str2num(erase(parts{2}, 'lvl'));
    if isempty(N)
        N = 6;
    end
else
    N = 6;
end
display([ 'Assuming ', num2str(N), ' levels' ]);

comp_data = {};
[filenames, folder] = uigetfile('../reference-data/*.csv', ...
    'Select compare trace', 'MultiSelect', 'on');
if iscell(filenames)
    for i = 1:length(filenames)
        comp_data{i} = csvread(fullfile(folder, filenames{i}));
    end
elseif (filenames == 0)
    % noop
else
    comp_data{1} = csvread(fullfile(folder, filenames));
end

% read global attributes
d_x = h5readatt(f, '/', 'gridpoint_size');
d_t = h5readatt(f, '/', 'timestep_size');
t_e = h5readatt(f, '/', 'sim_endtime');
L_x = h5readatt(f, '/', 'dev_length');

% time grid
t = 0:d_t:t_e;

% data
e = h5read(f, '/e/real');
for i = 1:N
   rho(i, :) = h5read(f, [ '/d', num2str(i), num2str(i), '/real' ]);
end

figure;
plot(t, e);
xlabel('Time/s');
ylabel('E-Field/Vm^{-1}');
xlim([0, t_e]);

figure;
plot(t, rho(1 ,:), 'DisplayName', '\rho_{11}');
hold on;
for i = 2:N
  plot(t, rho(i, :), 'DisplayName', ...
      [ '\rho_{', num2str(i), num2str(i), '}' ]);
end

for i = 1:length(comp_data)
    plot(comp_data{i}(:, 1) * 1e-12, comp_data{i}(:, 2), ...
        'DisplayName', 'Ref');
end
xlabel('Time/s');
ylabel('Populations');
xlim([0, t_e]);
ylim([0, 1]);
legend show;
