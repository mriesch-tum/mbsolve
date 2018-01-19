clear all;
close all;

%  choose hdf5 file
[filename, folder] = uigetfile('../../*.hdf', 'Select result data');
f = fullfile(folder, filename);

% read global attributes
d_x = h5readatt(f, '/', 'gridpoint_size');
d_t = h5readatt(f, '/', 'timestep_size');
t_e = h5readatt(f, '/', 'sim_endtime');
L_x = h5readatt(f, '/', 'dev_length');

% complete grid
x = 0:d_x:L_x;
t = 0:d_t:t_e;



%plot(x, e(:, 76));
%plot(x, inv12(:, 76));

if 0
    %% display "video" of e-field pattern in cavity
    %% takes a while
    e = h5read(f, '/e/real');
    for i = 1:size(e, 2)

    figure(1);
    plot(e(:, i));
    ylabel('E-Field');
    %xlim([0, XDim]);

    pause(0.1)

    end
else
    e = h5read(f, '/e0/real');
    plot(t, e(1, :));
    
    
end