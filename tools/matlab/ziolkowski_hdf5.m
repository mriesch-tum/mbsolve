clear all;
close all;

% read out hdf5 file

% TODO: unify with original ziolkowski in order to handle both file formats
% with one file

d_x = h5readatt('~/CPH/Work/mbsolve/build/Ziolkowski_Basic.hdf', '/', 'gridpoint_size');
d_t = h5readatt('~/CPH/Work/mbsolve/build/Ziolkowski_Basic.hdf', '/', 'timestep_size');

x = 0:d_x:150e-6;

e = h5read('~/CPH/Work/mbsolve/build/Ziolkowski_Basic.hdf', '/e');
inv12 = h5read('~/CPH/Work/mbsolve/build/Ziolkowski_Basic.hdf', '/inv12');

%plot(x, e(:, 76));
plot(x, inv12(:, 76));
