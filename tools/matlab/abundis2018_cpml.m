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
    plot(x, e(:, i));
    ylabel('E-Field');
    %xlim([0, XDim]);

    pause(0.1)

    end
else
    e = h5read(f, '/e0/real');
    plot(t, e(1, :));
    
    
    I = abs(e).^2;

figure;
plot(t, I);
xlim([1.8e-9, 1.84e-9]);

interval = t > 2e-9 ;

e = e(interval);
t = t(interval);

NFFT = 2^nextpow2(length(e));
win = hanning(length(e)).';

nconstant = normalize_fourier(NFFT,win,@ifft);
Y = ifft(e.*win,NFFT)/nconstant/2; Y = Y(1:NFFT/2);
%Y_p = ifft(P.*win,NFFT)/nconstant/2; Y_p = Y_p(1:NFFT/2);
f = [0:NFFT/2-1]/NFFT*1/d_t;

halfmax = max(abs(Y))*0.3162;

figure;
plot(f, abs(Y));
hold on;
plot(f, ones(size(f)) * halfmax);
xlim([1.5e12, 3.5e12]);
xlabel('Frequency/Hz');
ylabel('Electric field at facet/a.u.');

end
