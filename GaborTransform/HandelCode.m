clc; close all; clear all;
load handel
v = y'/2;
plot((1:length(v))/Fs,v);
xlabel('Time [sec]');
ylabel('Amplitude');
title('Signal of Interest, v(n)');

a = 500;
ShannonWindow = 0.25;
b = 100;

L = length(v)/Fs;
k = (2*pi/(2*L))*[0:(length(v)-1)/2 -(length(v)-1)/2:-1];
ks = fftshift(k);

tfinal = length(v)/Fs;
t = (1:length(v))/Fs;

Sgtvector = [];
t_step = 0:0.01:5;

figure(2)
for i = 1:length(t_step)
    g = exp(-a*(t-t_step(i)).^2);%Gaussian
    g = (1-b*(t-t_step(i).^2).*exp(-a*(t-t_step(i)).^2); %Mexican Hat
    %Shannon window
    if t_step(i)+ShannonWindow > 8.92
        break
    end
    g = 0*[1:length(v)];
    if t_step(i) <2*ShannonWindow
        tend = floor((t_step(i)+ShannonWindow)*8192;
        g(1:tend) = 1;
    elseif tfinal -t_step(i) < 2*ShannonWindow
        tstart = ceil((t_step(i) *8192))+1;
        g(tstart:end) = 1;
    else
        tstart = floor((t_step(i)-ShannonWindow)*8192);
        tend = floor((t_step(i)+ShannonWindow)*8192);
        g(tstart:tend) = 1;
    end
    Sg = g.*v;
    Sgt = fft(Sg);
    Sgtvector = [Sgtvector; abs(fftshift(Sgt))];
    subplot(3,1,1), plot(t, v, t, g, 'r');
    axis([0 length(v)/Fs -0.5 1]);
    xlabel('Time [sec]'); 
    ylabel('Amplitude');

    subplot(3,1,2), plot(t, Sg);
    axis([0 length(v)/Fs -0.5 1]);
    xlabel('Time [sec]'); 
    ylabel('Amplitude');

    subplot(3,1,3), plot(ks/(2*pi), abs(fftshift(Sgt)));
    xlabel('Frequency(hz)');
    ylabel('Amplitude');
    set(gca, 'Fontsize', 14);

    pause(0.0000001)
end

figure(3) 
pcolor(t_step, ks/(2*pi), Sgtvector.'), shading interp
set(gca, 'Ylim', [0 800], 'Fontsize', (20));
xlabel('Time(sec)');
ylabel('Frequency(hz)');