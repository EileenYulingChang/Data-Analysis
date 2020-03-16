clear all; clc; close all;
%[y,Fs] = audioread('music1.wav');
[y,Fs] = audioread('music2.wav');
tr_piano=length(y)/Fs;  % record time in seconds
plot((1:length(y))/Fs,y);
xlabel('Time [sec]'); 
ylabel('Amplitude');
title('Mary had a little lamb (piano)');

v = y'/2;
a = 100;

L = length(v)/Fs;
k = (2*pi/(2*L))*[0:length(v)/2-1 -length(v)/2:-1];
ks = fftshift(k);

tfinal = length(v)/Fs;
t = (1:length(v))/Fs;

Sgtvector = [];
t_step = 0:0.05:10;

figure(2)
for i = 1:length(t_step)
    g = exp(-a*(t-t_step(i)).^2);
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

    subplot(3,1,3), plot(ks, abs(fftshift(Sgt)));
    xlabel('Frequency(hz)');
    ylabel('Amplitude');

    pause(0.0000001)
end

figure(3) 
pcolor(t_step, ks/(2*pi), Sgtvector.'), shading interp
set(gca, 'Ylim', [100 600], 'Fontsize', (14));