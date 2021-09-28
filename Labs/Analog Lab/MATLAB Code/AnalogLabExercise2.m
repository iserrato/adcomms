clear
fs = 300e3; % this is the sample rate
fc = 107.9e6; % this is the center frequency

x = zeros(3e6,1); % empty vector to store data

% create object for RTL-SDR receiver
rx = comm.SDRRTLReceiver('CenterFrequency',fc, 'EnableTunerAGC', false, 'TunerGain', 35,  'SampleRate', fs);

counter = 1; % initialize a counter
while(counter < length(x)) % while the buffer for data is not full
    rxdata = rx();   % read from the RTL-SDR
    x(counter:counter + length(rxdata)-1) = rxdata; % save the samples returned
    counter = counter + length(rxdata); % increment counter
end
% the data are returned as complex numbers
% separate real and imaginary part, and remove any DC offset
y_I = real(x)-mean(real(x));
y_Q = imag(x)-mean(imag(x));

%% Exercise 2
clear
load('ex2_jeep.mat')

%%

t = [1:length(y_I)]*(1/fs);

figure(1)
plot(t, y_I);
title('Time domain signal of y_I')
xlabel('Time (seconds)')
ylabel('Magnitude')

figure(2)
plot(t, y_Q);
title('Time domain signal of y_Q')
xlabel('Time (seconds)')
ylabel('Magnitude')

t_short = t(:, 1:end-1);
t_sshort = t(:, 2:end-1);

diff_yi = diff(y_I);
diff_yq = diff(y_Q);

figure(3)
plot(t_short, diff_yi);
title('Derivative of y_I')
xlabel('Time (seconds)')
ylabel('Magnitude')

figure(4)
plot(t_short, diff_yq);
title('Derivative of y_Q')
xlabel('Time (seconds)')
ylabel('Magnitude')


diff1 = diff_yq .* y_I(1:end-1,:);
diff2 = diff_yi .* y_Q(1:end-1,:);
diffs = diff1 - diff2;

diffs_nm = diffs - mean(diffs);
diffs_norm = diffs * 0.1 / max(abs(diffs));
diffs_dec = decimate(diffs_norm, 4);

figure(5)

plot(t_short, diffs_norm);
title('Demodulated Signal')
xlabel('Time (seconds)')
ylabel('Magnitude')

%(122000:122300)

% figure(8)
% plot(t_short, diffs_norm);

% sound(diffs_dec, 300000/4)
