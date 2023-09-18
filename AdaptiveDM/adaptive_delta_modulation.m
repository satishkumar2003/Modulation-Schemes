%% 2021BEC0014 S J Satish Kumar

% This code demonstrates the working of Adaptive Delta Modulation to transmit
% signals

clear all;
close all;
clc;

% Parameters
f_message = 3;                 % Frequency of the message signal 
fs = 40*f_message;              % Sampling Frequency
t = 0:1/fs:1;                          % Time scale
A_message = 1;                     % Amplitude of the message signal
message_signal = A_message * sin(2*pi*f_message*t);

% DM Parameters
base_delta = 0.2;
curr_delta = 0;                    % Step size or quantization step

%% Transmitter

output_signal = zeros(size(t));
quantized_signal = zeros(size(t));
prev_value = 0;

for i = 1:length(t)
    if message_signal(i) > prev_value
        quantized_signal(i) = 1;
        curr_delta = curr_delta+base_delta;
        prev_value = prev_value + abs(curr_delta);
    else
        quantized_signal(i) = 0;
        curr_delta = curr_delta-base_delta;
        prev_value = prev_value - abs(curr_delta);
    end
    output_signal(i) = prev_value;
    disp(abs(curr_delta));
end

disp("Quantized Signal : ");
disp(quantized_signal);

% Plots

figure("Name","Transmitter");
subplot(3,1,1);
plot(t, message_signal);
title('Message Signal');
xlabel('Time (s)');
ylabel('Amplitude');
grid on;

subplot(3,1,2);
stem(t, quantized_signal,"filled");
title('Quantized Signal');
xlabel('Time (s)');
ylabel('Amplitude');
ylim([-2 2]);
grid on;

subplot(3,1,3);
stairs(t, output_signal,"Color","red");
hold on
plot(t,message_signal,'Color','blue');
title('Delta Modulated Signal');
xlabel('Time (s)');
ylabel('Amplitude');
legend("Delta Modulated Signal","Message signal");
ylim([-2 2]);
grid on;



%% Receiver 

% Received quantized signal

reconstructed_DM_signal = zeros(1,length(quantized_signal));
prev_value = 0;

curr_delta = 0;

for i=1:length(quantized_signal)
    if quantized_signal(i)==1
        curr_delta = curr_delta+base_delta;
        prev_value = prev_value + abs(curr_delta);
        reconstructed_DM_signal(i) = prev_value;
    else
        curr_delta = curr_delta-base_delta;
        prev_value = prev_value - abs(curr_delta);
        reconstructed_DM_signal(i) = prev_value;
    end
end

f_cutoff = 2*f_message;
filtered_signal = lowpass(reconstructed_DM_signal,f_cutoff,fs);

error_signal = message_signal - filtered_signal;

figure("Name","Receiver")
subplot(3,1,1);
stairs(t,reconstructed_DM_signal);
hold on;
plot(t,filtered_signal);
xlabel("Time (s)");
ylabel("Amplitude");
title("Reconstructed DM signal from Quantized Signal");
legend("Reconstructed DM signal","Low pass filtered signal");
ylim([-2 2]);
grid on;

subplot(3,1,2);
plot(t,message_signal);
hold on;
plot(t,filtered_signal);
xlabel("Time (s)");
ylabel("Amplitude");
title("Message/Received Signal");
legend("Message signal","Received signal");
ylim([-2 2]);
grid on;

subplot(3,1,3);
plot(t,error_signal);
xlabel("Time(s)");
ylabel("Amplitude");
title("Error signal");
ylim([-2 2]);
grid on;