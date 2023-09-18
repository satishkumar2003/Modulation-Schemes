%% 2021BEC0014 S J SATISH KUMAR
%% This code demonstrates the generation of a Pulse Code Modulated signal from a given analog signal
clear all
close all
clc

%% Generation of the "Analog Signal"
f_signal = 5000;
A_m = 10;
f=1;
f_s = 50;
L=8;
t_b = 30;
t = 0 : 1/f_signal : 1-1/f_signal; % 1000 samples
analog_signal = A_m*sin(2*pi*f.*t);

%% Sampling of the analog_signal
t_sampled = linspace(0,1-(1/f_signal),f_s);
sampled_signal = A_m*sin(2*pi*f.*t_sampled);
no_of_samples = length(sampled_signal);

%% Generating Quantization and Decision Levels
delta = 2*A_m/(L-1); % Step size for quantization levels 
curr = A_m;
quantization_levels = zeros(1,L);
for i = 1:1:L
    quantization_levels(i) = curr;
    curr = curr - delta;
end

disp("Quantization Levels");
disp(quantization_levels);

decision_levels = zeros(1,L-1);

curr = A_m;
for i=1:1:L-1
    decision_levels(i) = curr-(delta/2);
    curr = curr-delta;
end

disp("Decision Levels");
disp(decision_levels);

%% Quantization process

quantized_samples = quantize(sampled_signal,quantization_levels,decision_levels,L);


bits_per_symbol = log2(L);
encoded_bits = [];

for i=1:1:no_of_samples
    q_level = find(quantization_levels==quantized_samples(i)); % Returns index
    pattern = de2bi(L-q_level,bits_per_symbol);
    pattern = flip(pattern);
    encoded_bits = [ encoded_bits pattern ];
end

pulses = UNRZ(encoded_bits,t_b);
sample_nos = linspace(1,no_of_samples,no_of_samples);
t_step = length(t)/no_of_samples;
SAH_quantized_signal = UNRZ(quantized_samples,t_step);
noise_signal = analog_signal-SAH_quantized_signal; 

% Theoretical SQNR and Practical SQNR
SQNR_theoretical = 6*A_m^2/(delta^2);
signal_power_practical = A_m^2/2;
noise_power_practical = var(sampled_signal-quantized_samples);
SQNR_practical = signal_power_practical/noise_power_practical;

disp("L value");
disp(L);
disp("Theoretical SQNR");
disp(SQNR_theoretical);
disp("Practical SQNR");
disp(SQNR_practical);

%% Plots

figure(1)

subplot(2,1,1);
plot(t,analog_signal);
title("Analog Continuous Signal");
xlabel("Time (s)");
ylabel("Amplitude (V)");

subplot(2,1,2);
stem(sample_nos,sampled_signal);
title("Sampled Signal");
xlim([1 no_of_samples]);
xlabel("Sample No.");
ylabel("Amplitude ( V )");




figure(2);

subplot(2,1,1);

hold on
stem(sample_nos,quantized_samples,"r","filled","DisplayName","Quantized Signal");
stem(sample_nos,sampled_signal,"b","DisplayName","Sampled Signal");
yline(quantized_samples,"-g","HandleVisibility","off");
yline(decision_levels,"--m","HandleVisibility","off");
plot([NaN NaN],"-g","DisplayName","Quantization Levels");
plot([NaN NaN],"--m","DisplayName","Decision Levels");
hold off

title("Quantization Process");
legend();
xlim([1 no_of_samples]);
ylim([-12 12]);
xlabel("Sample No.");
ylabel("Amplitude ( V )");


subplot(2,1,2);
stem(encoded_bits);
title("Encoded Bits");
xlim([1 no_of_samples*bits_per_symbol]);
ylim([0 1.1]);
xlabel("Sample No.");
ylabel("Logic");

disp("Encoded Bits");
disp(encoded_bits);




figure(3)

subplot(2,1,1);
stem(encoded_bits);
title("Encoded Bits");
xlim([1 no_of_samples*bits_per_symbol]);
ylim([0 1.1]);
xlabel("Sample No.");
ylabel("Logic");

subplot(2,1,2);
plot(pulses);
title("Pulse Code Modulated Waveform");
xlim([0 t_b*f_s*bits_per_symbol])
ylim([0 1.1]);
xlabel("Sample No.");
ylabel("Logic");

figure(4)
subplot(2,1,1);

hold on
plot(t,analog_signal,"DisplayName","Analog Signal");
plot(t,SAH_quantized_signal,"DisplayName","Quantized Signal");
yline(quantized_samples,"-.m","HandleVisibility","off");
plot([NaN NaN],"-.m","DisplayName","Quantization levels");
hold off

title("Signals");
legend();

subplot(2,1,2);

plot(t,noise_signal);
title("Quantization Noise");


%% Plot between L and SQNR

L_values = [2,4,8,16,32,64];

noise_power = zeros(1,length(L_values));
Th_SQNR = zeros(1,length(L_values));
Pr_SQNR = zeros(1,length(L_values));

for i=1:length(L_values)
    [noise_power(i),Th_SQNR(i),Pr_SQNR(i)] = returnSQNR(A_m,sampled_signal,L_values(i));
end

figure(5)

hold on
plot(L_values,Th_SQNR);
plot(L_values,Pr_SQNR);
plot(L_values,log10(noise_power));
legend("Theoretical SQNR","Practical SQNR","Noise Power");
hold off
title("Relation between L and SQNR");
xlabel("L");
ylabel("SQNR (dB)");

%% Functions

% Quantization function
% if above threshold, then above level
% if equal or lower, then lower level

function q_signal = quantize(signal,q_levels,d_levels,L)
    len = length(signal);
    q_signal = zeros(1,len);
    
    for i=1:1:len
        sample = signal(i);
        qvalue = 0;
        for j=1:1:L-1
            if(sample >= d_levels(j))
                qvalue = q_levels(j);
                break
            end
        end
        if(sample < d_levels(L-1))
            qvalue = q_levels(L);
        end
        q_signal(i) = qvalue;
    end
end

% Function to generate pulses from encoded sequence of bits
function pulses = UNRZ(bit_pattern,tb)
    pulses = [];
    for i = 1:1:length(bit_pattern)
        temp = repmat(bit_pattern(i),1,tb);
        pulses = cat(2,pulses,temp);
    end
end

function [noise_power_practical,SQNR_theoretical,SQNR_practical] = returnSQNR(A_m,sampled_signal,L)
    delta = 2*A_m/(L-1); % Step size for quantization levels 
    curr = A_m;
    quantization_levels = zeros(1,L);
    for i = 1:1:L
        quantization_levels(i) = curr;
        curr = curr - delta;
    end

    decision_levels = zeros(1,L-1);

    curr = A_m;
    for i=1:1:L-1
        decision_levels(i) = curr-(delta/2);
        curr = curr-delta;
    end

    quantized_samples = quantize(sampled_signal,quantization_levels,decision_levels,L);
    
    SQNR_theoretical = log10(6*A_m^2/(delta^2));
    signal_power_practical = A_m^2/2;
    noise_power_practical = var(sampled_signal-quantized_samples);
    SQNR_practical = log10(signal_power_practical/noise_power_practical);
end