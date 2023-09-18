%% S J Satish Kumar 2021BEC0014

%% This code represents the generation, transmission and recovery of a Binary Frequency Shift Key ( BFSK )

clear all;
close all;
clc;


%% Code

bit_pattern = randi([0,1],1,10);
disp("Bit Pattern :- ");
disp(bit_pattern);

inv_bit_pattern = ~bit_pattern;

% We have to decide the two frequencies for BFSK as f1 and f2
f1 = 100;
f2 = 20;

% Sampling terms
fs = 30*max(f1,f2);
ts = 1/fs;

% Samples for time
t = 0:ts:1-ts;

% To get bit length
tb = length(t)/length(bit_pattern);

% Two carrier signals

phi_1 = sqrt(2/tb)*cos(2*pi*f1.*t);
phi_2 = sqrt(2/tb)*cos(2*pi*f2.*t);

LC_bits = UNRZ(bit_pattern,tb);
inv_LC_bits = UNRZ(inv_bit_pattern,tb);

% Modulated signal
s_t = phi_1.*LC_bits + phi_2.*inv_LC_bits;


% Demodulation

y1_t = s_t.*phi_1;
y2_t = s_t.*phi_2;

threshold = max(phi_2)/2;

x1 = recover(y1_t,tb,threshold);
x2 = recover(y2_t,tb,threshold);

bit_pattern_received = x1>x2;
disp("Bit Pattern Received :- ");
disp(bit_pattern_received);

result = biterr(bit_pattern,bit_pattern_received);
disp("Presence of bit error?");
disp(boolean(result));

temp = UNRZ(bit_pattern_received,tb);

%% Plots

% Inputs
figure(1);


subplot(3,1,1);
plot(t,phi_1);
xlabel("Time ( ms )");
ylabel("Amplitude ( V ) ");
title("High Frequency Carrier");

subplot(3,1,2);
plot(t,phi_2);
xlabel("Time ( ms )");
ylabel("Amplitude ( V ) ");
title("Low Frequency Carrier");

subplot(3,1,3);
plot(t,LC_bits);
xlabel("Time ( ms )");
ylabel("Amplitude ( V ) ");
title("UNRZ Bit Pattern");


% Modulated signal
figure(2);

subplot(3,1,1);
plot(t,LC_bits);
xlabel("Time ( ms )");
ylabel("Amplitude ( V ) ");
title("UNRZ Bit Pattern");

subplot(3,1,2);
plot(t,inv_LC_bits);
xlabel("Time ( ms )");
ylabel("Amplitude ( V ) ");
title("UNRZ Inverted Bit Pattern");

subplot(3,1,3);
plot(t,s_t);
xlabel("Time ( ms )");
ylabel("Amplitude ( V ) ");
title("Modulated Signal");

% Demodulation
figure(3);

subplot(4,1,1);
plot(t,s_t);
xlabel("Time ( ms )");
ylabel("Amplitude ( V ) ");
title("Modulated Signal");

subplot(4,1,2);
plot(t,y1_t);
xlabel("Time ( ms )");
ylabel("Amplitude ( V ) ");
title("Upper Hand of Demodulator y1(t)");

subplot(4,1,3);
plot(t,y2_t);
xlabel("Time ( ms )");
ylabel("Amplitude ( V ) ");
title("Lower Hand of Demodulator y2(t)");

subplot(4,1,4);
plot(t,temp);
xlabel("Time ( ms )");
ylabel("Amplitude ( V ) ");
title("Received Bit Pattern");

%% Function generated for Line Coding with UNRZ method

function line_coded_bits = UNRZ(bit_pattern,tb)
    % We have to perform Line Coding ( UniPolar Non Return to Zero ) 
    line_coded_bits = [];
    for i = 1:1:length(bit_pattern)
        if(bit_pattern(i)==1)
            temp = ones(1,tb);
        else
            temp = zeros(1,tb);
        end
        line_coded_bits = cat(2,line_coded_bits,temp);
    end
end

%% Function generated to act as integrator and return recovered binary bit pattern

function bit_arr = recover(x_t,tb,threshold)
    no_of_bits = length(x_t)/tb;
    bit_arr = zeros(1,no_of_bits);
    for i = 1:no_of_bits
        bit_arr(i) = sum(x_t( (i-1)*tb+1 : i*tb )) > threshold;
    end
end