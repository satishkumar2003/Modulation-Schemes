%% S J Satish Kumar 2021BEC0014

%% This code represents the generation, tranmission and recovery of a Quadrature Phase Shift Keyed ( QPSK ) Signal
clear all;
close all;
clc;

%% Code

% Generation of bit pattern

bit_pattern = randi([0,1],1,10);
disp("Bit Pattern : ");
disp(bit_pattern);

% Carrier signal information
f_c = 30;

% Sampling frequency information
f_s = 50*f_c;
t_s = 1/f_s;
t = 0:t_s:(1-t_s);

% Generate bit length from t_s
t_b = length(t)/length(bit_pattern);

% Generate the carrier signal

phi1_t = sqrt(2/t_b)*cos(2*pi*f_c.*t);

phi2_t = sqrt(2/t_b)*sin(2*pi*f_c.*t);


%% Modulation
% modulated signal is of the format s(t) = aI * phi1_t + aQ * phi2_t
% where aI is the inphase coeffecient and aQ is the quadrature coefficient
% Separate into odd { in phase } and even { quadrature } bits
% Perform BPSK on both corresponding parititons and sum them at the end

odd_bits = [];
even_bits = [];

for i =1:length(bit_pattern)
    if (mod(i,2)==0)
        even_bits = [even_bits bit_pattern(i)];
    else
        odd_bits = [odd_bits bit_pattern(i)];
    end
end

% PNRZ is done to both the corresponding bit sequences

EVEN_LC = PNRZ(even_bits,2*t_b); % 2*t_b because the bit sequence length is divided by 2, so to accomodate more samples for dot multiplication 
ODD_LC = PNRZ(odd_bits,2*t_b);

s1_t = ODD_LC.*phi1_t;
s2_t = EVEN_LC.*phi2_t;

% Modulated signal
s_t = s1_t+s2_t;

figure(1)

subplot(5,1,1);
plot(t.*10,PNRZ(bit_pattern,t_b));
xlabel("Bit number");
ylabel("Bit Value");
title("Bit pattern");

subplot(5,1,2);
plot(t.*10,ODD_LC);
xlabel("Bit number");
ylabel("Bit Value");
title("Odd Bits");

subplot(5,1,3);
plot(t.*10,EVEN_LC);
xlabel("Bit number");
ylabel("Bit Value");
title("Even Bits");

subplot(5,1,4);
plot(t,phi1_t);
xlabel("Time ( ms )");
ylabel("Amplitude (V)");
title("Carrier 1 cosine");

subplot(5,1,5);
plot(t,phi2_t);
xlabel("Time ( ms )");
ylabel("Amplitude (V)");
title("Carrier 2 Sine");


figure(2)

subplot(5,1,1);
plot(t.*10,ODD_LC);
xlabel("Bit number");
ylabel("Bit Value");
title("Odd bits");

subplot(5,1,3);
plot(t,s1_t);
xlabel("Time ( ms )");
ylabel("Amplitude (V)");
title("Odd bits modulated \{ Cosine \}");

subplot(5,1,2);
plot(t.*10,EVEN_LC);
xlabel("Bit number");
ylabel("Bit Value");
title("Even bits");

subplot(5,1,4);
plot(t,s2_t);
xlabel("Time ( ms )");
ylabel("Amplitude (V)");
title("Even bits modulated \{ Sine \}");

subplot(5,1,5);
plot(t,s_t);
xlabel("Time ( ms )");
ylabel("Amplitude (V)");
title("Modulated QPSK signal");




%% Demodulation
y1_t = s_t.*phi1_t; % Cosine carrier { top }
y2_t = s_t.*phi2_t; % Sine carrier { bottom }

y1 = integrator(y1_t,2*t_b);
y2 = integrator(y2_t,2*t_b);

received_bit_pattern = zeros(1,length(bit_pattern));
for i = 1:10
    if(mod(i,2))
        received_bit_pattern(i) = y1((i+1)/2);
    else
        received_bit_pattern(i) = y2(i/2);
    end
end


disp("Received Bit Pattern:");
disp(received_bit_pattern);

figure(3);

subplot(4,2,1:2);
plot(t,s_t);
xlabel("Time ( ms )");
ylabel("Amplitude (V)");
title("Modulated QPSK signal");

subplot(4,2,3);
plot(t,y1_t);
xlabel("Time ( ms )");
ylabel("Amplitude (V)");
title("y1(t) \{ Upper product signal \}");

subplot(4,2,4);
plot(t,y2_t);
xlabel("Time ( ms )");
ylabel("Amplitude (V)");
title("y2(t) \{Lower product signal \}");

subplot(4,2,5);
plot(t.*10,PNRZ(y1,2*t_b));
xlabel("Bit number");
ylabel("Bit Value");
title("Recovered odd bits");

subplot(4,2,6);
plot(t.*10,PNRZ(y2,2*t_b));
xlabel("Bit number");
ylabel("Bit Value");
title("Recovered even bits");

subplot(4,2,7:8);
plot(t.*10,PNRZ(received_bit_pattern,t_b));
xlabel("Bit number");
ylabel("Bit Value");
title("Received Bit Pattern");

figure(4)
subplot(2,1,1);
plot(t.*10,PNRZ(bit_pattern,t_b));
xlabel("Bit number");
ylabel("Bit value");
title("Transmitted Bit Pattern");

subplot(2,1,2);
plot(t.*10,PNRZ(received_bit_pattern,t_b));
xlabel("Bit number");
ylabel("Bit value");
title("Received Bit Pattern");

disp("Error in Transmission: ");
disp(biterr(bit_pattern,received_bit_pattern));

%% Utility function for performing Line Coding ( Polar Non Return to Zero )
function line_coded_bits = PNRZ(bit_pattern,t_b)
    line_coded_bits = [];
    for i=1:length(bit_pattern)
        if bit_pattern(i)==1
            x = ones(1,t_b);
        else
            x = (-1).*ones(1,t_b);
        end
        line_coded_bits = cat(2,line_coded_bits,x);
    end
end


%% Utility function for performing integration and recovering bits
function arr = integrator(signal,t_b)
    no_of_bits = length(signal)/t_b;
    arr = zeros(1,no_of_bits);
    for i = 0:(no_of_bits-1)
        start = i*t_b + 1;
        stop = (i+1)*t_b;
        arr(i+1) = sum(signal(start:stop)) > 0; % Threshold is 0 as we have used PNRZ
    end
end