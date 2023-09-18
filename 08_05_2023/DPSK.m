%% S J Satish Kumar 2021BEC0014

%% This code represents the generation, tranmission and recovery of a Differential Phase Shift Key ( DPSK ) signal

clear all;
close all;
clc;

%% Code

% Generation of bit pattern

bit_pattern = randi([0,1],1,9);
disp("Bit Pattern : ");
disp(bit_pattern);

% Carrier signal information
f_c = 30;

% Sampling frequency information
f_s = 32*f_c;
t_s = 1/f_s;
t = 0:t_s:(1-t_s);

% Generate bit length from t_s
t_b = length(t)/(length(bit_pattern)+1);

% Generate the carrier signal

phi_t = sqrt(2/t_b)*cos(2*pi*f_c.*t);

% Modulator

reference_bit = 1;

encoded_bit_pattern = zeros(1,length(bit_pattern)+1);
encoded_bit_pattern(1) = reference_bit;

for i=2:length(bit_pattern)+1 % First bit is reference
    encoded_bit_pattern(i) = ~xor(encoded_bit_pattern(i-1),bit_pattern(i-1));
end

disp("Encoded Bit Pattern: ");
disp(encoded_bit_pattern);

ebp_t = PNRZ(encoded_bit_pattern,t_b);
s_t = phi_t.*ebp_t;


% Demodulator
y_t = s_t.*phi_t;
recovered_bits = recoverBPSK(y_t,t_b);

disp("Recovered Bits: ");
disp(recovered_bits);

decoded_bits = zeros(1,length(bit_pattern));

%   Recovered_bits(1) is the reference bit
for i=1:length (bit_pattern)
    decoded_bits(i) = ~xor(recovered_bits(i),recovered_bits(i+1));
end

disp("Decoded Bits: ");
disp(decoded_bits);

disp("Reference bit used");
disp(reference_bit);

disp("Operation used for encoding :- ");
disp("XNOR");

fprintf("Error in transmission :- %d\n",biterr(bit_pattern,decoded_bits));

figure(1)

subplot(4,1,1);
stem(bit_pattern);
xlabel("Bit Number");
ylabel("Logic");
title("Bit Pattern");

subplot(4,1,2);
plot(t,ebp_t);
xlabel("Time ( ms )");
ylabel("Amplitude ( V ) ");
title("Line-coded encoded bits { PNRZ }");


subplot(4,1,3);
plot(t,phi_t);
xlabel("Time ( ms )");
ylabel("Amplitude ( V ) ");
title("Carrier Signal");


subplot(4,1,4);
plot(t,s_t);
xlabel("Time ( ms )");
ylabel("Amplitude ( V ) ");
title("Modulated Signal");


figure(2)

subplot(3,1,1);
plot(y_t);
xlabel("Time ( ms )");
ylabel("Amplitude ( V ) ");
title("Product Modulator output");


subplot(3,1,2);
stem(recovered_bits);
xlabel("Bit number");
ylabel("Logic");
title("Recovered bits");


subplot(3,1,3);
stem(decoded_bits);
xlabel("Bit number");
ylabel("Logic");
title("Decoded Bits");


figure(3)
subplot(2,1,1);
stem(bit_pattern);
xlabel("Bit Number");
ylabel("Logic");
title("Transmitted Bits");


subplot(2,1,2);
stem(decoded_bits);
xlabel("Bit Number");
ylabel("Logic");
title("Received Bits");


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
function arr = recoverBPSK(signal,t_b)
    no_of_bits = length(signal)/t_b;
    arr = zeros(1,no_of_bits);
    for i = 0:(no_of_bits-1)
        start = i*t_b + 1;
        stop = (i+1)*t_b;
        arr(i+1) = sum(signal(start:stop)) > 0; % Threshold is 0 as we have used PNRZ
    end
end