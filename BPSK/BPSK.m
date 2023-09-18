%% S J Satish Kumar 2021BEC0014

%% This code represents the generation, tranmission and recovery of a Binary Phase Shift Key ( BPSK ) signal

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
f_s = 32*f_c;
t_s = 1/f_s;
t = 0:t_s:(1-t_s);

% Generate bit length from t_s
t_b = length(t)/length(bit_pattern);

% Generate the carrier signal

phi_t = sqrt(2/t_b)*cos(2*pi*f_c.*t);

% We use Polar Non Return to Zero as the line coding mechanism
LC_bits = PNRZ(bit_pattern,t_b);

% Modulated signal
s_t = phi_t.*LC_bits;

% Modulated plots
figure(1);

subplot(3,1,2);
plot(t,phi_t);
title("Carrier signal");
xlabel("Time ( ms )");
ylabel("Amplitude ( V )");

subplot(3,1,1);
plot(t,LC_bits);
title("Line coded bit pattern");
xlabel("Time ( ms )");
ylabel("Amplitude ( V )");

subplot(3,1,3);
plot(t,s_t);
title("Modulated BPSK Signal");
xlabel("Time ( ms )");
ylabel("Amplitude ( V )");


% Demodulation ( Coherent detection )
y_t = phi_t.*s_t;

% We then recover the bits by integrating in bit intervals and then using a
% decision device

recovered_bit_array = recoverBPSK(y_t,t_b);
disp("Recovered bit pattern : ");
disp(recovered_bit_array);

% Demodulated plots
figure(2);

subplot(3,1,1);
plot(t,s_t);
title("Modulated signal");
xlabel("Time ( ms )");
ylabel("Amplitude ( V )");

subplot(3,1,2);
plot(t,y_t);
title("Signal after product modulator");
xlabel("Time ( ms )");
ylabel("Amplitude ( V )");

subplot(3,1,3);
plot(t,PNRZ(recovered_bit_array,t_b)); % Only for visual purposes we are passing it to the line coded function again
title("Recovered bits");
xlabel("Time ( ms )");
ylabel("Amplitude ( V )");


% Error checking
result = biterr(bit_pattern,recovered_bit_array);
disp("Presence of error : ");
disp(result);

% Final bit patterns plots
figure(3);

subplot(2,1,1);
plot(t,LC_bits);
title("Transmitted bits");
xlabel("Time ( ms )");
ylabel("Amplitude ( V )");

subplot(2,1,2);
plot(t,PNRZ(recovered_bit_array,t_b));
title("Recovered bits");
xlabel("Time ( ms )");
ylabel("Amplitude ( V )");


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