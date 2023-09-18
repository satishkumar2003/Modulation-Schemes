%% S J Satish Kumar 2021BEC0014

%% This code represents the generation, tranmission and recovery of a Direct Sequence Spread Spectrum signal under BPSK { DS-BPSK }

clear all;
close all;
clc;

%% Code

% Generation of bit pattern

message_1 = [1 0 1];
disp("Message 1 : ");
disp(message_1);

spread_code_1 = [ 1 1 0 0 ];
disp("Spread Code 1: ");
disp(spread_code_1);

% Carrier signal information
f_c = 30;

% Sampling frequency information
f_s = 50*f_c;
t_s = 1/f_s;
t = 0:t_s:(1-t_s);

% Generate bit length from t_s
N = length(spread_code_1);
t_b = length(t)/length(message_1);
t_c = t_b/N;

% Generate the carrier signal

phi_t = sqrt(2/t_c)*cos(2*pi*f_c.*t);

% We use Polar Non Return to Zero as the line coding mechanism
line_coded_message_1 = PNRZ(message_1,t_b);
line_coded_sp_1 = PNRZ(spread_code_1,t_c);

m_t = zeros(1,length(line_coded_message_1));

for i=1:length(m_t)
    relative_index = (mod(i,t_b));
    if mod(relative_index,t_b)==0
        relative_index = mod(i,t_b+1);
    end
    m_t(i) = line_coded_message_1(i)*line_coded_sp_1(relative_index);
end

figure(1)
subplot(2,1,1);
stem(message_1,'filled');
title("Message signal");

subplot(2,1,2);
stem(spread_code_1,'filled');
title("Spread code");


figure(2)

subplot(5,1,1);
plot(t,line_coded_message_1);
title("Line Coded Message signal");

subplot(5,1,2);
plot(line_coded_sp_1);
title("Line Coded Spread code ");

subplot(5,1,3);
plot(t,m_t);
title("Spreaded message signal");


% Modulated signal
s_t = phi_t.*m_t;

% Modulated plots
subplot(5,1,4);
plot(t,phi_t);
title("Carrier Signal");

subplot(5,1,5);
plot(t,s_t);
title("DS-CDMA signal for 1 sender");

% Demodulation ( Coherent detection )
y_t = phi_t.*s_t;

% We then recover the bits by integrating in bit intervals and then using a
% decision device

recovered_spreaded_code = recoverBPSK(y_t,t_c);
disp("Recovered spreaded code: ");
disp(recovered_spreaded_code);

spread_code_polarity = zeros(1,length(spread_code_1));
for i = 1:length(spread_code_polarity)
    if (spread_code_1(i)>0) 
        spread_code_polarity(i) = 1;
    else 
        spread_code_polarity(i) = -1;
    end
end

recovered_decoded_bits = recovered_spreaded_code.*repmat(spread_code_polarity,1,3);

received_bits = [];
for i=1:length(message_1)
    start = (i-1)*N+1;
    stop = i*N;
    received_bits = [received_bits sum(recovered_decoded_bits(start:stop))>0];
end

disp("Received Bits :");
disp(received_bits);

figure(3)
subplot(4,1,1);
plot(t,y_t);
title("Product Modulator output");

subplot(4,1,2);
stem(recovered_spreaded_code,'filled');
title("Recovered bit pattern");

subplot(4,1,3);
stem(recovered_decoded_bits,'filled');
title("Recovered decoded bits");

subplot(4,1,4);
stem(received_bits,'filled');
title("Received/Despreaded bits");

disp("Error in Transmission :");
disp(biterr(message_1,received_bits));

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
function arr = recoverBPSK(signal,t_c)
    no_of_bits = length(signal)/t_c;
    arr = zeros(1,no_of_bits);
    for i = 0:(no_of_bits-1)
        start = i*t_c + 1;
        stop = (i+1)*t_c;
        temp = sum(signal(start:stop)) > 0; % Threshold is 0 as we have used PNRZ
        if(temp)
            arr(i+1) = 1;
        else
            arr(i+1) = -1;
        end
    end
end