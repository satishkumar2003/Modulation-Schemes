%% S J Satish Kumar 2021BEC0014

%% This code represents the generation, tranmission and recovery of a DS-BPSK signal from two different users using Code Division Multiple Access by using orthogonal spread codes

clear all;
close all;
clc;

%% Code

% Generation of bit pattern

message_1 = [1 0 1];
disp("Message 1 : ");
disp(message_1);

message_2 = [0 0 1];
disp("Message 2: ");
disp(message_2);

spread_code_1 = [ 1 1 0 0 ];
disp("Spread Code 1: ");
disp(spread_code_1);

spread_code_2 = [1 0 1 0];
disp("Spread Code 2: ");
disp(spread_code_2);

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

line_coded_message_2 = PNRZ(message_2,t_b);
line_coded_sp_2 = PNRZ(spread_code_2,t_c);

m1_t = zeros(1,length(line_coded_message_1));
m2_t = zeros(1,length(line_coded_message_2));

for i=1:length(m1_t)
    relative_index = (mod(i,t_b));
    if mod(relative_index,t_b)==0
        relative_index = mod(i,t_b+1);
    end
    m1_t(i) = line_coded_message_1(i)*line_coded_sp_1(relative_index);
    m2_t(i) = line_coded_message_2(i)*line_coded_sp_2(relative_index);
end

figure(1)
subplot(4,1,1);
stem(message_1,'filled');
title("Message signal 1");

subplot(4,1,2);
stem(spread_code_1,'filled');
title("Spread code 1");

subplot(4,1,3);
stem(message_2,'filled');
title("Message signal 2");

subplot(4,1,4);
stem(spread_code_2,'filled');
title("Spread code 2");


figure(2)

subplot(4,1,1);
plot(t,line_coded_message_1);
title("Line Coded Message signal 1");

subplot(4,1,2);
plot(t,line_coded_message_2);
title("Line Coded Message signal 2");

subplot(4,1,3);
plot(line_coded_sp_1);
title("Line Coded Spread code 1");

subplot(4,1,4);
plot(line_coded_sp_2);
title("Line Coded Spread code 2");


% Modulated signal
s1_t = phi_t.*m1_t;
s2_t = phi_t.*m2_t;
s_t = s1_t + s2_t;

% Modulated plots
figure(3)

subplot(4,2,1:2);
plot(t,phi_t);
title("Carrier Signal");

subplot(4,2,3);
plot(t,m1_t);
title("Message signal 1");

subplot(4,2,4);
plot(t,m2_t);
title("Message signal 2");

subplot(4,2,5);
plot(t,s1_t);
title("Transmission signal user 1");

subplot(4,2,6);
plot(t,s2_t);
title("Transmission signal user 2")

subplot(4,2,7:8);
plot(t,s_t);
title("Transmitted signal");

% Demodulation ( Coherent detection )
y_t = phi_t.*s_t; % Same for both users as they have the same local carrier

% We then recover the bits by integrating in bit intervals and then using a
% decision device

recovered_spreaded_code = recoverBPSK(y_t,t_c); % Both users will recover the same spreaded spread code as all blocks are same till the receiving part
disp("Recovered spreaded code: ");
disp(recovered_spreaded_code);

spread_code_polarity_1 = zeros(1,length(spread_code_1));
for i = 1:length(spread_code_polarity_1)
    if (spread_code_1(i)>0) 
        spread_code_polarity_1(i) = 1;
    else 
        spread_code_polarity_1(i) = -1;
    end
end

spread_code_polarity_2 = zeros(1,length(spread_code_2));
for i = 1:length(spread_code_polarity_2)
    if (spread_code_2(i)>0)
        spread_code_polarity_2(i) = 1;
    else
        spread_code_polarity_2(i) = -1;
    end
end


recovered_decoded_bits_user_1 = recovered_spreaded_code.*repmat(spread_code_polarity_1,1,3);
recovered_decoded_bits_user_2 = recovered_spreaded_code.*repmat(spread_code_polarity_2,1,3);

disp("Recovered decoded bits user 1");
disp(recovered_decoded_bits_user_1);
disp("Recovered decoded bits user 2");
disp(recovered_decoded_bits_user_2);

received_bits_user_1 = [];
for i=1:length(message_1)
    start = (i-1)*N+1;
    stop = i*N;
    received_bits_user_1 = [received_bits_user_1 sum(recovered_decoded_bits_user_1(start:stop))>0];
end

received_bits_user_2 = [];
for i = 1:length(message_2)
    start = (i-1)*N+1;
    stop = i*N;
    received_bits_user_2 = [received_bits_user_2 sum(recovered_decoded_bits_user_2(start:stop))>0];
end

disp("Received Bits for user 1:");
disp(received_bits_user_1);

disp("Received Bits for user 2:");
disp(received_bits_user_2);

figure(4)
subplot(6,1,1);
plot(t,y_t);
title("Product Modulator output");

subplot(6,1,2);
stem(recovered_spreaded_code,'filled');
title("Recovered bit pattern");

subplot(6,1,3);
stem(recovered_decoded_bits_user_1,'filled');
title("Recovered decoded bits 1");

subplot(6,1,4);
stem(recovered_decoded_bits_user_2,'filled');
title("Recovered decoded bits user 2");

subplot(6,1,5);
stem(received_bits_user_1,'filled');
title("Received/Despreaded bits user 1");

subplot(6,1,6);
stem(received_bits_user_2,'filled');
title("Received/Despreaded bits user 2");

disp("Error in Transmission user 1:");
disp(biterr(message_1,received_bits_user_1));

disp("Error in Tranmission user 2:");
disp(biterr(message_2,received_bits_user_2));

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
        temp = sum(signal(start:stop)); % Threshold is 0 as we have used PNRZ
        if(temp>0)
            arr(i+1) = 1;
        elseif(temp<0)
            arr(i+1) = -1;
        end
    end
end