%   (7; 4) Hamming code over BSC
%
%
%    Autor: WANG Yuqing
%    Last Modify:2021.11.23
%    Runtime:MATLAB(R) 2016a

clc,clear,close all



G=[1,0,0,1,0,1,1;
   0,1,0,1,0,1,0;
   0,0,1,1,0,0,1;
   0,0,0,0,1,1,1];

H=[1,1,1,1,0,0,0;
   1,1,0,0,1,1,0;
   1,0,1,0,1,0,1];
HT=H';

%You can modify the number of the sample here
no_of_samples = 1e4;

%generate transmit signals +1/-1 randomly
Tx_signal = sign(randn(no_of_samples,1));

for i = 1:no_of_samples
    if Tx_signal(i) == -1
        Tx_signal(i)=0;
    end
end

Encode = zeros(no_of_samples/4,7);

for i=1:no_of_samples/4
    temp1 = [Tx_signal((i-1)*4+1),Tx_signal((i-1)*4+2),Tx_signal((i-1)*4+3),Tx_signal((i-1)*4+4)];
    Encode(i,:) = temp1*G;
end

Encode = mod(Encode,2);
Encode_information_bits = [Encode(:,1),Encode(:,2),Encode(:,3),Encode(:,5)];

p=0:0.1:0.5;
[m,n]=size(p);
bit_error_rate_decoded = zeros(1,n);
information_bit_error_rate_decoded = zeros(1,n);
bit7_error_rate_decoded = zeros(1,n);

for a=1:n
    Rx_signal = Encode; 
    p_s=p(a);
    BSC = sign(rand(no_of_samples/4,7)-p_s);

    %pass all bits to AWGN channel

    for i = 1:no_of_samples/4
        for j = 1:7
            if BSC(i,j) == -1  %no change of state
                if Rx_signal(i,j) == 1
                    Rx_signal(i,j) = 0 ;
                else
                    Rx_signal(i,j) = 1 ;
                end
            end
        end
    end
    
    %calculate the raw bit error rate
    no_error_bit_raw = 0;
    for i=1:no_of_samples/4
        for j=1:7
            if Rx_signal(i,j) ~= Encode(i,j)
                no_error_bit_raw = no_error_bit_raw+1;
            end
        end
    end

    bit_error_rate_raw(a) = no_error_bit_raw/((no_of_samples/4)*7);
    
     %calculate the raw block code error rate
    
    no_error_7bit_raw = 0;
    for i=1:no_of_samples/4
        temp1 = Rx_signal(i,:);
        temp2 = Encode(i,:);
        for j=1:7     
            if temp1(j)~=temp2(j)
                 no_error_7bit_raw = no_error_7bit_raw+1;
                 break
            end
        end
    end

    bit7_error_rate_raw(a) = no_error_7bit_raw/(no_of_samples/4);
    

    Error_Syndrome = [];
    for i=1:no_of_samples/4
        temp1 = Rx_signal(i,:);
        Error_Syndrome(i,:) = temp1*HT;
    end

    Error_Syndrome = mod(Error_Syndrome,2);
    Error_Syndrome = num2str(Error_Syndrome);

    %find the error bit's location

    for i=1:no_of_samples/4    
        temp1 = Error_Syndrome(i,:);
        loc_wrong_bit(i) = 8-bin2dec(temp1);
        if loc_wrong_bit(i)==8
            loc_wrong_bit(i)=0;
        end
    end

    %correct the error
    Decode = Rx_signal;
    for i=1:no_of_samples/4

        if loc_wrong_bit(i) ~= 0
            if Decode(i,loc_wrong_bit(i)) == 1
                Decode(i,loc_wrong_bit(i)) = 0;
            else
                Decode(i,loc_wrong_bit(i)) = 1;
            end
        end
    end

    Decode_information_bits = [Decode(:,1),Decode(:,2),Decode(:,3),Decode(:,5)];
    
    %calculate the bit_error_rate
    no_error_bit = 0;
    for i=1:no_of_samples/4
        for j=1:7
            if Decode(i,j) ~= Encode(i,j)
                no_error_bit = no_error_bit+1;
            end
        end
    end

    bit_error_rate_decoded(a) = no_error_bit/((no_of_samples/4)*7);

    %calculate the 7bit_error_rate
    no_error_7bit = 0;
    for i=1:no_of_samples/4
        temp1 = Decode(i,:);
        temp2 = Encode(i,:);
        for j=1:7     
            if temp1(j)~=temp2(j)
                 no_error_7bit = no_error_7bit+1;
                 break
            end
        end
    end

    bit7_error_rate_decoded(a) = no_error_7bit/(no_of_samples/4);

    %calculate the information bit_error_rate
    information_bit_error = 0;
    for i=1:no_of_samples/4
        for j=1:4
            if Decode_information_bits(i,j) ~= Encode_information_bits(i,j)
                information_bit_error = information_bit_error+1;
            end
        end
    end
    information_bit_error_rate_decoded(a) = information_bit_error/no_of_samples;
end

figure
scatter(p,bit_error_rate_decoded,'b','LineWidth',2)
xlabel('ps'),ylabel('error rate'),title('error rate for (7,4)Hamming code in BSC')
hold on
plot(p,bit7_error_rate_decoded,'r','LineWidth',2)
plot(p,information_bit_error_rate_decoded,'g','LineWidth',2)
plot(p,bit_error_rate_raw,'y','LineWidth',2)
plot(p,bit7_error_rate_raw,'c','LineWidth',2)
legend('bit error rate decoded','bit7 error rate decoded','information bit error rate decoded','bit error rate raw','bit7 error rate raw','Location','NorthWest')
grid


%   end

