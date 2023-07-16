            %% Section 2
load('hw3-2.mat');
[M,T,N] = size(data);
estimated_freqs = zeros(size(label));
for i=1:N
   Y = data(:,:,i);
   f = CCA(Y,freq);
   estimated_freqs(i) = f;
end

correct_estimation_num = sum(label==estimated_freqs);
accuracy = correct_estimation_num/N;
disp("Accuracy:");
disp(accuracy);

    %% Local Necessary Functions

function X = Harmonic_Template(f)
    fmax = 40;
    T = 1250;
    fs = 250;
    t = 0: 1/fs: T/fs - 1/fs;
    K = 2 * ceil(fmax/f);
    X = zeros(K,T);
    for i=1:K/2
       ftemp = f*i;
       temp_signal = [sin(2*pi*ftemp*t);cos(2*pi*ftemp*t)];
       X(2*i-1:2*i,:) = temp_signal;
    end
end

function u = GEVD(A,B)
    C = inv(B)*A;
    [U,D] = eig(C);
    [~,i] = max(max(D));
    u = U(:,i);
end

function [a,b] = ab_Calc(Y,f)
    X = Harmonic_Template(f);
    Rx = X*transpose(X);
    Ry = Y*transpose(Y);
    Rxy = X*transpose(Y);
    Ryx = Y*transpose(X);
    [U,D] = eig(Rx);
    Rx_inv_half = inv(U*(D.^(1/2))*transpose(U));
    [U,D] = eig(Ry);
    Ry_inv_half = inv(U*(D.^(1/2))*transpose(U));
    Rx_inv = inv(Rx);
    Ry_inv = inv(Ry);
    S1 = Rx_inv_half * Rxy * Ry_inv * Ryx * Rx_inv_half;
    S2 = Ry_inv_half * Ryx * Rx_inv * Rxy * Ry_inv_half;
    c = GEVD(S1,eye(length(S1)));
    d = GEVD(S2,eye(length(S2)));
    a = Rx_inv_half * c;
    b = Ry_inv_half * d;
end

function p = corr_val(a,b)
    p = (transpose(a)*b)/(norm(a)*norm(b));
end

function f_e = CCA(Y,freq)
    Corr_Values = zeros(length(freq),1);
    for i=1:length(freq)
        f = freq(i);
        X = Harmonic_Template(f);
        [a,b] = ab_Calc(Y,f);
        z1_T = transpose(a)*X;
        z2_T = transpose(b)*Y;
        p = corr_val(transpose(z1_T),transpose(z2_T));
        Corr_Values(i) = p;
    end
    [~,idx] = max(Corr_Values);
    f_e = freq(idx);
end