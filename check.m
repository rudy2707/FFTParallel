% Axel Fahy - 01.12.2015
% Verification of matrix multiplication
% for parallele programming.

out;    % Import out.m

C = fft(A);

if C == B
    disp("CORRECT")
else
    disp("FALSE")
end

