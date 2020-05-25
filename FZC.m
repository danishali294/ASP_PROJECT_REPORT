function [ FZC ] = FZC(Num_BS_Antennas,Root)

if (nargin ~= 2 && Num_BS_Antennas == 256) 
    PRIME = [3 5 7 11 13 17 19 23 29 31 37 41 43 47 53 59 61];
    Root  = randsample(PRIME,1);
elseif (nargin ~= 2 && Num_BS_Antennas == 128) 
    PRIME = [3 5 7 11 13 17 19 23 29 31 37 41 43 47 53 59 61];
    Root  = randsample(PRIME,1);
elseif (nargin ~= 2 && Num_BS_Antennas == 64) 
    PRIME = [3 5 7 11 13 17 19 23 29 31 37 41 43 47 53 59 61];
    Root  = randsample(PRIME,1);
elseif(nargin ~= 2 && Num_BS_Antennas == 32) 
    PRIME = [3 5 7 11 13 17 19 23 29 31];
    Root  = randsample(PRIME,1);
elseif(nargin ~= 2 && Num_BS_Antennas == 16) 
    PRIME = [3 5 7 11 13];
    Root  = randsample(PRIME,1);
end

FZC = zeros(Num_BS_Antennas,1);
for n=0:Num_BS_Antennas-1
    FZC(n+1,1)=exp(-1j*pi*Root*n^2/Num_BS_Antennas);
end

end