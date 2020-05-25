function [DFT] = DFT_Codebook(Num_Antennas,Range)

Num_Antennas_Index=0:1:Num_Antennas-1; 

for g=1:1:Num_Antennas
    codebook(:,g)=sqrt(1/Num_Antennas)*exp(-1j*(2*pi)*Num_Antennas_Index*((g-1)/Num_Antennas));
end

DFT = codebook(:,Range);

end
