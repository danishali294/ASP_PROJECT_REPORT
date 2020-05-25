function [Channel,Abh,Amh,alpha] = Channel_Generation(Num_BS_Antennas,Num_MS_Antennas,Num_clusters,Num_rays,Sector)

BSAntennas_Index=0:1:Num_BS_Antennas-1; % Indices of the BS Antennas
MSAntennas_Index=0:1:Num_MS_Antennas-1; % Indices of the MS Antennas

Channel=zeros(Num_MS_Antennas,Num_BS_Antennas);


AoD = (2*rand(1,Num_clusters)/(Sector+0.2)) - 1/(Sector+0.2);
AoA = rand(1,Num_clusters);

for  n=1:1:Num_clusters
   Abh(:,n) = sqrt(1/Num_BS_Antennas)*exp(-1j*(2*pi)*BSAntennas_Index*AoD(n));
   Amh(:,n) = sqrt(1/Num_MS_Antennas)*exp(-1j*(2*pi)*MSAntennas_Index*AoA(n));
   alpha(n)    = (sqrt(1/2)*(randn(1,Num_rays)+1j*randn(1,Num_rays)));  
   Channel = Channel+ alpha(n)*Amh(:,n)*Abh(:,n)';
end

Channel = sqrt((Num_BS_Antennas*Num_MS_Antennas)/(Num_clusters*Num_rays))*Channel;
end
