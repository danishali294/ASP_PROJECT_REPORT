function Rate = Rate_Caculation(Channel,Est_Channel,NumofStream,No)

[U, ~, V]=svd(Est_Channel);
G = U(:,1:NumofStream)'*Channel*V(:,1:NumofStream);
Rate = abs(log2(det(eye(size(G,1))+1/(No*size(G,1))*(G*G'))));