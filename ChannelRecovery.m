function H_Est = ChannelRecovery (y,eq_mea_mtx,AbG,AmG,sparsity,stop_thr)

supp = OMP(y,eq_mea_mtx,sparsity,stop_thr);

X_Est          = zeros(size(eq_mea_mtx,2),1);
X_Est(supp,1)  = (eq_mea_mtx(:,supp))\y; 

Ha_Est = reshape(X_Est,size(AmG,2),size(AbG,2));
H_Est  = AmG*Ha_Est*AbG';

end