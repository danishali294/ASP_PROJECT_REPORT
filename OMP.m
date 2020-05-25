function [supp] = OMP(Y,A,s,Tol)
% Simulated Orthogona matching pursuit

% Normalize the columns of A
C = sum(abs(A).^2,1).^0.5; 
A = A./C(ones(size(A,1),1),:);

[m,n] = size(A);

supp = [];
supp_c = 1:n;
res = Y;

if nargin < 4
    while length(supp) < s
       obs = sum(abs((res)'*A(:,supp_c)).^2,1);
       [mval,midx] = max(obs);
       supp = [supp supp_c(midx)];
       supp_c(midx) = [];
       [tmpU,R] = qr(A(:,supp),0);
       tmpPperp = eye(m) - tmpU*tmpU';
       res = tmpPperp*Y;
    end
    
    supp = sort(supp(1:end),'ascend');
else
    norm_now = norm(res,2)^2 ;
    while  length(supp) < 0.8*m 
       norm_past = norm_now;
       obs = sum(abs((res)'*A(:,supp_c)).^2,1);
       [mval,midx] = max(obs);
       supp = [supp supp_c(midx)];
       supp_c(midx) = [];
       [tmpU,R] = qr(A(:,supp),0);
       tmpPperp = eye(m) - tmpU*tmpU';
       res = tmpPperp*Y;
       
       norm_now = norm(res,2)^2;
       if norm_past - norm_now < Tol 
          supp = sort(supp(1:end-1),'ascend');
          break;
       end
    end
end

