function [Hn,Hn1] = Hn_evaluation2(u,Hn1,Bu,lambda,mu,dsctrx,ii,jj,elementCounter)
%% history functional evaluation
epsilon_gauss = reshape(Bu'*u(dsctrx),3,1);
eps_tr_plus = operator_trace_plus(epsilon_gauss);
eps_plus = zeros(3,1);
  eps_plus(:,1) = operator_plus(epsilon_gauss(:,1));
Psi_plus = (lambda/2).*eps_tr_plus.^2 + mu.*(eps_plus(1,:).^2+eps_plus(2,:).^2+0.5*eps_plus(3,:).^2); % 
if Psi_plus>Hn1(ii+(jj-1)*4,elementCounter)
   Hn1(ii+(jj-1)*4,elementCounter) = Psi_plus;
   Hn = Psi_plus;
else
    Hn = Hn1(ii+(jj-1)*4,elementCounter);
end

% evaluation of the trace of the positive epsilon
function eps_tr_plus = operator_trace_plus(eps)
eps_tr_plus = ((eps(1,:) + eps(2,:)) + abs(eps(1,:)+eps(2,:)))/2;

% evaluation of the positive epsilon
function eps_plus = operator_plus(eps)
[V,D] = eig([eps(1),eps(3)/2; eps(3)/2,eps(2)]);
eps_plus = (D(1,1)+abs(D(1,1)))/2*kron(V(:,1),V(:,1)')...
  + (D(2,2)+abs(D(2,2)))/2*kron(V(:,2),V(:,2)');
eps_plus = [eps_plus(1,1); eps_plus(2,2); 2*eps_plus(1,2)];