function Hn = Hn_evaluation(u,Hn,Bu,edofMat,lambda,mu,nel)
epsilon_gauss=zeros(3,16*nel); %高斯点上的应变值
U=u(edofMat);
for iel=1:nel
    epsilon_gauss(:,16*(iel-1)+1:16*iel)=reshape(Bu(:,32*(iel-1)+1:32*iel)*U(iel,:)',3,16);
end
eps_tr_plus = operator_trace_plus(epsilon_gauss);
eps_plus=zeros(3,16*nel);
parfor i = 1:size(eps_plus,2)
  eps_plus(:,i) = operator_plus(epsilon_gauss(:,i)); %分解出的正应变计算
end
Psi_plus = (lambda/2).*eps_tr_plus.^2 + mu.*(eps_plus(1,:).^2+eps_plus(2,:).^2+2*((eps_plus(3,:)/2).^2));%分解出的正存储能量计算
%Psi = (abs(Psi_plus-gc)+(Psi_plus-gc))/2;
Hn = max(reshape(Psi_plus,16,nel),Hn);