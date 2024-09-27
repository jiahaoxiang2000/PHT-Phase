function u = precomutation(u,nel,Cmat,loaddofs,freedofs,Bu,iK,jK,detJacobi)
%stiff = (stiff+stiff')/2
KK=zeros(32,32*nel);
for iel=1:nel
    KK(:,32*(iel-1)+1:32*iel)=Bu(:,32*(iel-1)+1:32*iel)'*blkdiag(Cmat,Cmat,Cmat,Cmat,Cmat,Cmat,Cmat,Cmat,Cmat,Cmat,Cmat,Cmat,Cmat,...
       Cmat,Cmat,Cmat)*(Bu(:,32*(iel-1)+1:32*iel).*kron(detJacobi(:,iel),ones(3,32)));%逗号之前表示取所有行，逗号之后表示取a到b的数
end
sK = reshape(KK,32*32*nel,1);
K = sparse(iK,jK,sK); K = (K+K')/2;
u(freedofs) = -K(freedofs,freedofs)\(K(freedofs,loaddofs)*u(loaddofs));
%u(freedofs) = -stiff(freedofs,freedofs)\(stiff(freedofs,loaddofs)*u(loaddofs));