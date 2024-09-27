function [Fload,Uload,u]=solve_u(ustiff ,u, ufree,uload,ij,u_inc,Fload,Uload)
%u(uload) = u_inc*ij;
u(uload) = sum(u_inc(1:ij));
u(ufree) =  -ustiff(ufree,ufree)\(ustiff(ufree,uload)*u(uload)); 
Fload(ij+1) = sum(ustiff(uload,:)*u);
%Uload(ij+1) = u_inc*ij;
Uload(ij+1) = sum(u_inc(1:ij));
end

