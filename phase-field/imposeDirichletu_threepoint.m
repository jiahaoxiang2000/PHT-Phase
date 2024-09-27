function [u,Fload,Uload ] = imposeDirichletu_threepoint(ustiff, PHTelem, p, q,sizeBasis,ij,u_inc,Fload,Uload,u)
%impose Dirichlet boundary conditions for elastic square
%fixed (homogeneous) boundary conditions on the left side
%displacement (bound_disp units) to the right on the right side
%supports multipatches. Assumes boundary conditions are imposed on the
%first and last patches

%detect which nodes have support on the left and right boundary
bcdof_mid = [];
bcdof_mid1 = [];
bcdof_mid2 = [];

%define side node indices
down_nodes = 1:p+1;
right_nodes = (p+1):(p+1):(p+1)*(q+1);
up_nodes = 1+(p+1)*q:(p+1)*(q+1);
left_nodes = 1:(p+1):(1+(p+1)*q);

for patchIndex=1:length(PHTelem)
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            
           if (patchIndex==90) && (PHTelem{patchIndex}(i).vertex(3)==1)  && (PHTelem{patchIndex}(i).vertex(4)==1)
               bcdof_mid=[bcdof_mid,PHTelem{patchIndex}(i).nodesGlobal(16)];
           end
           if (patchIndex==91) && (PHTelem{patchIndex}(i).vertex(1)==0)  && (PHTelem{patchIndex}(i).vertex(4)==1) 
               bcdof_mid=[bcdof_mid,PHTelem{patchIndex}(i).nodesGlobal(13)];
           end
           
           if (patchIndex==20) && (PHTelem{patchIndex}(i).vertex(3)==1) && (PHTelem{patchIndex}(i).vertex(2)==0) 
               bcdof_mid1=[bcdof_mid1,PHTelem{patchIndex}(i).nodesGlobal(4)];
           end
           
           if (patchIndex==1) && (PHTelem{patchIndex}(i).vertex(1)==0) && (PHTelem{patchIndex}(i).vertex(2)==0) 
               bcdof_mid2=[bcdof_mid2,PHTelem{patchIndex}(i).nodesGlobal(1)];
           end
           

        end
    end
end

%remove duplicated entries
bcdof_mid = unique(bcdof_mid);
bcdof_mid1 = unique(bcdof_mid1);
bcdof_mid2 = unique(bcdof_mid2);
%bcdof_mid = unique(bcdof_mid);
%take into account that there are 2 global indices (d) per node
alldofs = (1:2*sizeBasis);
loaddofs = 2*bcdof_mid;
fixeddofs1=2*bcdof_mid2-1;
fixeddofs2=2*bcdof_mid2;
%fixeddofs3=2*bcdof_right;
%fixeddofs4=2*bcdof_left;
fixeddofs12=union(fixeddofs1,fixeddofs2);%f1和f2的并集
%fixeddofs34=union(fixeddofs3,fixeddofs4);
%fixeddofs1234=union(fixeddofs12,fixeddofs34);
%fixeddofs5=2*bcdof_mid;
fixeddofs3=2*bcdof_mid1;
%fixeddofs12345=union(fixeddofs5,fixeddofs1234);
fixeddofs=union(fixeddofs3,fixeddofs12);
%loaddofs = 2*bcdof_up;
%fixeddofs1 = 2*bcdof_left-1;
%fixeddofs2 = 2*bcdof_right-1;
%fixeddofs3 = 2*bcdof_down-1;
%fixeddofs4 = 2*bcdof_down;
%fixeddofs5 = 2*bcdof_up-1;
%fixeddofs12 = union(fixeddofs1,fixeddofs2);
%fixeddofs34 = union(fixeddofs3,fixeddofs4);
%fixeddofs125 = union(fixeddofs12,fixeddofs5);
%fixeddofs = union(fixeddofs125,fixeddofs34);


uload = loaddofs(:);
ufixed = fixeddofs(:);
ufree = setdiff(alldofs,union(ufixed,uload));

[Fload,Uload,u]=solve_u(ustiff ,u, ufree,uload,ij,u_inc,Fload,Uload);