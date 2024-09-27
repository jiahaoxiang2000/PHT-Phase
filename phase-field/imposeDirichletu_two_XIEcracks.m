function [u,Fload,Uload ] = imposeDirichletu_two_XIEcracks(ustiff, PHTelem, p, q,sizeBasis,ij,u_inc,Fload,Uload,u)
%impose Dirichlet boundary conditions for elastic square
%fixed (homogeneous) boundary conditions on the left side
%displacement (bound_disp units) to the right on the right side
%supports multipatches. Assumes boundary conditions are imposed on the
%first and last patches

%detect which nodes have support on the left and right boundary
bcdof_left = [];
bcdof_up = [];
bcdof_down = [];
bcdof_right = [];
%define side node indices
down_nodes = 1:p+1;
right_nodes = (p+1):(p+1):(p+1)*(q+1);
up_nodes = 1+(p+1)*q:(p+1)*(q+1);
left_nodes = 1:(p+1):(1+(p+1)*q);

for patchIndex=1:length(PHTelem)
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            
           if (patchIndex==1) && (PHTelem{patchIndex}(i).vertex(1)==0) 
               bcdof_left=[bcdof_left,PHTelem{patchIndex}(i).nodesGlobal(left_nodes)];
           end
           if (patchIndex==5) && (PHTelem{patchIndex}(i).vertex(1)==0) 
               bcdof_left=[bcdof_left,PHTelem{patchIndex}(i).nodesGlobal(left_nodes)];
           end
           if (patchIndex==9) && (PHTelem{patchIndex}(i).vertex(1)==0) 
               bcdof_left=[bcdof_left,PHTelem{patchIndex}(i).nodesGlobal(left_nodes)];
           end
           if (patchIndex==13) && (PHTelem{patchIndex}(i).vertex(1)==0) 
               bcdof_left=[bcdof_left,PHTelem{patchIndex}(i).nodesGlobal(left_nodes)];
           end
           
           
           
           if  (patchIndex==1) && (PHTelem{patchIndex}(i).vertex(2)==0)
               bcdof_down=[bcdof_down,PHTelem{patchIndex}(i).nodesGlobal(down_nodes)];
           end
           if (patchIndex==2) && (PHTelem{patchIndex}(i).vertex(2)==0)
               bcdof_down=[bcdof_down,PHTelem{patchIndex}(i).nodesGlobal(down_nodes)];
           end
           if  (patchIndex==3) && (PHTelem{patchIndex}(i).vertex(2)==0)
               bcdof_down=[bcdof_down,PHTelem{patchIndex}(i).nodesGlobal(down_nodes)];
           end
           if (patchIndex==4) && (PHTelem{patchIndex}(i).vertex(2)==0)
               bcdof_down=[bcdof_down,PHTelem{patchIndex}(i).nodesGlobal(down_nodes)];
           end

          
           
           if (patchIndex==4) && (PHTelem{patchIndex}(i).vertex(3)==1)
               bcdof_right=[bcdof_right,PHTelem{patchIndex}(i).nodesGlobal(right_nodes)];
           end
           if (patchIndex==8) && (PHTelem{patchIndex}(i).vertex(3)==1) 
               bcdof_right=[bcdof_right,PHTelem{patchIndex}(i).nodesGlobal(right_nodes)];
           end
           if (patchIndex==12) && (PHTelem{patchIndex}(i).vertex(3)==1)
               bcdof_right=[bcdof_right,PHTelem{patchIndex}(i).nodesGlobal(right_nodes)];
           end
           if (patchIndex==16) && (PHTelem{patchIndex}(i).vertex(3)==1) 
               bcdof_right=[bcdof_right,PHTelem{patchIndex}(i).nodesGlobal(right_nodes)];
           end
          
           
           if (patchIndex==13) && (PHTelem{patchIndex}(i).vertex(4)==1) 
               bcdof_up=[bcdof_up,PHTelem{patchIndex}(i).nodesGlobal(up_nodes)];
           end
           if (patchIndex==14) && (PHTelem{patchIndex}(i).vertex(4)==1)
               bcdof_up=[bcdof_up,PHTelem{patchIndex}(i).nodesGlobal(up_nodes)];
           end
            if (patchIndex==15) && (PHTelem{patchIndex}(i).vertex(4)==1) 
               bcdof_up=[bcdof_up,PHTelem{patchIndex}(i).nodesGlobal(up_nodes)];
           end
           if (patchIndex==16) && (PHTelem{patchIndex}(i).vertex(4)==1)
               bcdof_up=[bcdof_up,PHTelem{patchIndex}(i).nodesGlobal(up_nodes)];
           end
          
        end
    end
end

%remove duplicated entries
bcdof_right = unique(bcdof_right);
bcdof_down = unique(bcdof_down);
bcdof_up = unique(bcdof_up);
bcdof_left = unique(bcdof_left);
%bcdof_mid = unique(bcdof_mid);
%take into account that there are 2 global indices (d) per node
alldofs = (1:2*sizeBasis);
loaddofs = 2*bcdof_up;
fixeddofs1=2*bcdof_down-1;
fixeddofs2=2*bcdof_down;
%fixeddofs3=2*bcdof_right;
%fixeddofs4=2*bcdof_left;
fixeddofs12=union(fixeddofs1,fixeddofs2);%f1和f2的并集
%fixeddofs34=union(fixeddofs3,fixeddofs4);
%fixeddofs1234=union(fixeddofs12,fixeddofs34);
%fixeddofs5=2*bcdof_mid;
fixeddofs8=2*bcdof_up-1;
%fixeddofs12345=union(fixeddofs5,fixeddofs1234);
fixeddofs=union(fixeddofs8,fixeddofs12);
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
