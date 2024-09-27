function [ ustiff, urhs, bcdof, bcval,bcdof_up_x  ] = imposeDirichletESMP_u(ustiff, urhs, PHTelem, p, q, ij,u_inc)
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
           if (patchIndex==4) && (PHTelem{patchIndex}(i).vertex(1)==0) 
               bcdof_left=[bcdof_left,PHTelem{patchIndex}(i).nodesGlobal(left_nodes)];
           end
          
          
           
           if  (patchIndex==1) && (PHTelem{patchIndex}(i).vertex(2)==0)
               bcdof_down=[bcdof_down,PHTelem{patchIndex}(i).nodesGlobal(down_nodes)];
           end
           if (patchIndex==2) && (PHTelem{patchIndex}(i).vertex(2)==0)
               bcdof_down=[bcdof_down,PHTelem{patchIndex}(i).nodesGlobal(down_nodes)];
           end

          
           
           if (patchIndex==3) && (PHTelem{patchIndex}(i).vertex(3)==1)
               bcdof_right=[bcdof_right,PHTelem{patchIndex}(i).nodesGlobal(right_nodes)];
           end
           if (patchIndex==2) && (PHTelem{patchIndex}(i).vertex(3)==1) 
               bcdof_right=[bcdof_right,PHTelem{patchIndex}(i).nodesGlobal(right_nodes)];
           end
          
           
           if (patchIndex==3) && (PHTelem{patchIndex}(i).vertex(4)==1) 
               bcdof_up=[bcdof_up,PHTelem{patchIndex}(i).nodesGlobal(up_nodes)];
           end
           if (patchIndex==4) && (PHTelem{patchIndex}(i).vertex(4)==1)
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
%take into account that there are 2 global indices (d) per node
%take into account that there are 2 global indices (displacements) per node
bcdof_down_x = 2*bcdof_down;
bcdof_down_y = 2*bcdof_down;

bcdof_right_y = 2*bcdof_right;
bcdof_left_y = 2*bcdof_left;
bcdof_up_x = 2*bcdof_up-1;
bcdof_up_y = 2*bcdof_up;

%impose the prescribed displacement for each global index
bcval_down_x = zeros(size(bcdof_down_x));
bcval_down_y = zeros(size(bcdof_down_y));
bcval_right_y = zeros(size(bcdof_right_y));
bcval_left_y = zeros(size(bcdof_left_y));
bcval_up_x = ij*u_inc*ones(size(bcdof_up_x));
bcval_up_y = zeros(size(bcdof_up_y));

bcdof = [bcdof_down_x, bcdof_down_y, bcdof_up_x, bcdof_up_y,bcdof_left_y,bcdof_right_y];
bcval = [bcval_down_x, bcval_down_y, bcval_up_x, bcval_up_y,bcval_left_y,bcval_right_y];
[ustiff,urhs]=feaplyc2sym(ustiff, urhs, bcdof, bcval);

