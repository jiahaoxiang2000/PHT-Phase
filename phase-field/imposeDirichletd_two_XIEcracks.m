function [ d ] = imposeDirichletd_two_XIEcracks(dstiff, drhs, PHTelem, p, q,sizeBasis)
%impose Dirichlet boundary conditions for elastic square
%fixed (homogeneous) boundary conditions on the left side
%displacement (bound_disp units) to the right on the right side
%supports multipatches. Assumes boundary conditions are imposed on the
%first and last patches

%detect which nodes have support on the left and right boundary
bcdof_mid = [];
%define side node indices
down_nodes = 1:p+1;
right_nodes = (p+1):(p+1):(p+1)*(q+1);
up_nodes = 1+(p+1)*q:(p+1)*(q+1);
left_nodes = 1:(p+1):(1+(p+1)*q);


for patchIndex=1:length(PHTelem)
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            
          if (patchIndex==6) && (PHTelem{patchIndex}(i).vertex(3)==0.5)  && (PHTelem{patchIndex}(i).vertex(4)<=0.5) 
               bcdof_mid=[bcdof_mid,PHTelem{patchIndex}(i).nodesGlobal(right_nodes)];
          end
          if (patchIndex==6) && (PHTelem{patchIndex}(i).vertex(1)==0.5)   && (PHTelem{patchIndex}(i).vertex(4)<=0.5) 
               bcdof_mid=[bcdof_mid,PHTelem{patchIndex}(i).nodesGlobal(left_nodes)];
            
          end
          if (patchIndex==11) && (PHTelem{patchIndex}(i).vertex(3)==0.5)  && (PHTelem{patchIndex}(i).vertex(2)>=0.5) 
               bcdof_mid=[bcdof_mid,PHTelem{patchIndex}(i).nodesGlobal(right_nodes)];
          end
          if (patchIndex==11) && (PHTelem{patchIndex}(i).vertex(1)==0.5)  && (PHTelem{patchIndex}(i).vertex(2)>=0.5) 
               bcdof_mid=[bcdof_mid,PHTelem{patchIndex}(i).nodesGlobal(left_nodes)];
          end
        end
    end
end

%remove duplicated entries

bcdof_mid = unique(bcdof_mid);
%take into account that there are 1 global indices (d) per node
d=zeros(sizeBasis,1);
dalldofs=(1:sizeBasis);
dfixed=bcdof_mid;
dfree=setdiff(dalldofs,dfixed);
d(dfixed)=1;

[d]=solve_d(dstiff, drhs, d, dfree, dfixed);