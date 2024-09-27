function [ ustiff, urhs, bcdof, bcval ] = imposeDirichletESNeuMP_u(ustiff, urhs, PHTelem, GIFTmesh, p, q, ij,u_inc)
%impose Dirichlet boundary conditions for elastic rectangle
%fixed (homogeneous) boundary conditions on the left side
%Neumann (traction) boundary conditions on the right side
%supports multipatches

ngauss_edge = q+2;
[gauss_weight_edge, gauss_coord_edge] = quadrature( ngauss_edge, 'GAUSS', 1 );

%numPatches = length(PHTelem);


%detect which nodes have support on the left and right boundary
bcdof_left = [];
bcdof_right = [];
%bcdof_up = [];
bcdof_down = [];


%for each neumann edge store the element index and orientation
%orientation: 1-bottom, 2-right, 3-top, 4-left
neumann_up = [];
%neumann_down = [];
%neumann_right = [];
%neumann_left = [];

%define side node indices
down_nodes = 1:p+1;
right_nodes = (p+1):(p+1):(p+1)*(q+1);
up_nodes = 1+(p+1)*q:(p+1)*(q+1);
left_nodes = 1:(p+1):(1+(p+1)*q);

%find the nodes corresponding to the Dirichlet boundary on patch 1
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
           
           %if (patchIndex==3) && (PHTelem{patchIndex}(i).vertex(4)==1) 
               %bcdof_up=[bcdof_up,PHTelem{patchIndex}(i).nodesGlobal(up_nodes)];
           %end
           %if (patchIndex==4) && (PHTelem{patchIndex}(i).vertex(4)==1)
              % bcdof_up=[bcdof_up,PHTelem{patchIndex}(i).nodesGlobal(up_nodes)];
          % end
        end
    end
end

%find the elements corresponding to the Neumann boundary on the last patch
for patchIndex=1:length(PHTelem)
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
           if (patchIndex==3) && (PHTelem{patchIndex}(i).vertex(4)==1) 
               neumann_up = [neumann_up; i, 3,patchIndex];
           end
           if (patchIndex==4) && (PHTelem{patchIndex}(i).vertex(4)==1)
               neumann_up = [neumann_up; i, 3,patchIndex];
           end
        end    
    end
end

%remove duplicated entries
%bcdof_right = unique(bcdof_right);
bcdof_left = unique(bcdof_left);
bcdof_right = unique(bcdof_right);
bcdof_down = unique(bcdof_down);

bcdof_left_y = 2*bcdof_left;
bcdof_right_y = 2*bcdof_right;
bcdof_down_x = 2*bcdof_down-1;
bcdof_down_y = 2*bcdof_down;

%impose the prescribed displacement for each global index
bcval_left_y = zeros(size(bcdof_left_y));
bcval_right_y = zeros(size(bcdof_right_y));
bcval_down_x = zeros(size(bcdof_down_x));
bcval_down_y = zeros(size(bcdof_down_y));

bcdof = [bcdof_down_x, bcdof_down_y,bcdof_left_y,bcdof_right_y];
bcval = [bcval_down_x, bcval_down_y,bcval_left_y,bcval_right_y];


%impose Neumann boundary conditons
neumann = neumann_up;
for i_neu=1:size(neumann,1)
    
    i = neumann(i_neu,1);
    corient = neumann(i_neu, 2);
    patchIndex = neumann(i_neu, 3);
    
   xmin = PHTelem{patchIndex}(i).vertex(1);
    xmax = PHTelem{patchIndex}(i).vertex(3);
    ymin = PHTelem{patchIndex}(i).vertex(2);
    ymax = PHTelem{patchIndex}(i).vertex(4);
    
    if (corient == 1) || (corient == 3)
        scalefac = (xmax-xmin)/2;
    else
        scalefac = (ymax-ymin)/2;
    end
    
    nument = size(PHTelem{patchIndex}(i).C,1);
    
    
    scrtx = PHTelem{patchIndex}(i).nodesGlobal(1:nument);
    
    if corient==1
        scrtx = scrtx(down_nodes);
        nument_edge = p+1;
    elseif corient==2
        scrtx = scrtx(right_nodes);
        nument_edge = q+1;
    elseif corient==3
        scrtx = scrtx(up_nodes);
        nument_edge = p+1;
    else
        scrtx = scrtx(left_nodes);
        nument_edge = q+1;
    end
    dscrtx = reshape([2*scrtx-1; 2*scrtx],1,2* nument_edge);
    localrhsed = zeros(2*nument_edge, 1);
    
    %loop over Gauss points and compute the integral
    for igauss = 1:ngauss_edge
        
        %find the evaluation point, taking into account the edge type
        if corient==1
            v_hat = -1;
            u_hat = gauss_coord_edge(igauss);
        elseif corient==2
            u_hat = 1;
            v_hat = gauss_coord_edge(igauss);
        elseif corient==3
            v_hat = 1;
            u_hat = gauss_coord_edge(igauss);
        else
            u_hat = -1;
            v_hat = gauss_coord_edge(igauss);
        end
        
        %evaluate the basis functions
        R = phtBasis(u_hat, v_hat, PHTelem{patchIndex}(i).C, p, q);
        
        %evaluate the derivatives of the mapping from parameter
        %space to physical space
        [coord, dxdxi] = paramMap( GIFTmesh{patchIndex}, u_hat, v_hat, xmin, ymin, xmax, ymax);
        dxdxi = dxdxi';
        %  jacobian of edge mapping;
        if((corient==1)||(corient==3))
            J = hypot(dxdxi(1,1), dxdxi(2,1));
        else
            J = hypot(dxdxi(1,2), dxdxi(2,2));
        end
        
        
        %consider on the values corresponding to the shape functions on the
        %boundary
        if corient==1
            R = R(down_nodes)';
        elseif corient==2
            R = R(right_nodes)';
        elseif corient==3
            R = R(up_nodes)';
        else
            R = R(left_nodes)';
        end
        
        %assume traction is constant
        taux = u_inc*ij*ones(nument_edge,1);
        tauy = 0*ones(nument_edge,1);
        
        localrhsed(1:2:end-1) = localrhsed(1:2:end-1) + R.*taux.*scalefac.*gauss_weight_edge(igauss).*J;
        localrhsed(2:2:end) = localrhsed(2:2:end) + R.*tauy.*scalefac.*gauss_weight_edge(igauss).*J;
        
    end
    urhs(dscrtx)=urhs(dscrtx)+localrhsed;
end

%update the stiffness matrix and rhs with Dirichlet boundary values
[ustiff,urhs]=feaplyc2(ustiff, urhs, bcdof, bcval);

