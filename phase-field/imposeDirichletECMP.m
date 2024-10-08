
function [ stiff, rhs, bcdof, bcval ] = imposeDirichletECMP(stiff, rhs, PHTelem, GIFTmesh, p, q, force)
%impose Dirichlet boundary conditions for elastic square
%fixed (homogeneous) boundary conditions on the left side
%displacement (bound_disp units) to the right on the right side
%supports multipatches. Assumes boundary conditions are imposed on the
%first and last patches
%将弹性方形固定(齐次)边界条件施加于左侧位移(bound_disp单位)到右侧支持多块的狄利克雷边界条件。
%假设边界条件施加在第一个和最后一个补丁

ngauss_edge = 2*q;
[gauss_weight_edge, gauss_coord_edge] = quadrature( ngauss_edge, 'GAUSS', 1 );

%for each neumann edge store the element index and orientation
%orientation: 1-down, 2-right, 3-up, 4-left
neumann_up =[];
neumann_down =[];
neumann_left =[];
neumann_right =[];

%detect which nodes have support on the left and right boundary
bcdof_down = [];
bcdof_up = [];
bcdof_1 = [];
bcdof_2 = [];
bcdof_3 = [];
bcdof_4 = [];
bcdof_5 = [];

%define side node indices
down_nodes = 1:p+1;
right_nodes = (p+1):(p+1):(p+1)*(q+1);
up_nodes = 1+(p+1)*q:(p+1)*(q+1);
left_nodes = 1:(p+1):(1+(p+1)*q);

numPatches = length(PHTelem);

for patchIndex=1:length(PHTelem)
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            
            if (patchIndex==3) && (isempty(PHTelem{patchIndex}(i).neighbor_up))
                 neumann_up = [neumann_up; i, 3, patchIndex];
            end
            
             if (patchIndex==4) && (isempty(PHTelem{patchIndex}(i).neighbor_up))
                 neumann_up = [neumann_up; i, 3, patchIndex];   
             end
            
             if (patchIndex==1) && (isempty(PHTelem{patchIndex}(i).neighbor_left))
                 neumann_left = [neumann_left; i, 4, patchIndex];  
             end
            
             if (patchIndex==4) && (isempty(PHTelem{patchIndex}(i).neighbor_left))
                 neumann_left = [neumann_left; i, 4, patchIndex];
             end
             
             if (patchIndex==2) && (isempty(PHTelem{patchIndex}(i).neighbor_right))
                 neumann_right = [neumann_right; i, 2, patchIndex];
               
             end
            
             if (patchIndex==3) && (isempty(PHTelem{patchIndex}(i).neighbor_right))
                 neumann_right = [neumann_right; i, 2, patchIndex];
             end
             
             
              if (patchIndex==1) && (isempty(PHTelem{patchIndex}(i).neighbor_down))
                 neumann_down = [neumann_down; i, 1, patchIndex];
             end
            
             if (patchIndex==2) && (isempty(PHTelem{patchIndex}(i).neighbor_down))
                 neumann_down = [neumann_down; i, 1, patchIndex];
             end
            
             %%%%该向量节点编号
             if (patchIndex==1) && (isempty(PHTelem{patchIndex}(i).neighbor_down))
                if (PHTelem{patchIndex}(i).vertex(2)==0) 
                    bcdof_1= PHTelem{patchIndex}(i).nodesGlobal(down_nodes);                    
                end
             end
             if (patchIndex==2) && (isempty(PHTelem{patchIndex}(i).neighbor_down))
                if (PHTelem{patchIndex}(i).vertex(2)==0) 
                    bcdof_1= [bcdof_1,PHTelem{patchIndex}(i).nodesGlobal(down_nodes)];                    
                end
             end
            
              if (patchIndex==3) && (isempty(PHTelem{patchIndex}(i).neighbor_up))
                if (PHTelem{patchIndex}(i).vertex(4)==1) 
                    bcdof_3 = PHTelem{patchIndex}(i).nodesGlobal(up_nodes);                    
                end
              end
              if (patchIndex==4) && (isempty(PHTelem{patchIndex}(i).neighbor_up))
                if (PHTelem{patchIndex}(i).vertex(4)==1) 
                    bcdof_3= [bcdof_3,PHTelem{patchIndex}(i).nodesGlobal(up_nodes)];                    
                end
              end
             if (patchIndex==1) && (isempty(PHTelem{patchIndex}(i).neighbor_up))
                if (PHTelem{patchIndex}(i).vertex(4)==1) 
                    bcdof_5 = PHTelem{patchIndex}(i).nodesGlobal(up_nodes);                    
                end
              end
              if (patchIndex==4) && (isempty(PHTelem{patchIndex}(i).neighbor_down))
                if (PHTelem{patchIndex}(i).vertex(2)==0) 
                    bcdof_5 = [bcdof_2,PHTelem{patchIndex}(i).nodesGlobal(down_nodes)];                    
                end
             end
                 
        end
    end
end



%remove duplicated entries
%take into account that there are 2 global indices (displacements) per node
bcdof_1_x = 2*bcdof_1-1;
bcdof_1_y = 2*bcdof_1;
bcdof_3_y = 2*bcdof_3;
bcdof_5_y = 2*bcdof_5 



%impose the prescribed displacement for each global index 对每个全球指数施加规定的位移
bcval_1_x = zeros(size(bcdof_1_x));
bcval_1_y = zeros(size(bcdof_1_y));
bcval_3_y = zeros(size(bcdof_3_y));
bcval_5_y = zeros(size(bcdof_5_y));

bcdof = [bcdof_1_x, bcdof_1_y, bcdof_3_y,bcdof_5_y];
bcval = [bcval_1_x, bcval_1_y, bcval_3_y,bcval_5_y];


%impose Neumann boundary conditons
neumann = [neumann_up; neumann_down; neumann_left; neumann_right];
num_neumann_elem = size(neumann,1);

for i_neu=1:num_neumann_elem
    
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
    dscrtx = reshape([2*scrtx-1; 2*scrtx],1,2* nument_edge);%单元边界节点编号
    localrhsed = zeros(2*nument_edge, 1);
    
    %loop over Gauss points and compute the integral ????????
    for igauss = 1:ngauss_edge
        
        %find the evaluation point, taking into account the edge type（映射单元的对应划分）
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
        
        %------------------------------------------------------------------------;
        %  computation of normal;
        if(corient==1)
            nor(1) = dxdxi(2,1);% dy/dxi
            nor(2) = -dxdxi(1,1);%dx/dxi
        elseif(corient==2)
            nor(1) = dxdxi(2,2);
            nor(2) = -dxdxi(1,2);
        elseif(corient==3)
            nor(1) = -dxdxi(2,1);
            nor(2) = dxdxi(1,1);
        else
            nor(1) = -dxdxi(2,2);
            nor(2) = dxdxi(1,2);
        end
        
        tmp = hypot(nor(1),nor(2));
        normal = nor/tmp; % normal vector in two dimensions
        
        
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
        
        [theta, r] = cart2pol(coord(1)-0.5,coord(2)-1);
        cracklength = 0.5;
        Kexact1=force*sqrt(pi*cracklength); Kexact2 = 0;
        exactStress=exact_stressesEC(r, theta, Kexact1, Kexact2);       
        
        taux = normal(1)*exactStress(1) + normal(2)*exactStress(3);
        tauy = normal(1)*exactStress(3) + normal(2)*exactStress(2);
                
        localrhsed(1:2:end-1) = localrhsed(1:2:end-1) + R.*taux.*scalefac.*gauss_weight_edge(igauss).*J;
        localrhsed(2:2:end) = localrhsed(2:2:end) + R.*tauy.*scalefac.*gauss_weight_edge(igauss).*J;
        
    end
    rhs(dscrtx)=rhs(dscrtx)+localrhsed;
end

[stiff,rhs]=feaplyc2sym(stiff, rhs, bcdof, bcval);
