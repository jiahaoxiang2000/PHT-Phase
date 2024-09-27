function [Nd,Nu,Ndp,Bd,Bu,Bdp] = shape_function(PHTelem,GIFTmesh,p,q,edofMat,dedofMat)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
nel=size(edofMat,1);
%Gauss points
if p>=3
    ngauss_x = p+1;
    ngauss_y = q+1;
else
    ngauss_x = p+2;
    ngauss_y = p+2;
end
[gauss_weight_x, gauss_coord_x] = quadrature( ngauss_x, 'GAUSS', 1 );
[gauss_weight_y, gauss_coord_y] = quadrature( ngauss_y, 'GAUSS', 1 );

%take the transpose so that they are in the format expected by
%bernstein_basis
gauss_coord_x = gauss_coord_x';
gauss_coord_y = gauss_coord_y';

%1D bernstein polynomials evaluated at the Gauss points on the master element
[B_u, dB_u] = bernstein_basis(gauss_coord_x,p);
[B_v, dB_v] = bernstein_basis(gauss_coord_y,q);

dBdu = zeros(ngauss_x, ngauss_y, (p+1)*(q+1));
dBdv = zeros(ngauss_x, ngauss_y, (p+1)*(q+1));
B = zeros(ngauss_x, ngauss_y, (p+1)*(q+1));


%the derivatives of the 2D Bernstein polynomials at Gauss points on the
%master element
basisCounter = 0;
for j=1:q+1
    for i=1:p+1
        basisCounter = basisCounter + 1;
        B(:,:,basisCounter) = B_u(:,i)*B_v(:,j)';
        dBdu(:,:,basisCounter) = dB_u(:,i)*B_v(:,j)';
        dBdv(:,:,basisCounter) = B_u(:,i)*dB_v(:,j)';
    end
end

nument=size(dedofMat,2);
%allocate memory for the triplet arrays

elementCounter = 0;
Nd = zeros(ngauss_x*ngauss_y,nument*nel);
Bd = zeros(2*ngauss_x*ngauss_y,nument*nel);
Nu = zeros(2*ngauss_x*ngauss_y,2*nument*nel);
Bu = zeros(3*ngauss_x*ngauss_y,2*nument*nel);

for patchIndex = 1:length(PHTelem)
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            elementCounter = elementCounter + 1;
            xmin = PHTelem{patchIndex}(i).vertex(1);
            xmax = PHTelem{patchIndex}(i).vertex(3);
            ymin = PHTelem{patchIndex}(i).vertex(2);
            ymax = PHTelem{patchIndex}(i).vertex(4);
            
            %the jacobian of the transformation from [-1,1]x[-1,1] to
            %[xmin, xmax]x [ymin, ymax]
            scalefac = (xmax - xmin)*(ymax - ymin)/4;
            nument = size(PHTelem{patchIndex}(i).C,1);
        
             for jj=1:ngauss_y
                 for ii=1:ngauss_x
              
                    
                    %evaluate the derivatives of the mapping from parameter
                    %space to physical space 求从参数空间到物理空间的映射的导数
                    
                    [~, dxdxi] = paramMap( GIFTmesh{patchIndex}, gauss_coord_x(ii), gauss_coord_y(jj), xmin, ymin, xmax, ymax);
                    J = det(dxdxi);
                   
                    R = (PHTelem{patchIndex}(i).C)*squeeze(B(ii,jj,:));
                   
                    dRdx = (PHTelem{patchIndex}(i).C)*squeeze(dBdu(ii,jj,:));
                    dRdy = (PHTelem{patchIndex}(i).C)*squeeze(dBdv(ii,jj,:));
                    
                    %multiply by the jacobian of the transformation from reference
                    %space to the parameter space 乘以从参考空间到参数空间的变换的雅可比矩阵
                    
                    dRdx = dRdx*2/(xmax-xmin); 
                    dRdy = dRdy*2/(ymax-ymin);
                    
                    % Solve for first derivatives in global coordinates 求解全局坐标中的一阶导数
                    dR = dxdxi\[dRdx';dRdy'];
                   
                    
                    %[~, dxdxi] = paramMap( GIFTmesh{patchIndex}, gauss_coord_x(ii), gauss_coord_y(jj), xmin, ymin, xmax, ymax);
                    % coord
                    
                    % plot(coord(1),coord(2),'+r')
                    % hold on
                    % drawnow
                    %J = det(dxdxi);
                    
                    localNd=zeros(1,nument);
                    localNd=R;
                    Nd((jj+(ii-1)*4-1)+1:(jj+(ii-1)*4),(elementCounter-1)*nument+1:elementCounter*nument)= localNd;
                    localBd=zeros(2,nument);
                    localBd(1,1:nument) = dR(1,:);
                    localBd(2,1:nument) = dR(2,:);
                    Bd(2*(jj+(ii-1)*4-1)+1:2*(jj+(ii-1)*4),(elementCounter-1)*nument+1:elementCounter*nument) = localBd;
                    localNu = zeros(2,2*nument);
                    localNu(1,1:2:2*nument-1) = R;
                    localNu(2,2:2:2*nument) = R;
                    Nu(2*(jj+(ii-1)*4-1)+1:2*(jj+(ii-1)*4),2*(elementCounter-1)*nument+1:2*elementCounter*nument) = localNu;
                    localBu = zeros(3,2*nument);
                    localBu(1,1:2:2*nument-1) = dR(1,:);
                    localBu(2,2:2:2*nument) = dR(2,:);
                    localBu(3,1:2:2*nument-1) =dR(2,:);
                    localBu(3,2:2:2*nument) = dR(1,:);
                    Bu(3*(jj+(ii-1)*4-1)+1:3*(jj+(ii-1)*4),2*(elementCounter-1)*nument+1:2*elementCounter*nument) = localBu;
                    %Bu1(3*(jj+(ii-1)*4-1)+1:3*(jj+(ii-1)*4),2*(elementCounter-1)*nument+1:2*elementCounter*nument) = localBu * scalefac * gauss_weight_x(ii).*gauss_weight_y(jj).*J;
                end
              end
        end
    end
end
Ndp = [];
Bdp = [];
for iel=1:nel
for i=1:ngauss_y*ngauss_x
Ndpx=Nd(i,16*(iel-1)+1:16*iel)'*Nd(i,16*(iel-1)+1:16*iel);
Ndp = [Ndp,Ndpx(:)];
Bdpx=Bd(2*(i-1)+1:2*i,16*(iel-1)+1:16*iel)'*Bd(2*(i-1)+1:2*i,16*(iel-1)+1:16*iel);
Bdp = [Bdp,Bdpx(:)];
end
end




