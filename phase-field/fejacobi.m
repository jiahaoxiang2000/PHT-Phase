function [ detJacobi ] = fejacobi( PHTelem, GIFTmesh, nel, p, q )
%uses GIFT mapping 动态映射
%supports multipatches

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


%assemble the stiffness matrix and RHS 装配刚度矩阵和RHS
detJacobi = zeros(ngauss_x*ngauss_y,nel);
%J0 = zeros(ngauss_y*ngauss_x,nel);
%Jacobi = sparse(2*sizeBasis,2*sizeBasis);
%fJacobi = sparse(sizeBasis,1);
%dJacobi = sparse(sizeBasis,sizeBasis);
%Dj = zeros(ngauss_x*ngauss_y,nel);
elementCounter = 0;
for patchIndex= 1:length(PHTelem)
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
            sctrx = PHTelem{patchIndex}(i).nodesGlobal(1:nument);
            dsctrx = reshape([2*sctrx-1; 2*sctrx],1,2*nument);
            %J1 = zeros(ngauss_y, ngauss_x);
            
            %loop over the ngauss_x x ngauss_y gauss points on each element
            %循环遍历每个元素上的ngauss_x x ngauss_y高斯点
            for jj=1:ngauss_y
                 for ii=1:ngauss_x
                    
                    %evaluate the derivatives of the mapping from parameter
                    %space to physical space 求从参数空间到物理空间的映射的导数
                    
                    [~, dxdxi] = paramMap( GIFTmesh{patchIndex}, gauss_coord_x(ii), gauss_coord_y(jj), xmin, ymin, xmax, ymax);
                    % coord
                    
                    % plot(coord(1),coord(2),'+r')
                    % hold on
                    % drawnow
                    J= det(dxdxi);
                    %J2(jj, ii) = det(dxdxi)*gauss_weight_x(ii)*gauss_weight_y(jj)*scalefac;
                    %J0(ii+(jj-1)*4,elementCounter) = det(dxdxi)*gauss_weight_x(ii)*gauss_weight_y(jj)*scalefac;
                    %localJacobi=localJacobi+scalefac * gauss_weight_x(ii).*gauss_weight_y(jj).*J;
                    %localfJacobi=localfJacobi+scalefac * gauss_weight_x(ii).*gauss_weight_y(jj).*J;
                    %localdJacobi = localdJacobi+scalefac * gauss_weight_x(ii).*gauss_weight_y(jj).*J;
                    detJacobi(ii+(jj-1)*4,elementCounter) = scalefac * gauss_weight_x(ii).*gauss_weight_y(jj).*J;
                    %Dj(ii+(jj-1)*5,elementCounter) = J
                    %TODO: implement non-zero volume force TODO:实施非零体积力
                    %localstiff = localstiff + Bu * Cmat * Bu' * scalefac * gauss_weight_x(ii).*gauss_weight_y(jj).*J;
                    
                end
            end
            
%Jacobi(dsctrx, dsctrx) = Jacobi(dsctrx, dsctrx) + localJacobi;
%fJacobi(sctrx,1) = fJacobi(sctrx, 1)+localfJacobi;
%dJacobi(sctrx, sctrx) = dJacobi(sctrx, sctrx) + localdJacobi;
        end
    end
end
