function [ ustiff] = assembleKu_phase( PHTelem, GIFTmesh, sizeBasis, p, q,lambda,mu,u,d,k)
%assembles the stiffness matrix and rhs (Galerkin method)
%uses GIFT mapping
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

%1D bernstein polynomials evaluated at the Gauss points on the master element
[B_u, dB_u] = bernstein_basis(gauss_coord_x,p);
[B_v, dB_v] = bernstein_basis(gauss_coord_y,q);

dBdu = zeros(ngauss_x, ngauss_y, (p+1)*(q+1));
dBdv = zeros(ngauss_x, ngauss_y, (p+1)*(q+1));
Buv = zeros(ngauss_x, ngauss_y, (p+1)*(q+1));


%the derivatives of the 2D Bernstein polynomials at Gauss points on the
%master element
basisCounter = 0;
for j=1:q+1
    for i=1:p+1
        basisCounter = basisCounter + 1;
        Buv(:,:,basisCounter) = B_u(:,i)*B_v(:,j)';
        dBdu(:,:,basisCounter) = dB_u(:,i)*B_v(:,j)';
        dBdv(:,:,basisCounter) = B_u(:,i)*dB_v(:,j)';
    end
end

%allocate memory for the triplet arrays
indexCounter = 0;
for patchIndex=1:length(PHTelem)
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            nument = size(PHTelem{patchIndex}(i).C,1);
            indexCounter = indexCounter + 4*nument^2;
        end
    end
end
%indexCounter
II = zeros(1,indexCounter);
JJ = zeros(1,indexCounter);
S = zeros(1, indexCounter);

indexCounter = 0;

%initialize LHS stiffness matrix and RHS vector
dim = 2; %the dimension of physical space
ustiff = sparse(dim*sizeBasis,dim*sizeBasis);
%urhs = zeros(dim*sizeBasis,1);

%assemble the stiffness matrix and RHS
elementCounter = 0;
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
            sctrx = PHTelem{patchIndex}(i).nodesGlobal(1:nument);
            dsctrx = reshape([2*sctrx-1; 2*sctrx],1,2*nument);
            
            localustiff = zeros(2*nument, 2*nument); %local stiffness
            
            %loop over the ngauss_x x ngauss_y gauss points on each element
            for jj=1:ngauss_y
                for ii=1:ngauss_x
                    
                    %evaluate the derivatives of the mapping from parameter
                    %space to physical space
                    
                    [~, dxdxi] = paramMap( GIFTmesh{patchIndex}, gauss_coord_x(ii), gauss_coord_y(jj), xmin, ymin, xmax, ymax);
                    % coord
                    
                    % plot(coord(1),coord(2),'+r')
                    % hold on
                    % drawnow
                    J = det(dxdxi);
                    
                    R = (PHTelem{patchIndex}(i).C)*squeeze(Buv(ii,jj,:));
                    dRdx = (PHTelem{patchIndex}(i).C)*squeeze(dBdu(ii,jj,:));
                    dRdy = (PHTelem{patchIndex}(i).C)*squeeze(dBdv(ii,jj,:));
                    
                    %multiply by the jacobian of the transformation from reference
                    %space to the parameter space
                    dRdx = dRdx*2/(xmax-xmin);
                    dRdy = dRdy*2/(ymax-ymin);
                    
                    % Solve for first derivatives in global coordinates
                    dR = dxdxi\[dRdx';dRdy'];
                    Nd = zeros(nument,1);
                    Nd = R;
                    Nu = zeros(2*nument,2);
                    Nu(1:2:2*nument,1) = R;
                    Nu(2:2:2*nument,2) = R;
                    Bu = zeros(2*nument,3);
                    Bu(1:2:2*nument-1,1) = dR(1,:);
                    Bu(2:2:2*nument,2) = dR(2,:);
                    Bu(1:2:2*nument-1,3) = dR(2,:);
                    Bu(2:2:2*nument,3) = dR(1,:);
                    epsilon_gauss = zeros(3,1);
                    d_gauss = zeros(1,1);
                    U = u(dsctrx);
                    DD = d(sctrx);
                    epsilon_gauss = Bu'*U;
                    d_gauss = Nd'*DD;
                    [R_plus,R_minus] = operator_R(epsilon_gauss);%R+ R-
                    P_plus = zeros(9,1); P_minus = P_plus;%  P+ P-
                    [P_plus,P_minus] = operator_P2(epsilon_gauss); %  P+ P-
                    I = [1;1;0]; Ip = reshape(I*I',9,1);
                    D = kron(((1-d_gauss).^2+k),ones(9,1)).*(kron(lambda*R_plus,Ip)+2*kron(mu,repmat([1;1;0.5],3,1)).*P_plus)...
                       + kron(lambda*R_minus,Ip) + 2*kron(mu,repmat([1;1;0.5],3,1)).*P_minus;
                    %TODO: implement non-zero volume force
                    localustiff = localustiff + (Bu*reshape(D,3,3)*Bu')*scalefac * gauss_weight_x(ii).*gauss_weight_y(jj).*J;
                end
            end
            %drhs(scrtx, scrtx) = drhs(scrtx, scrtx) + localdrhs;
            II(indexCounter+1:indexCounter+4*nument^2) = repmat(dsctrx,1,2*nument);
            JJ(indexCounter+1:indexCounter+4*nument^2) = reshape(repmat(dsctrx,2*nument,1),1,4*nument^2);
            S(indexCounter+1:indexCounter+4*nument^2) = reshape(localustiff,1,4*nument^2);
            indexCounter = indexCounter + 4*nument^2;
        end
    end
end

ustiff = sparse(II,JJ,S,2*sizeBasis,2*sizeBasis);
disp(['The mesh has ', num2str(elementCounter), ' active elements.'])