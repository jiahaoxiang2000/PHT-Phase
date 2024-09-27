function  plotSolPHTElasticMPVM( d, PHTelem, GIFTmesh, p, q )
%supports multipatches
%plots the von Misses stresses

numPts = 2;

%if nargin<7
    %set a default magnification factor for the displacement
    %factor=10;
%end

%plots the deformed shape + stresses
fudge =0;
xi = linspace(-1+fudge,1-fudge,numPts);
eta = linspace(-1+fudge,1-fudge,numPts);


%calculate the number of actual elements (i.e., non-refined, without children)
numElem = 0;
numPatches = length(PHTelem);
for patchIndex=1:numPatches
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            numElem = numElem+1;
        end
    end
end


%Displaying the displacements
element4 = zeros(numElem, 4);
physcoord = zeros(4*numElem, 2);
dcoord = zeros(4*numElem, 1);

%define the 2D Bernstein polynomials
[B_u, dB_u] = bernstein_basis(xi,p);
[B_v, dB_v] = bernstein_basis(eta,q);

Buv = zeros(numPts, numPts, numPts*numPts);
dBdu = zeros(numPts, numPts, numPts*numPts);
dBdv = zeros(numPts, numPts, numPts*numPts);

basisCounter = 0;
for j=1:q+1
    for i=1:p+1
        basisCounter = basisCounter + 1;
        Buv(:,:,basisCounter) = B_u(:,i)*B_v(:,j)';
        dBdu(:,:,basisCounter) = dB_u(:,i)*B_v(:,j)';
        dBdv(:,:,basisCounter) = B_u(:,i)*dB_v(:,j)';
    end
end

elementCounter = 0;
for indexPatch = 1:numPatches
    for i=1:length(PHTelem{indexPatch})
        if isempty(PHTelem{indexPatch}(i).children)
            
            elementCounter = elementCounter + 1;
            element4(elementCounter, :) = (elementCounter-1)*4+1:(elementCounter-1)*4+4;
            xmin = PHTelem{indexPatch}(i).vertex(1);
            xmax = PHTelem{indexPatch}(i).vertex(3);
            ymin = PHTelem{indexPatch}(i).vertex(2);
            ymax = PHTelem{indexPatch}(i).vertex(4);
            coord = cell(numPts,numPts);
            
            nument = size(PHTelem{indexPatch}(i).C,1); %number of basis functions with support on current knotspan
            
            %initialize matrices to store the displacement, strain and stress
            %values at each plot point
            dmat = zeros(numPts,numPts);
            %dispmatx = zeros(numPts,numPts);
          
            
            scrt = PHTelem{indexPatch}(i).nodesGlobal(1:nument);
            
            for jj=1:numPts
                for ii=1:numPts
                    
                    %compute the mapping from parameter space to physical space
                    [ coord_pt, ~]  = paramMap( GIFTmesh{indexPatch}, xi(ii), eta(jj), xmin, ymin, xmax, ymax);
                    
                    %evaluate the basis functions
                    
                    cR = PHTelem{indexPatch}(i).C * squeeze(Buv(ii,jj,:));

                    
                    coord{jj,ii} = [coord_pt(1), coord_pt(2)];
                    
                    % Solve for first derivatives in global coordinates
                    
                    %calculate displacement values
                    dmat(jj,ii) = dmat(jj,ii) + cR'*d(scrt);
                
                    
                    %calculate the stress values
                end
            end
            
            physcoord((elementCounter-1)*4+1, :) = coord{1,1};
            physcoord((elementCounter-1)*4+2, :) = coord{1,end};
            physcoord((elementCounter-1)*4+3, :) = coord{end,end};
            physcoord((elementCounter-1)*4+4, :) = coord{end,1};
            
            dcoord((elementCounter-1)*4+1, :) = dmat(1,1);
            dcoord((elementCounter-1)*4+2, :) = dmat(1,2);
            dcoord((elementCounter-1)*4+3, :) = dmat(2,2);
            dcoord((elementCounter-1)*4+4, :) = dmat(2,1);
            
            hold on
        end
    end
end


figure
trisurf(element4,physcoord(:,1), physcoord(:,2), zeros(size(physcoord,1),1),dcoord, 'EdgeColor','none','facecolor','interp')
view(0,90)
title(' phase field d ');
colorbar('vert')
drawnow
