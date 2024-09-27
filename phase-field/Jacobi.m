function  [J0]=Jacobi( PHTelem, GIFTmesh, p, q,nel)
%calculate the actual error in the computed solution for the edge crack
%example
%supports multipatches

numPatches = length(PHTelem);


numGaussX = p+1;
numGaussY = q+1;

[gwx, gpx]=quadrature(numGaussX, 'GAUSS', 1);
[gwy, gpy]=quadrature(numGaussY, 'GAUSS', 1);

gpx=gpx';
gpy=gpy';
J0 = zeros(numGaussY*numGaussX,nel);
asj=0;
for patchIndex = 1:numPatches
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            asj=asj+1;
            xmin = PHTelem{patchIndex}(i).vertex(1);
            xmax = PHTelem{patchIndex}(i).vertex(3);
            ymin = PHTelem{patchIndex}(i).vertex(2);
            ymax = PHTelem{patchIndex}(i).vertex(4);
            
            %jacobian of the transformation from reference [-1,1]x[-1,1]
            %element to the local element in parameter
            scalefac = (xmax - xmin)*(ymax - ymin)/4;
            
            
           % J1 = zeros(numGaussY, numGaussX);
            
            
            
            nument = size(PHTelem{patchIndex}(i).C,1); %number of basis functions with support on current knotspan
            
            scrt = PHTelem{patchIndex}(i).nodesGlobal(1:nument);
            scrt_x = 2*scrt-1;
            scrt_y = 2*scrt;
            dscrtx = reshape([2*scrt-1; 2*scrt],1,2*nument);
            
            for jj=1:numGaussY
                for ii=1:numGaussX
                    
                    [ coord, dxdxi]  = paramMap( GIFTmesh{patchIndex}, gpx(ii), gpy(jj), xmin, ymin, xmax, ymax);
                    
                   
                   % stress_ex = stress_ex';
                    
                   
                    J0((jj-1)*4+ii,asj) = det(dxdxi)*gwx(ii)*gwy(jj)*scalefac;
                   
                end
            end

        end
    end
end

