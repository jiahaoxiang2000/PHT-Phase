function [GIFTmesh] = init2DGeometryGIFTMPrectangle(L,W,numPatches)
 GIFTmesh = cell(numPatches,1);
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %divide the patches along the x direction
    xVertices = linspace(0,L,numPatches+1);
    
    for patchIndex = 1:numPatches
        
        %set the dimensions of the patch
        patchMinX = xVertices(patchIndex);
        patchMaxX = xVertices(patchIndex+1);
        patchMinY = 0;
        patchMaxY = W;
        
        %initialize geometry on coarsest mesh
        coefs(1:3,1,1) = [patchMinX; patchMinY; 0];
        coefs(1:3,1,2) = [patchMinX; patchMaxY; 0];
        coefs(1:3,2,1) = [patchMaxX; patchMinY; 0];
        coefs(1:3,2,2) = [patchMaxX; patchMaxY; 0];
        coefs(4,1,1) = 1;
        coefs(4,1,2) = 1;
        coefs(4,2,1) = 1;
        coefs(4,2,2) = 1;
        
        knotU = [0 0 1 1];
        knotV = [0 0 1 1];
        
        GIFTmesh{patchIndex} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
        
   
end

