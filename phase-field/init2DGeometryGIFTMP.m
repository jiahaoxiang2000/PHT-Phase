function GIFTmesh = init2DGeometryGIFTMP(numPatches,L,H)
  GIFTmesh = cell(numPatches,1);
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %patch 1
    
    %set the dimensions of the patch
    patchMinX = 0;
    patchMaxX = 0.5;
    patchMinY = 0;
    patchMaxY = 0.5;
    
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
    
    GIFTmesh{1} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %patch 2
    
    %set the dimensions of the patch
    patchMinX = 0.5;
    patchMaxX = 1;
    patchMinY = 0;
    patchMaxY = 0.5;
    
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
    
    GIFTmesh{2} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %patch 3
    
    %set the dimensions of the patch
    patchMinX = 0.5;
    patchMaxX = 1;
    patchMinY = 0.5;
    patchMaxY = 1;
    
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
    
    GIFTmesh{3} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
    
    %patch 4
    
    %set the dimensions of the patch
    patchMinX = 0;
    patchMaxX = 0.5;
    patchMinY = 0.5;
    patchMaxY = 1;
    
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
    
    GIFTmesh{4} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
   
   
   
end


