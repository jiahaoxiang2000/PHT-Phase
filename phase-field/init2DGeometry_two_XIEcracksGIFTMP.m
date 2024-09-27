function GIFTmesh = init2DGeometry_two_XIEcracksGIFTMP(numPatches,L,H)
   
    GIFTmesh = cell(numPatches,1);
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %divide the patches along the x direction
    xVertices = [0,10,30,40,60];
    yVertices = [0,20,30,40,60];
    
    patchCounter = 0;
    
    for patchIndexY = 1
        for patchIndexX = 1:4
            %set the dimensions of the patch
            patchCounter = patchCounter + 1;
            patchMinX = xVertices(patchIndexX);
            patchMaxX = xVertices(patchIndexX+1);
            patchMinY = yVertices(patchIndexY);
            patchMaxY = yVertices(patchIndexY+1);
            
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
            
            GIFTmesh{patchCounter} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
            
            
        end
    end
%5th
            coefs(1:3,1,1) = [0; 20; 0];
            coefs(1:3,1,2) = [0; 30; 0];
            coefs(1:3,2,1) = [10; 20; 0];
            coefs(1:3,2,2) = [20; 30; 0];
            coefs(4,1,1) = 1;
            coefs(4,1,2) = 1;
            coefs(4,2,1) = 1;
            coefs(4,2,2) = 1;
            
            knotU = [0 0 1 1];
            knotV = [0 0 1 1];
            
            GIFTmesh{5} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
%6th
            coefs(1:3,1,1) = [10; 20; 0];
            coefs(1:3,1,2) = [20; 30; 0];
            coefs(1:3,2,1) = [30; 20; 0];
            coefs(1:3,2,2) = [30; 30; 0];
            coefs(4,1,1) = 1;
            coefs(4,1,2) = 1;
            coefs(4,2,1) = 1;
            coefs(4,2,2) = 1;
            
            knotU = [0 0 1 1];
            knotV = [0 0 1 1];
            
            GIFTmesh{6} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
            %7th
            coefs(1:3,1,1) = [30; 20; 0];
            coefs(1:3,1,2) = [30; 30; 0];
            coefs(1:3,2,1) = [40; 20; 0];
            coefs(1:3,2,2) = [40; 30; 0];
            coefs(4,1,1) = 1;
            coefs(4,1,2) = 1;
            coefs(4,2,1) = 1;
            coefs(4,2,2) = 1;
            
            knotU = [0 0 1 1];
            knotV = [0 0 1 1];
            
            GIFTmesh{7} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
            %8th
            coefs(1:3,1,1) = [40; 20; 0];
            coefs(1:3,1,2) = [40; 30; 0];
            coefs(1:3,2,1) = [60; 20; 0];
            coefs(1:3,2,2) = [60; 30; 0];
            coefs(4,1,1) = 1;
            coefs(4,1,2) = 1;
            coefs(4,2,1) = 1;
            coefs(4,2,2) = 1;
            
            knotU = [0 0 1 1];
            knotV = [0 0 1 1];
            
            GIFTmesh{8} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
             %9th
            coefs(1:3,1,1) = [0; 30; 0];
            coefs(1:3,1,2) = [0; 40; 0];
            coefs(1:3,2,1) = [20; 30; 0];
            coefs(1:3,2,2) = [20; 40; 0];
            coefs(4,1,1) = 1;
            coefs(4,1,2) = 1;
            coefs(4,2,1) = 1;
            coefs(4,2,2) = 1;
            
            knotU = [0 0 1 1];
            knotV = [0 0 1 1];
            
            GIFTmesh{9} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
            %10th
            coefs(1:3,1,1) = [20; 30; 0];
            coefs(1:3,1,2) = [20; 40; 0];
            coefs(1:3,2,1) = [30; 30; 0];
            coefs(1:3,2,2) = [30; 40; 0];
            coefs(4,1,1) = 1;
            coefs(4,1,2) = 1;
            coefs(4,2,1) = 1;
            coefs(4,2,2) = 1;
            
            knotU = [0 0 1 1];
            knotV = [0 0 1 1];
            
            GIFTmesh{10} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
            %11th
            coefs(1:3,1,1) = [30; 30; 0];
            coefs(1:3,1,2) = [30; 40; 0];
            coefs(1:3,2,1) = [40; 30; 0];
            coefs(1:3,2,2) = [50; 40; 0];
            coefs(4,1,1) = 1;
            coefs(4,1,2) = 1;
            coefs(4,2,1) = 1;
            coefs(4,2,2) = 1;
            
            knotU = [0 0 1 1];
            knotV = [0 0 1 1];
            
            GIFTmesh{11} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
             %12th
            coefs(1:3,1,1) = [40; 30; 0];
            coefs(1:3,1,2) = [50; 40; 0];
            coefs(1:3,2,1) = [60; 30; 0];
            coefs(1:3,2,2) = [60; 40; 0];
            coefs(4,1,1) = 1;
            coefs(4,1,2) = 1;
            coefs(4,2,1) = 1;
            coefs(4,2,2) = 1;
            
            knotU = [0 0 1 1];
            knotV = [0 0 1 1];
            
            GIFTmesh{12} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
            
 patchCounter=12 ;        
 xVertices = [0,20,30,50,60];
   for patchIndexY = 4
        for patchIndexX = 1:4
            %set the dimensions of the patch
            patchCounter = patchCounter + 1;
            patchMinX = xVertices(patchIndexX);
            patchMaxX = xVertices(patchIndexX+1);
            patchMinY = yVertices(patchIndexY);
            patchMaxY = yVertices(patchIndexY+1);
            
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
            
            GIFTmesh{patchCounter} = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV );
            
            
        end
    end
            
end