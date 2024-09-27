function GIFTmesh = init2DGeometry_two_cracksGIFTMP(numPatches,L,H)
   
    GIFTmesh = cell(numPatches,1);
    numberElementsU = 1;
    numberElementsV = 1;
    p = 1;
    q = 1;
    
    %divide the patches along the x direction
    xVertices = linspace(0,L,5);
    yVertices = linspace(0,H,5);
    
    patchCounter = 0;
    
    for patchIndexY = 1:4
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
