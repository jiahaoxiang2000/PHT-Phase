function [PHTmesh] =init2DGeometryPHTMP(L,H,numPatches)
    PHTmesh = cell(numPatches,1);
    numberElementsU = 1;
    numberElementsV = 1;
    p=1;
    q=1;
%divide the patches along the x direction
xVertices = linspace(0,L,numPatches-1);
yVertices = linspace(0,H,numPatches-1);

patchCounter=0;

for patchIndex = 1:(numPatches-2)
   for patchIndey = 1:(numPatches-2)
    patchCounter = patchCounter + 1
    %set the dimensions of the patch
    patchMinX = xVertices(patchIndex);
    patchMaxX = xVertices(patchIndex+1);
    patchMinY = yVertices(patchIndey);
    patchMaxY = yVertices(patchIndey+1);
    
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
   end
end

tempPHTmesh = PHTmesh{4};
    PHTmesh{4} = PHTmesh{3};
    PHTmesh{3} = tempPHTmesh;