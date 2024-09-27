function [MP_quadRef,level] = Phase_mesh(numPatches,PHTelem,quadList,d,fai)
d_mesh = cell(1,numPatches);
MP_nodes = find(d>=fai);
MP_quadRef = cell(1, numPatches);
numQuads = 0;
for Indexpatch = 1:numPatches
    d_mesh{Indexpatch}=zeros(1,length(PHTelem{Indexpatch}));
    MP_quadRef{Indexpatch} = zeros(1, size(quadList{Indexpatch},1));
    for j=1:size(quadList{Indexpatch},1)
        for e=1:4
            d_glob = quadList{Indexpatch}(j,e);
            nument = size(PHTelem{Indexpatch}(d_glob).C,1);
            scrtx = PHTelem{Indexpatch}(d_glob).nodesGlobal(1:nument);
            if ~isempty(intersect(MP_nodes,scrtx))
                d_mesh{Indexpatch}(d_glob) = 1;
            end
        end
    end
    numQuads = numQuads+size(quadList{Indexpatch},1);
end
quadCounter = 0;
quad_d = zeros(1,numQuads);
quad_level = zeros(1,numQuads);
quadIndex = zeros(numQuads,2);
for indexPatch = 1:numPatches
    for j = 1:size(quadList{indexPatch,1})
        quadCounter = quadCounter +1;
        quad_d(quadCounter) = sum(d_mesh{indexPatch}(quadList{indexPatch}(j,:)));
        quad_level(quadCounter) = PHTelem{indexPatch}(quadList{indexPatch}(j,1)).level;
        quadIndex(quadCounter,:) = [indexPatch,j];
    end
end
level =max(quad_level)+1;
if level>4
    level=4;
end
[quad_d,index]=sort(quad_d,'descend');
 quadIndex = quadIndex(index,:);
 quad_level = quad_level(:,index);
for i=1:numQuads
    if quad_d(i)>=1  &&  quad_level(i)<4
       MP_quadRef{quadIndex(i,1)}(quadIndex(i,2))=1;
    else
       MP_quadRef{quadIndex(i,1)}(quadIndex(i,2))=0;
    end
end     
