function [Mp_mesh] = Mesh_Phase(PHTelem,faiT,d)
elementCounter = 0;
Mp_mesh=zeros(1-length(PHTelem),2);
P_nodes=find(d>=faiT)';
for patchIndex = 1:length(PHTelem)
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            nument = size(PHTelem{patchIndex}(i).C,1);
            sctrx = PHTelem{patchIndex}(i).nodesGlobal(1:nument);
            nodes_P = intersect(P_nodes,sctrx);
         if ~isempty(nodes_P)
            elementCounter = elementCounter + 1;
            Mp_mesh(elementCounter,1)=patchIndex;
            Mp_mesh(elementCounter,2)=i;
        end
    end
end
end
end


