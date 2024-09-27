function [iK,jK,iKd,jKd,iFd,jFd,dedofMat,edofMat] = edofMatconnect(PHTelem)
elementCounter = 0;
for patchIndex = 1:length(PHTelem)
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            elementCounter = elementCounter + 1;
        end
    end
end
nel=elementCounter;
edofMat = zeros(nel,32);
dedofMat = zeros(nel,16);
elementCounter = 0;
for patchIndex = 1:length(PHTelem)
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            elementCounter = elementCounter+1;
            nument = size(PHTelem{patchIndex}(i).C,1);
            sctrx = PHTelem{patchIndex}(i).nodesGlobal(1:nument);
            dsctrx = reshape([2*sctrx-1; 2*sctrx],1,32);
            dedofMat(elementCounter,:) = sctrx;
            edofMat(elementCounter,:) = dsctrx;
        end
    end
end
iK = reshape(kron(edofMat,ones(32,1))',32*32*nel,1);
jK = reshape(kron(edofMat,ones(1,32))',32*32*nel,1);
iKd = reshape(kron(dedofMat,ones(16,1))',16*16*nel,1);
jKd = reshape(kron(dedofMat,ones(1,16))',16*16*nel,1);
iFd = reshape(dedofMat',16*nel,1);
jFd = ones(16*nel,1);

