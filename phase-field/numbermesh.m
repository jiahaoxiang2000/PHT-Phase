function [nel] = numbermesh(PHTelem)
nel = 0;
for patchIndex = 1:length(PHTelem)
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            nel = nel + 1;
        end
    end
end
end

