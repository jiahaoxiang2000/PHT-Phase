function [ GIFTmesh] = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV )
%generates a GIFTmesh object from the knot vectors, control
%points/weights���ýڵ������������Ƶ㡢Ȩ������һ��T����
%and polynomial degrees given
%uses the same format as NURBS toolbox��ʹ����NURBS��ͬ����ʽ��

dim = 2; %number of physical dimensions��ά����

%tolerance for equality tests����
toleq = 1e-10;

%the number of control points in the u and v directions��u��v����Ŀ��Ƶ�����
lenU = length(knotU)-p-1;  %number of basis functions in the u direction��u�ϵĻ�����������
lenV = length(knotV)-q-1;  %number of basis functions in the v direction��v�ϵĻ�����������
numnodes = lenU*lenV;      %number of control points�����Ƶ������
coordinates = zeros(numnodes, dim+1);  %allocate one more column for the weights��Ϊ����һά����Ȩ�أ�
for j=1:lenV %for each node in the y direction(��y�����ϵ�ÿ���ڵ�)
    for i=1:lenU % for each node in the x direction����x�����ϵ�ÿ���ڵ㣩
        index = (j-1)*lenU + i; %the index of the coordinate array����������������
        coordinates(index,:) = [coefs(1,i,j)./coefs(4,i,j), coefs(2,i,j)./coefs(4,i,j), coefs(4,i,j)]; %put the (i,j) node in the coordinate array����i,j������ؾ����У�
    end %for i
end %for j

%do Bezier extraction��Bezier������ȡ��

C = zeros((p+1)*(q+1),(p+1)*(q+1),numberElementsU*numberElementsV);
[C_u, ~] = bezierExtraction(knotU,p);
[C_v, ~] = bezierExtraction(knotV,q);

%calcuate the Bezier extraction operators������Bezier��ȡ���ӣ�
elementCounter = 0;%(�������ĵ�)
for j=1:numberElementsV
    for i=1:numberElementsU
        elementCounter = elementCounter + 1;
        C(:,:,elementCounter) =  kron(C_v(:,:,j),C_u(:,:,i));%��������
    end
end

%make elementVertex and elementNode arrays
elementCounter = 0;
elementVertex = zeros(numberElementsU*numberElementsV,4);
elementNode = zeros(numberElementsU*numberElementsV, (p+1)*(q+1));

for j=1:length(knotV)-1
    for i=1:length(knotU)-1
        if (abs(knotU(i+1)-knotU(i))>toleq) && (abs(knotV(j+1)-knotV(j))>toleq)
            elementCounter = elementCounter + 1;
            elementVertex(elementCounter, :) = [knotU(i), knotV(j), knotU(i+1), knotV(j+1)];
            tcount = 0;
            currow = zeros(1, (p+1)*(q+1));
            %now we add the nodes from i-p...i in the u direction and
            %j-q...j in the v direction
            for t2=j-q:j
                for t1 = i-p:i
                    tcount = tcount + 1;
                    currow(tcount) = t1+(t2-1)*lenU;
                end
            end
            elementNode(elementCounter,:)=currow;
        end
    end
end

GIFTmesh.numberElements = numberElementsU*numberElementsV;
GIFTmesh.numberElementsU = numberElementsU;
GIFTmesh.numberElementsV = numberElementsV;

GIFTmesh.p = p;
GIFTmesh.q = q;
GIFTmesh.c_net = coordinates;%��ֵ�����ڵ����
GIFTmesh.C = C;%bezier����
GIFTmesh.elementNode = elementNode;
GIFTmesh.elementVertex =elementVertex;

end