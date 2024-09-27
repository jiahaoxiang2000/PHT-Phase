function [ GIFTmesh] = genGIFTmesh( knotU, knotV, coefs, p, q, numberElementsU, numberElementsV )
%generates a GIFTmesh object from the knot vectors, control
%points/weights（用节点向量、、控制点、权重生成一个T网格）
%and polynomial degrees given
%uses the same format as NURBS toolbox（使用与NURBS相同的形式）

dim = 2; %number of physical dimensions（维数）

%tolerance for equality tests（误差）
toleq = 1e-10;

%the number of control points in the u and v directions（u和v方向的控制点数）
lenU = length(knotU)-p-1;  %number of basis functions in the u direction（u上的基函数个数）
lenV = length(knotV)-q-1;  %number of basis functions in the v direction（v上的基函数个数）
numnodes = lenU*lenV;      %number of control points（控制点个数）
coordinates = zeros(numnodes, dim+1);  %allocate one more column for the weights（为增加一维分配权重）
for j=1:lenV %for each node in the y direction(在y方向上的每个节点)
    for i=1:lenU % for each node in the x direction（在x方向上的每个节点）
        index = (j-1)*lenU + i; %the index of the coordinate array（坐标数组索引）
        coordinates(index,:) = [coefs(1,i,j)./coefs(4,i,j), coefs(2,i,j)./coefs(4,i,j), coefs(4,i,j)]; %put the (i,j) node in the coordinate array（将i,j放入相关矩阵中）
    end %for i
end %for j

%do Bezier extraction（Bezier曲线提取）

C = zeros((p+1)*(q+1),(p+1)*(q+1),numberElementsU*numberElementsV);
[C_u, ~] = bezierExtraction(knotU,p);
[C_v, ~] = bezierExtraction(knotV,q);

%calcuate the Bezier extraction operators（计算Bezier提取算子）
elementCounter = 0;%(坐标中心点)
for j=1:numberElementsV
    for i=1:numberElementsU
        elementCounter = elementCounter + 1;
        C(:,:,elementCounter) =  kron(C_v(:,:,j),C_u(:,:,i));%求张量积
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
GIFTmesh.c_net = coordinates;%最粗的网格节点矩阵
GIFTmesh.C = C;%bezier算子
GIFTmesh.elementNode = elementNode;
GIFTmesh.elementVertex =elementVertex;

end