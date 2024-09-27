function [ PHTelem, dimBasis ] = initPHTmesh( p,q )
%initialize the PHT geometry on coarse mesh 初始化粗网格上的PHT几何

%start with num_elements_u, num_elements_v uniformly distributed knot spans
%in each parametric direction 在每个参数方向均匀分布节点
num_elements_u = 1;
num_elements_v = 1;
dimBasis = (p+1)*(q+1); %coarsest 1x1 mesh has dimension (p+1)*(q+1) 最粗的网格尺寸

knotU = [zeros(1,p+1), (1:num_elements_u)/num_elements_u, ones(1,p)];
knotV = [zeros(1,q+1), (1:num_elements_v)/num_elements_v, ones(1,q)];

%repeat the interior knots p-1 times 重复内部节点p-1次
rep_knotU = linspace(0,1,num_elements_u+1);%=[0,1]
rep_knotU = rep_knotU(2:end-1);%空矩阵预留储存位置
rep_knotU = repmat(rep_knotU,1,p-2);

rep_knotV = linspace(0,1,num_elements_v+1);
rep_knotV = rep_knotV(2:end-1);
rep_knotV = repmat(rep_knotV,1,q-2);

knotU = sort([knotU, rep_knotU]);%节点向量是一个单调不减的整数数组
knotV = sort([knotV, rep_knotV]);


[C_u, ~] = bezierExtraction(knotU,p);
[C_v, ~] = bezierExtraction(knotV,q);

PHTelem = struct;

%define the 1st element
PHTelem.parent = [];
PHTelem.children = [];
PHTelem.C = kron(C_v(:,:,1),C_u(:,:,1)); %bezier算子
PHTelem.vertex = [0, 0, 1, 1];
PHTelem.nodes = 1:(p+1)*(q+1);
PHTelem.neighbor_left = [];%左边界
PHTelem.neighbor_right = [];%右边界
PHTelem.neighbor_down = [];%下边界
PHTelem.neighbor_up = [];%上边界
PHTelem.level = 0;%水平

[ PHTelem, dimBasis ] = crossInsert( PHTelem, 1, dimBasis, p, q );
%[ PHTelem, dimBasis ] = crossInsert( PHTelem, 5, dimBasis, p, q );
