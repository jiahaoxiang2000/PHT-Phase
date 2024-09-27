function [C, nb] = bezierExtraction(knot,p)
% Bezier extraction
% Based on Algorithm 1, from Borden - Isogeometric finite element data
% structures based on Bezier extraction(从有限元数据结构中提取Bezier曲线数据结构)


m=length(knot)-p-1;%（基函数，控制点数量）
a  = p+1;%poly多项式阶数
b  = a+1;
nb = 1;
C(:,:,1) = eye(p+1);%(单位矩阵)

while b <= m
    C(:,:,nb+1) = eye(p+1);%初始化下一个提取操作符
    i=b;
    %Count multiplicity of the knot at location b 计算b点的多重度
    while b <= m && knot(b+1) == knot(b)
        b=b+1;
    end
    
    multiplicity = b-i+1;%（重复度）
    if multiplicity < p
        numerator=knot(b)-knot(a);
        for j=p:-1:multiplicity+1
            alphas(j-multiplicity)=numerator/(knot(a+j)-knot(a));%节点插入的公式
        end
        r=p-multiplicity;%（更新r个新节点的矩阵系数）
        for j=1:r
            save = r-j+1;
            s = multiplicity + j;
            for k=p+1:-1:s+1
                alpha=alphas(k-s);
                % Form extraction operator 提取算子
                C(:,k,nb)=alpha*C(:,k,nb)+(1-alpha)*C(:,k-1,nb);
            end
            if b <= m%（更新下一个操作符的重叠系数）
                C(save:save+j,save,nb+1)=C(p-j+1:p+1,p+1,nb);
            end
        end
        %（完成当前操作符）
        %（更新下一个操作符的索引）
        nb=nb+1;
        if b <= m
            a=b;
            b=b+1;
        end
    elseif multiplicity==p %In case multiplicity of knot is already p,如果重复度为p
        if b <= m
            nb=nb+1; a=b; b=b+1;
        end
    end
end

