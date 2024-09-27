function [C, nb] = bezierExtraction(knot,p)
% Bezier extraction
% Based on Algorithm 1, from Borden - Isogeometric finite element data
% structures based on Bezier extraction(������Ԫ���ݽṹ����ȡBezier�������ݽṹ)


m=length(knot)-p-1;%�������������Ƶ�������
a  = p+1;%poly����ʽ����
b  = a+1;
nb = 1;
C(:,:,1) = eye(p+1);%(��λ����)

while b <= m
    C(:,:,nb+1) = eye(p+1);%��ʼ����һ����ȡ������
    i=b;
    %Count multiplicity of the knot at location b ����b��Ķ��ض�
    while b <= m && knot(b+1) == knot(b)
        b=b+1;
    end
    
    multiplicity = b-i+1;%���ظ��ȣ�
    if multiplicity < p
        numerator=knot(b)-knot(a);
        for j=p:-1:multiplicity+1
            alphas(j-multiplicity)=numerator/(knot(a+j)-knot(a));%�ڵ����Ĺ�ʽ
        end
        r=p-multiplicity;%������r���½ڵ�ľ���ϵ����
        for j=1:r
            save = r-j+1;
            s = multiplicity + j;
            for k=p+1:-1:s+1
                alpha=alphas(k-s);
                % Form extraction operator ��ȡ����
                C(:,k,nb)=alpha*C(:,k,nb)+(1-alpha)*C(:,k-1,nb);
            end
            if b <= m%��������һ�����������ص�ϵ����
                C(save:save+j,save,nb+1)=C(p-j+1:p+1,p+1,nb);
            end
        end
        %����ɵ�ǰ��������
        %��������һ����������������
        nb=nb+1;
        if b <= m
            a=b;
            b=b+1;
        end
    elseif multiplicity==p %In case multiplicity of knot is already p,����ظ���Ϊp
        if b <= m
            nb=nb+1; a=b; b=b+1;
        end
    end
end

