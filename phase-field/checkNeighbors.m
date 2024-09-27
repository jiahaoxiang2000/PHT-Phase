function [ newBasisVert,  Tjunct, RTjunct, knotUl, knotUr, knotVd, knotVu, warningFlag ] = checkNeighbors( PHTelem, e )
%returns information about neighboring elements, to be used by the
%deCasteljau2dai routine �������Ԫ�ص���Ϣ����deCasteljau2daiʹ��

%collect information about new basis vertices and T junctions �ռ������µĻ������T������Ϣ
    
newBasisVert = []; %new basis vertices �µĻ���������
Tjunct = [];  %stores the sides that will become T junctions   ���潫��ΪT�ͽ��紦�ı�
RTjunct = []; %stores the sides where T junctions will be removed �洢T���ӽ����Ƴ��ı�
warningFlag = 0; 

%set default values for knotUl, knotUr, knotVd, knotVu ����Ĭ��ֵ
knotUl = 0;
knotUr = 1;
knotVd = 0;
knotVu = 1;

switch length(PHTelem(e).neighbor_left)
    case 0
        newBasisVert = [newBasisVert, 4];
    case 1
        Tjunct = [Tjunct, 4];
    case 2
        RTjunct = [RTjunct, 4];
        knotUl = PHTelem(PHTelem(e).neighbor_left(1)).vertex(1);        
    otherwise
        disp(['Error on element index ', num2str(e)])
        warning('Left neighbor is too refined... method not implemented')
        warningFlag = 1;
end

switch length(PHTelem(e).neighbor_right)
    case 0
        newBasisVert = [newBasisVert, 2];
    case 1
        Tjunct = [Tjunct, 2];
    case 2
        RTjunct = [RTjunct, 2];
        knotUr = PHTelem(PHTelem(e).neighbor_right(1)).vertex(3);        
    otherwise
        disp(['Error on element index ', num2str(e)])
        warning('Right neighbor is too refined... method not implemented')
        warningFlag = 1;
end        

 switch length(PHTelem(e).neighbor_down)
    case 0
        newBasisVert = [newBasisVert, 1];
    case 1
        Tjunct = [Tjunct, 1];
    case 2
        RTjunct = [RTjunct, 1];
        knotVd = PHTelem(PHTelem(e).neighbor_down(1)).vertex(2);                        
    otherwise
        disp(['Error on element index ', num2str(e)])
        warning('Down neighbor is too refined... method not implemented')
        warningFlag = 1;
 end

 switch length(PHTelem(e).neighbor_up)
    case 0
        newBasisVert = [newBasisVert, 3];
    case 1
        Tjunct = [Tjunct, 3];
    case 2
        RTjunct = [RTjunct, 3];
        knotVu = PHTelem(PHTelem(e).neighbor_up(1)).vertex(4);                
    otherwise
        disp(['Error on element index ', num2str(e)])
        warning('Up neighbor is too refined... method not implemented')
        warningFlag = 1;
 end