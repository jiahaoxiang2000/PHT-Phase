function post_treatment(coor,numD,connect)
A=[coor numD(:,1) numD(:,2)]; %A中含有四列，前两列分别为节点的X,Y坐标，后两列分别为节点上X,Y方向的位移;A中的行数等于节点个数
matrix=[A;connect]; %在A矩阵后面添加单元信息矩阵，形成matrix矩阵
fid=fopen('C:\Users\Dell\Desktop\shuju\21.txt','at'); %按目录新建并打开1.txt文件
StringVariable='TITLE="2Dmodel" VARIABLES="X","Y","d1","d2" ZONE N=54477,E=54151,F=FEPOINT,ET=QUADRILATERAL, '; %Tecplot文件的表头
fprintf(fid, '%s',StringVariable); %把表头存储到1.txt文件中
%以下目的是把matrix按矩阵形式存储到1.txt文件中
[m,n]=size(matrix);
for i=1:m  
    for j=1:n      
        if j==n        
            fprintf(fid,'%g\n',matrix(i,j));    
        else
            fprintf(fid,'%g\t',matrix(i,j)); 
        end
    end
end
fclose(fid);%关闭1.txt文件

A=[coor numD(:,3) numD(:,4)];
matrix=[A;connect];
fid=fopen('C:\Users\Dell\Desktop\shuju\22.txt','at');
StringVariable='TITLE="2Dmodel" VARIABLES="X","Y","d3","d4" ZONE N=54477,E=54151,F=FEPOINT,ET=QUADRILATERAL, ';
fprintf(fid, '%s',StringVariable);
[m,n]=size(matrix);
for i=1:m  
    for j=1:n      
        if j==n        
            fprintf(fid,'%g\n',matrix(i,j));    
        else
            fprintf(fid,'%g\t',matrix(i,j)); 
        end
    end
end
fclose(fid);

A=[coor numD(:,5) numD(:,6)];
matrix=[A;connect];
fid=fopen('C:\Users\Dell\Desktop\shuju\23.txt','at');
StringVariable='TITLE="2Dmodel" VARIABLES="X","Y","d5","d6" ZONE N=54477,E=54151,F=FEPOINT,ET=QUADRILATERAL, ';
fprintf(fid, '%s',StringVariable);
[m,n]=size(matrix);
for i=1:m  
    for j=1:n      
        if j==n        
            fprintf(fid,'%g\n',matrix(i,j));    
        else
            fprintf(fid,'%g\t',matrix(i,j)); 
        end
    end
end
fclose(fid);

A=[coor numD(:,7) numD(:,8)];
matrix=[A;connect];
fid=fopen('C:\Users\Dell\Desktop\shuju\24.txt','at');
StringVariable='TITLE="2Dmodel" VARIABLES="X","Y","d7","d8" ZONE N=54477,E=54151,F=FEPOINT,ET=QUADRILATERAL, ';
fprintf(fid, '%s',StringVariable);
[m,n]=size(matrix);
for i=1:m  
    for j=1:n      
        if j==n        
            fprintf(fid,'%g\n',matrix(i,j));    
        else
            fprintf(fid,'%g\t',matrix(i,j)); 
        end
    end
end
fclose(fid);

A=[coor numD(:,9) numD(:,10)];
matrix=[A;connect];
fid=fopen('C:\Users\Dell\Desktop\shuju\25.txt','at');
StringVariable='TITLE="2Dmodel" VARIABLES="X","Y","d9","d10" ZONE N=54477,E=54151,F=FEPOINT,ET=QUADRILATERAL, ';
fprintf(fid, '%s',StringVariable);
[m,n]=size(matrix);
for i=1:m  
    for j=1:n      
        if j==n        
            fprintf(fid,'%g\n',matrix(i,j));    
        else
            fprintf(fid,'%g\t',matrix(i,j)); 
        end
    end
end
fclose(fid);
