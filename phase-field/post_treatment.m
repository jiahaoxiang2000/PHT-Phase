function post_treatment(coor,numD,connect)
A=[coor numD(:,1) numD(:,2)]; %A�к������У�ǰ���зֱ�Ϊ�ڵ��X,Y���꣬�����зֱ�Ϊ�ڵ���X,Y�����λ��;A�е��������ڽڵ����
matrix=[A;connect]; %��A���������ӵ�Ԫ��Ϣ�����γ�matrix����
fid=fopen('C:\Users\Dell\Desktop\shuju\21.txt','at'); %��Ŀ¼�½�����1.txt�ļ�
StringVariable='TITLE="2Dmodel" VARIABLES="X","Y","d1","d2" ZONE N=54477,E=54151,F=FEPOINT,ET=QUADRILATERAL, '; %Tecplot�ļ��ı�ͷ
fprintf(fid, '%s',StringVariable); %�ѱ�ͷ�洢��1.txt�ļ���
%����Ŀ���ǰ�matrix��������ʽ�洢��1.txt�ļ���
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
fclose(fid);%�ر�1.txt�ļ�

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
