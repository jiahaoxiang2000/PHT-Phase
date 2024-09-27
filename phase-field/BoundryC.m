function [loaddofs,fixeddofs,freedofs] =BoundryC(PHTelem, p, q,sizeBasis)
%detect which nodes have support on the left and right boundary
bcdof_left = [];
bcdof_right = [];
bcdof_up = [];
bcdof_down = [];
bcdof_mid = [];
%for each neumann edge store the element index and orientation
%orientation: 1-down, 2-right, 3-up, 4-left
%neumann_left = [];
%neumann_right = [];
%neumann_up = [];
%neumann_down =[];
%neumann_mid = [];

%define side node indices
down_nodes = 1:p+1;
right_nodes = (p+1):(p+1):(p+1)*(q+1);
up_nodes = 1+(p+1)*q:(p+1)*(q+1);
left_nodes = 1:(p+1):(1+(p+1)*q);



for patchIndex=1:length(PHTelem)
    for i=1:length(PHTelem{patchIndex})
        if isempty(PHTelem{patchIndex}(i).children)
            
           if (patchIndex==1) && (PHTelem{patchIndex}(i).vertex(1)==0) 
               bcdof_left=[bcdof_left,PHTelem{patchIndex}(i).nodesGlobal(left_nodes)];
           end
           if (patchIndex==4) && (PHTelem{patchIndex}(i).vertex(1)==0) 
               bcdof_left=[bcdof_left,PHTelem{patchIndex}(i).nodesGlobal(left_nodes)];
           end
          
          
           
           if  (patchIndex==1) && (PHTelem{patchIndex}(i).vertex(2)==0)
               bcdof_down=[bcdof_down,PHTelem{patchIndex}(i).nodesGlobal(down_nodes)];
           end
           if (patchIndex==2) && (PHTelem{patchIndex}(i).vertex(2)==0)
               bcdof_down=[bcdof_down,PHTelem{patchIndex}(i).nodesGlobal(down_nodes)];
           end

          
           
           if (patchIndex==3) && (PHTelem{patchIndex}(i).vertex(3)==1)
               bcdof_right=[bcdof_right,PHTelem{patchIndex}(i).nodesGlobal(right_nodes)];
           end
           if (patchIndex==2) && (PHTelem{patchIndex}(i).vertex(3)==1) 
               bcdof_right=[bcdof_right,PHTelem{patchIndex}(i).nodesGlobal(right_nodes)];
           end
          
           
           if (patchIndex==3) && (PHTelem{patchIndex}(i).vertex(4)==1) 
               bcdof_up=[bcdof_up,PHTelem{patchIndex}(i).nodesGlobal(up_nodes)];
           end
           if (patchIndex==4) && (PHTelem{patchIndex}(i).vertex(4)==1)
               bcdof_up=[bcdof_up,PHTelem{patchIndex}(i).nodesGlobal(up_nodes)];
           end

          
           
           %if (patchIndex==1) && (PHTelem{patchIndex}(i).vertex(4)==1)  && (PHTelem{patchIndex}(i).vertex(3)==1) 
              % bcdof_mid=[bcdof_mid,PHTelem{patchIndex}(i).nodesGlobal(up_nodes)];
               % bcdof_mid=[bcdof_mid,PHTelem{patchIndex}(i).nodesGlobal(right_nodes)];
          % end
           %bcdof_mid1=[16];
           
           % if (patchIndex==4) && (PHTelem{patchIndex}(i).vertex(2)==0)  && (PHTelem{patchIndex}(i).vertex(3)==1) 
               %bcdof_mid=[bcdof_mid,PHTelem{patchIndex}(i).nodesGlobal(down_nodes)];
                %bcdof_mid=[bcdof_mid,PHTelem{patchIndex}(i).nodesGlobal(right_nodes)];
            %end
          %if (patchIndex==1) && (PHTelem{patchIndex}(i).vertex(4)==1)
               %bcdof_mid=[bcdof_mid,PHTelem{patchIndex}(i).nodesGlobal(up_nodes)];
          %end
           %if (patchIndex==4) && (PHTelem{patchIndex}(i).vertex(2)==0)  
               %bcdof_mid=[bcdof_mid,PHTelem{patchIndex}(i).nodesGlobal(down_nodes)];
              % bcdof_mid=[bcdof_mid,PHTelem{patchIndex}(i).nodesGlobal(left_nodes)];
           %end
        end
    end
end

%remove duplicated entries
bcdof_right = unique(bcdof_right);
bcdof_down = unique(bcdof_down);
bcdof_up = unique(bcdof_up);
bcdof_left = unique(bcdof_left);
%bcdof_mid = unique(bcdof_mid);


%d=zeros(sizeBasis,1);
%dalldofs=(1:sizeBasis);
%dfixed=bcdof_mid;
%dfree=setdiff(dalldofs,dfixed);
%d(dfixed)=1;



alldofs = (1:2*sizeBasis);
loaddofs = 2*bcdof_up-1;
fixeddofs1=2*bcdof_down-1;
fixeddofs2=2*bcdof_down;
fixeddofs3=2*bcdof_right;
fixeddofs4=2*bcdof_left;
fixeddofs12=union(fixeddofs1,fixeddofs2);%f1和f2的并集
fixeddofs34=union(fixeddofs3,fixeddofs4);
fixeddofs1234=union(fixeddofs12,fixeddofs34);
%fixeddofs5=2*bcdof_mid;
%fixeddofs12345=union(fixeddofs1234,fixeddofs5);
fixeddofs8=2*bcdof_up;
fixeddofs=union(fixeddofs8,fixeddofs1234);

loaddofs = loaddofs(:);
fixeddofs = fixeddofs(:);
freedofs = setdiff(alldofs,union(fixeddofs,loaddofs));
end

