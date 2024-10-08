function [P,Q] = operator_P(A)
[V,D] = eig([A(1),A(3)/2; A(3)/2,A(2)]);
P1=V(:,1)*V(:,1)'; P2=V(:,2)*V(:,2)'; 

P11=0.5*(1+abs(D(1,1))/D(1,1))*P1(1,1)*P1(1,1)+0.5*(1+abs(D(2,2))/D(2,2))*P2(1,1)*P2(1,1)+...
    0.5*(D(1,1)+abs(D(1,1)))*(2*P1(1,1)*P2(1,1)/(D(1,1)-D(2,2)))+0.5*(D(2,2)+abs(D(2,2)))*(2*P1(1,1)*P2(1,1)/(D(2,2)-D(1,1)));
P12=0.5*(1+abs(D(1,1))/D(1,1))*P1(1,1)*P1(2,2)+0.5*(1+abs(D(2,2))/D(2,2))*P2(1,1)*P2(2,2)+...
    0.5*(D(1,1)+abs(D(1,1)))*(2*P1(1,2)*P2(1,2)/(D(1,1)-D(2,2)))+0.5*(D(2,2)+abs(D(2,2)))*(2*P1(1,2)*P2(1,2)/(D(2,2)-D(1,1)));  
P13=0.5*(1+abs(D(1,1))/D(1,1))*P1(1,1)*P1(1,2)+0.5*(1+abs(D(2,2))/D(2,2))*P2(1,1)*P2(1,2)+...
    0.5*(D(1,1)+abs(D(1,1)))*(P1(1,1)*P2(1,2)+P1(1,2)*P2(1,1))/(D(1,1)-D(2,2))+0.5*(D(2,2)+abs(D(2,2)))*(P1(1,1)*P2(1,2)+P1(1,2)*P2(1,1))/(D(2,2)-D(1,1));
P22=0.5*(1+abs(D(1,1))/D(1,1))*P1(2,2)*P1(2,2)+0.5*(1+abs(D(2,2))/D(2,2))*P2(2,2)*P2(2,2)+...
    0.5*(D(1,1)+abs(D(1,1)))*(2*P1(2,2)*P2(2,2)/(D(1,1)-D(2,2)))+0.5*(D(2,2)+abs(D(2,2)))*(2*P1(2,2)*P2(2,2)/(D(2,2)-D(1,1)));
P23=0.5*(1+abs(D(1,1))/D(1,1))*P1(2,2)*P1(1,2)+0.5*(1+abs(D(2,2))/D(2,2))*P2(2,2)*P2(1,2)+...
    0.5*(D(1,1)+abs(D(1,1)))*(P1(2,1)*P2(2,2)+P1(2,2)*P2(2,1))/(D(1,1)-D(2,2))+0.5*(D(2,2)+abs(D(2,2)))*(P1(2,1)*P2(2,2)+P1(2,2)*P2(2,1))/(D(2,2)-D(1,1));
P33=0.5*(1+abs(D(1,1))/D(1,1))*P1(1,2)*P1(1,2)+0.5*(1+abs(D(2,2))/D(2,2))*P2(1,2)*P2(1,2)+...
    0.5*(D(1,1)+abs(D(1,1)))*0.5*(P1(1,1)*P2(2,2)+P1(1,2)*P2(2,1)+P1(2,2)*P2(1,1)+P1(2,1)*P2(1,2))/(D(1,1)-D(2,2))+0.5*(D(2,2)+abs(D(2,2)))*0.5*(P1(1,1)*P2(2,2)+P1(1,2)*P2(2,1)+P1(2,2)*P2(1,1)+P1(2,1)*P2(1,2))/(D(2,2)-D(1,1));
P=[P11 P12 P13;P12 P22 P23;P13 P23 P33];
%P=reshape(P,9,1);

Q11=0.5*(1-abs(D(1,1))/D(1,1))*P1(1,1)*P1(1,1)+0.5*(1-abs(D(2,2))/D(2,2))*P2(1,1)*P2(1,1)+...
    0.5*(D(1,1)-abs(D(1,1)))*(2*P1(1,1)*P2(1,1)/(D(1,1)-D(2,2)))+0.5*(D(2,2)-abs(D(2,2)))*(2*P1(1,1)*P2(1,1)/(D(2,2)-D(1,1)));
Q12=0.5*(1-abs(D(1,1))/D(1,1))*P1(1,1)*P1(2,2)+0.5*(1-abs(D(2,2))/D(2,2))*P2(1,1)*P2(2,2)+...
    0.5*(D(1,1)-abs(D(1,1)))*(2*P1(1,2)*P2(1,2)/(D(1,1)-D(2,2)))+0.5*(D(2,2)-abs(D(2,2)))*(2*P1(1,2)*P2(1,2)/(D(2,2)-D(1,1)));  
Q13=0.5*(1-abs(D(1,1))/D(1,1))*P1(1,1)*P1(1,2)+0.5*(1-abs(D(2,2))/D(2,2))*P2(1,1)*P2(1,2)+...
    0.5*(D(1,1)-abs(D(1,1)))*(P1(1,1)*P2(1,2)+P1(1,2)*P2(1,1))/(D(1,1)-D(2,2))+0.5*(D(2,2)-abs(D(2,2)))*(P1(1,1)*P2(1,2)+P1(1,2)*P2(1,1))/(D(2,2)-D(1,1));
Q22=0.5*(1-abs(D(1,1))/D(1,1))*P1(2,2)*P1(2,2)+0.5*(1-abs(D(2,2))/D(2,2))*P2(2,2)*P2(2,2)+...
    0.5*(D(1,1)-abs(D(1,1)))*(2*P1(2,2)*P2(2,2)/(D(1,1)-D(2,2)))+0.5*(D(2,2)-abs(D(2,2)))*(2*P1(2,2)*P2(2,2)/(D(2,2)-D(1,1)));
Q23=0.5*(1-abs(D(1,1))/D(1,1))*P1(2,2)*P1(1,2)+0.5*(1-abs(D(2,2))/D(2,2))*P2(2,2)*P2(1,2)+...
    0.5*(D(1,1)-abs(D(1,1)))*(P1(2,1)*P2(2,2)+P1(2,2)*P2(2,1))/(D(1,1)-D(2,2))+0.5*(D(2,2)-abs(D(2,2)))*(P1(2,1)*P2(2,2)+P1(2,2)*P2(2,1))/(D(2,2)-D(1,1));
Q33=0.5*(1-abs(D(1,1))/D(1,1))*P1(1,2)*P1(1,2)+0.5*(1-abs(D(2,2))/D(2,2))*P2(1,2)*P2(1,2)+...
    0.5*(D(1,1)-abs(D(1,1)))*0.5*(P1(1,1)*P2(2,2)+P1(1,2)*P2(2,1)+P1(2,2)*P2(1,1)+P1(2,1)*P2(1,2))/(D(1,1)-D(2,2))+0.5*(D(2,2)-abs(D(2,2)))*0.5*(P1(1,1)*P2(2,2)+P1(1,2)*P2(2,1)+P1(2,2)*P2(1,1)+P1(2,1)*P2(1,2))/(D(2,2)-D(1,1));
Q=[Q11 Q12 Q13;Q12 Q22 Q23;Q13 Q23 Q33];

%Q=reshape(Q,9,1);
%Q=[1;0;0;0;1;0;0;0;1] - P;
