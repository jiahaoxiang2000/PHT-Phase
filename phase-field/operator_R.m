function [B,C] = operator_R(A)
B = (sign(A(1,:)+A(2,:))+1)./2; % R+
C = (sign(-A(1,:)-A(2,:))+1)./2; % R-