function A = operator_trace_plus(A)
A = ((A(1,:) + A(2,:)) + abs(A(1,:)+A(2,:)))/2;