function [c,w]=cuad_T0

%nodos de cuadratura
c1=[1 1]/3;
c2=(6+sqrt(15))*[1 1]/21;
c3=[(9-2*sqrt(15)),(6+sqrt(15))]/21;
c4=[c3(2) c3(1)];
c5=(6-sqrt(15))*[1 1]/21;
c6=[(9+2*sqrt(15)),(6-sqrt(15))]/21;
c7=[ c6(2) c6(1)];
c=[c1;c2;c3;c4;c5;c6;c7]';
cx=c(1,:);
cy=(c(2,:));

% pesos en los nodos de cuadratura
w1=0.1125;
w2=(155+sqrt(15))/2400;
w3=w2;
w4=w2;
w5=(155-sqrt(15))/2400;
w6=w5;
w7=w5;

w=[w1 w2 w3 w4 w5 w6 w7]';

end