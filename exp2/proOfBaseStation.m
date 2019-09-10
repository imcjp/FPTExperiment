function [ z ] = proOfBaseStation(P,r0,d1,rm,lamda0)
v1=P*d1*lamda0;
if r0<=rm-d1
    z=4*v1*((d1*d1+r0*r0)/rm^4-1/rm^2)*exp(v1*d1*(d1*d1+2*r0*r0)/rm^4-2*v1*d1/rm^2);
else
    z=-2*P*lamda0*exp(-2*P*lamda0*(atan((r0^2 - d1^2 + rm^2)/(2*r0*(rm^2 - (r0^2 - d1^2 + rm^2)^2/(4*r0^2))^(1/2)))/pi + d1^2*(2/(pi*rm^2) - (d1^2 + 2*r0^2)/(pi*rm^4))*(pi/2 - atan((r0 - (r0^2 - d1^2 + rm^2)/(2*r0))/(d1^2 - (r0 - (r0^2 - d1^2 + rm^2)/(2*r0))^2)^(1/2))) - (2*(r0 - (r0^2 - d1^2 + rm^2)/(2*r0))^4 - 7*d1^2*(r0 - (r0^2 - d1^2 + rm^2)/(2*r0))^2 + 5*d1^4)/(6*pi*rm^4*(d1^2 - (r0 - (r0^2 - d1^2 + rm^2)/(2*r0))^2)^(1/2)) - ((2/(pi*rm^2) - ((r0^2 - d1^2 + rm^2)^2/(2*r0^2) + rm^2)/(3*pi*rm^4))*(rm^2 - (r0^2 - d1^2 + rm^2)^2/(4*r0^2))^(1/2)*(r0^2 - d1^2 + rm^2))/(2*r0) - ((d1^2 - (r0 - (r0^2 - d1^2 + rm^2)/(2*r0))^2)^(1/2)*(2*r0 - (r0^2 - d1^2 + rm^2)/r0))/(pi*rm^2) + ((d1^2 - (r0 - (r0^2 - d1^2 + rm^2)/(2*r0))^2)^(1/2)*(d1^2*(13*r0 + (3*(r0^2 - d1^2 + rm^2))/(2*r0)) + (r0^2 - d1^2 + rm^2)^2/(2*r0) - (3*(r0^2 - d1^2 + rm^2)^3)/(4*r0^3) + 2*r0^2 + r0*(r0^2 - d1^2 + rm^2)))/(6*pi*rm^4) + 1/2))*((2*d1^3*(pi/2 - atan((r0 - (r0^2 - d1^2 + rm^2)/(2*r0))/(d1^2 - (r0 - (r0^2 - d1^2 + rm^2)/(2*r0))^2)^(1/2))))/(pi*rm^4) - 2*d1*(2/(pi*rm^2) - (d1^2 + 2*r0^2)/(pi*rm^4))*(pi/2 - atan((r0 - (r0^2 - d1^2 + rm^2)/(2*r0))/(d1^2 - (r0 - (r0^2 - d1^2 + rm^2)/(2*r0))^2)^(1/2))) - (14*d1*(r0 - (r0^2 - d1^2 + rm^2)/(2*r0))^2 - 20*d1^3 - (8*d1*(r0 - (r0^2 - d1^2 + rm^2)/(2*r0))^3)/r0 + (14*d1^3*(r0 - (r0^2 - d1^2 + rm^2)/(2*r0)))/r0)/(6*pi*rm^4*(d1^2 - (r0 - (r0^2 - d1^2 + rm^2)/(2*r0))^2)^(1/2)) - (d1/(r0*(rm^2 - (r0^2 - d1^2 + rm^2)^2/(4*r0^2))^(1/2)) + (d1*(r0^2 - d1^2 + rm^2)^2)/(4*r0^3*(rm^2 - (r0^2 - d1^2 + rm^2)^2/(4*r0^2))^(3/2)))/(pi*((r0^2 - d1^2 + rm^2)^2/(4*r0^2*((r0^2 - d1^2 + rm^2)^2/(4*r0^2) - rm^2)) - 1)) + (d1^2*(((r0 - (r0^2 - d1^2 + rm^2)/(2*r0))*(2*d1 - (2*d1*(r0 - (r0^2 - d1^2 + rm^2)/(2*r0)))/r0))/(2*(d1^2 - (r0 - (r0^2 - d1^2 + rm^2)/(2*r0))^2)^(3/2)) - d1/(r0*(d1^2 - (r0 - (r0^2 - d1^2 + rm^2)/(2*r0))^2)^(1/2)))*(2/(pi*rm^2) - (d1^2 + 2*r0^2)/(pi*rm^4)))/((r0 - (r0^2 - d1^2 + rm^2)/(2*r0))^2/((r0 - (r0^2 - d1^2 + rm^2)/(2*r0))^2 - d1^2) - 1) - (d1*(2/(pi*rm^2) - ((r0^2 - d1^2 + rm^2)^2/(2*r0^2) + rm^2)/(3*pi*rm^4))*(rm^2 - (r0^2 - d1^2 + rm^2)^2/(4*r0^2))^(1/2))/r0 + ((d1^2 - (r0 - (r0^2 - d1^2 + rm^2)/(2*r0))^2)^(1/2)*(2*d1*r0 - 2*d1*(13*r0 + (3*(r0^2 - d1^2 + rm^2))/(2*r0)) + (3*d1^3)/r0 - (9*d1*(r0^2 - d1^2 + rm^2)^2)/(2*r0^3) + (2*d1*(r0^2 - d1^2 + rm^2))/r0))/(6*pi*rm^4) - ((2*d1 - (2*d1*(r0 - (r0^2 - d1^2 + rm^2)/(2*r0)))/r0)*(2*(r0 - (r0^2 - d1^2 + rm^2)/(2*r0))^4 - 7*d1^2*(r0 - (r0^2 - d1^2 + rm^2)/(2*r0))^2 + 5*d1^4))/(12*pi*rm^4*(d1^2 - (r0 - (r0^2 - d1^2 + rm^2)/(2*r0))^2)^(3/2)) + (d1*(2/(pi*rm^2) - ((r0^2 - d1^2 + rm^2)^2/(2*r0^2) + rm^2)/(3*pi*rm^4))*(r0^2 - d1^2 + rm^2)^2)/(4*r0^3*(rm^2 - (r0^2 - d1^2 + rm^2)^2/(4*r0^2))^(1/2)) + (2*d1*(d1^2 - (r0 - (r0^2 - d1^2 + rm^2)/(2*r0))^2)^(1/2))/(pi*r0*rm^2) + ((2*r0 - (r0^2 - d1^2 + rm^2)/r0)*(2*d1 - (2*d1*(r0 - (r0^2 - d1^2 + rm^2)/(2*r0)))/r0))/(2*pi*rm^2*(d1^2 - (r0 - (r0^2 - d1^2 + rm^2)/(2*r0))^2)^(1/2)) - ((2*d1 - (2*d1*(r0 - (r0^2 - d1^2 + rm^2)/(2*r0)))/r0)*(d1^2*(13*r0 + (3*(r0^2 - d1^2 + rm^2))/(2*r0)) + (r0^2 - d1^2 + rm^2)^2/(2*r0) - (3*(r0^2 - d1^2 + rm^2)^3)/(4*r0^3) + 2*r0^2 + r0*(r0^2 - d1^2 + rm^2)))/(12*pi*rm^4*(d1^2 - (r0 - (r0^2 - d1^2 + rm^2)/(2*r0))^2)^(1/2)) + (d1*(rm^2 - (r0^2 - d1^2 + rm^2)^2/(4*r0^2))^(1/2)*(r0^2 - d1^2 + rm^2)^2)/(3*pi*r0^3*rm^4));
end
z=abs(z);
end
