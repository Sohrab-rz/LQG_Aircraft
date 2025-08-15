function xdot=xfunc(t,x)
% Function for plant dynamics 
global u A B G Qp
%       beta=x(1);      sideslip
%       phi=x(2);       bank angle
%       p=x(3);         roll rate
%       r=x(4);         yaw rate
%       sigma_a=x(5);   aileron deflection
%       sigma_r=x(6);   rudder  deflection
w=0.01*sqrt(Qp)*rand(8,1);
xdot=(A*x+B*u+G*w);
end