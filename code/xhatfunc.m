function dxhat=xhatfunc(t,xhat)
% Function for estimate states (Kalman filter)
global A B C u L yy r
dxhat=A*xhat+B*u+L*(yy-C*xhat-r);
end