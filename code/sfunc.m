function ds=sfunc(t,s)
% Function for Ricatti equation solution
global A B R Q
s=reshape(s,[8,8]);
ds_m=A'*s+s*A-s*B*inv(R)*B'*s+Q;
dss=ds_m;
ds=reshape(dss,[64,1]);
end