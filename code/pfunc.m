function dp=pfunc(t,p)
% Function for Ricatti equation solution
global  A C G  Rp Qp
p=reshape(p,[8,8]);
ddp=A*p+p*A'+G*Qp*G'-p*C'*inv(Rp)*C*p;
dp=reshape(ddp,[64,1]);
end