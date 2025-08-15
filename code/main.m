%% Simulation of stochastic control system
% F-16 lateral LQG controller design
% Sohrab Rezaei
% To run this code, the following function files must be located in the 
% same directory as this main script:
%   - xhatfunc
%   - xfunc
%   - sfunc
%   - pfunc
%
% The simulation uses the following state variables:
%   beta      = x(1);  % Sideslip angle
%   phi       = x(2);  % Bank angle
%   p         = x(3);  % Roll rate
%   r         = x(4);  % Yaw rate
%   sigma_a   = x(5);  % Aileron deflection
%   sigma_r   = x(6);  % Rudder deflection
% References:
%   [1] Lewis, F. L., Xie, L., & Popa, D. (2017). Optimal and robust estimation: 
%       with an introduction to stochastic control theory. CRC Press.
%   [2] Stevens, B. L., Lewis, F. L., & Johnson, E. N. (2015). Aircraft control and 
%       simulation: dynamics, controls design, and autonomous systems. John Wiley & Sons.
clc;
close all; 
clear;
global u i A B C G R Q Qp L yy Rp r
%% time
t_max=10;
step_size=0.001;
%% system 
A=[-0.3220,0.064,0.0364,-0.9917,0.0003,0.0008,0,0
    0,0,1,0.0037,0,0,0,0
    -30.6492,0,-3.6784,0.6646,-0.7333,0.1315,0,0
    8.5395,0,-0.0254,-0.4765,-0.0319,-0.062,0,0
    0,0,0,0,-20.2,0,-0.01,-5.47
    0,0,0,0,0,-20.2,-0.168,51.71
    0,0,0,0,0,0,0,0
    0,0,0,0,0,0,0,0];
B=[0,0,0,0,0,0,1,0
   0,0,0,0,0,0,0,1]';
C=[0,57.2958,0,0,0,0,0,0
   57.2958,0,0,0,0,0,0,0];
G=eye(8);
D=zeros(2,2);
%The system is controllable if Co has full rank n.
Co = ctrb(A,B);
Corank=rank(Co);
unco= length(A)- rank(Co);
 if (unco == 0)
    disp('system is Controllable')
 end
 %% initial values
init=[0 0 0 0 0 0 0 0];
r=[1;0];                % demand 
Rp=1e-5.*[1 0;0 1]; % measurement covariance
N=0.1*ones(2,8);        % E(wv')=N 
Qp=diag([0.1,0.1,0.1,0.1,0,0,1,1]); 
Q=C'*C;
R=1e-12*eye(2);         % cost of u
n=t_max/step_size+1;    
x1=init;
x1hat=zeros(1,8);
i=2;p=zeros(8,8);s=zeros(8,8);
L=zeros(8,2);
u=[0;0];

for t=0:step_size:t_max
  %% modified kalman filter
  v=sqrt(Rp)*randn(2,1); 
  pinit=reshape(p,[1,64]);
  [~,pp]=ode23(@pfunc,[t t+step_size],pinit);
  p=reshape(pp(end,:),[8,8]);
  L=(p*C'+G*N')*inv(Rp);
  yy=C*x1(end,:)'+v;
  output_y(i,:)=yy;
  [~,xxhat]=ode23(@xhatfunc,[t t+step_size],x1hat(end,:));
  x1hat(i,:)=xxhat(end,:);
%% Optimal control input
  sinit=reshape(s,[1,64]);
  [~,ss]=ode23(@sfunc,[t t+step_size],s);
  s=reshape(ss(end,:),[8,8]);
  k=inv(R)*B'*s;
  u=-k*x1hat(end,:)';
  c_out(i,:)=u;
[td,x]=ode23(@xfunc,[t t+step_size],x1(end,:));
  x1(i,:)=x(end,:); t1(i)=td(end);
  i=i+1;
end
% [x1,x1hat]=f16lat(tmax,x_init);
%% plots
figure (1);
plot(t1,output_y(:,1));
hold on
plot(t1,output_y(:,2));
legend('bank angle','sideslip');
title('F-16 LQG Controller')
xlabel('time')
ylabel('Magnitude (deg)')
grid

figure (2);
subplot(2,1,1)
plot(t1,x1(:,1))
hold on
plot(t1,x1hat(:,1),'-.')
title('sideslip')
xlabel('time')
ylabel('Magnitude (deg)')
grid
legend('actual','estimate');
subplot(2,1,2)
plot(t1,x1(:,2))
hold on
plot(t1,x1hat(:,2),'-.')
legend('actual','estimate');
title('bank angle')
xlabel('time')
ylabel('Magnitude (deg)')
grid

figure(3);
plot(t1,x1(:,5));
title('F-16 LQG Controller')
xlabel('time')
hold on
plot(t1,x1(:,6));
legend('aileron','rudder');
grid





