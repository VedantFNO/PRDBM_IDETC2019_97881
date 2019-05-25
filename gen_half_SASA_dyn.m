% Generate 2n-DOF planar dynamic model using Lagrangian method
% Author: Vedant
% Based on work of Prof. Hae-Won Park

%function [B,C,D,G] = gen_half_SASA_dyn(n)

clear
close all
clc
n=4;
qR = sym('qR', [1,n], 'real'); % Joints of Right array
dqR = sym('dqR', [1,n], 'real'); %Angular rates
ddqR = sym('ddqR', [1,n], 'real'); %Angular rates
M = sym('M', [1,n], 'real'); %mass of each link
L = sym('L', [1,n], 'real'); %length of each link
J = sym('J', [1,n], 'real'); %inertia of each link
K = sym('K', [1,n], 'real'); %inertia of each link
MR_COM_x = sym('Mr_com_x', [1,n], 'real'); % COM position for each link
MR_COM_y = sym('Mr_com_y', [1,n], 'real');

% % For Center Body
syms LCenter JCenter MCenter real % center body mass, intertia and length
syms MB_COM_x MB_COM_y real
syms qCenter dqCenter real
% syms xCenter dxCenter ddxCenter real
% syms yCenter dyCenter ddyCenter real

% if(stance_flag_fore)
q = [qCenter qR];
dq = [dqCenter dqR];

Design_Parameters = [M,L,J,K,MR_COM_x,MR_COM_y,LCenter,JCenter,MCenter,MB_COM_x,MB_COM_y];
Control_States = [q,dq];
% end
disp('Initialization done')

%% forward kinematics
rot_2d = @(theta) [cos(theta) -sin(theta); sin(theta) cos(theta)]; %counterclockwise rotation matrix for [x;y] vectors

abs_qR = sym('abs_qR',[1 n]);

abs_dqR = sym('abs_dqR',[1 n]);

for(i=1:n)
    %absolute angle
    abs_qR(i) = sum(qR(1:i));
    abs_dqR(i) = sum(dqR(1:i));
end

pFR = sym('pFR',[2 n+1]);

pFR(:,1) = rot_2d(qCenter)*[LCenter;0];
for(i=2:n+1)
    pFR(:,i) = [(pFR(1,(i-1)));(pFR(2,(i-1)))]+rot_2d(abs_qR(i-1))*[L(i-1);0];
end

% calculate COM position for each link
pFR_COM = sym('pFR_COM',[2 n]);
pFB_COM = [0;0]+rot_2d(qCenter)*[ MB_COM_x;MB_COM_y];% body
for(i=1:n)
    pFR_COM(:,i) = [LCenter;0]+rot_2d(abs_qR(i))*[ MR_COM_x(i);MR_COM_y(i)];
    if(i>1)
        pFR_COM(:,i) = pFR_COM(:,i)+pFR(:,(i-1));
    end
end
disp('Sat frame done')

%% Velocity kinematics

for(i=1:n)
    VR_COM(:,i)=jacobian(pFR_COM(:,i),q)*dq';
    
end

%% Kinetic Energy, Potential Energy, and Lagrangiandq2abs

for(i=1:n)
    KE_R(i) = 0.5*M(i)*VR_COM(:,i).'*VR_COM(:,i)+0.5*J(i)*abs_qR(i)^2;
    PE_R(i) = 0.5*K(i)*(qR(i)^2);
end

KE = sum(KE_R);
PE = sum(PE_R);


Upsilon = [qR] ; % where control torques go
dimensionparameters = [L,LCenter];
[D, C, G, B] = std_dynamics(KE,PE,q,dq, Upsilon);

matlabFunction(B,'File','B','Vars',{Design_Parameters ,Control_States})
matlabFunction(C,'File','C','Vars',{Design_Parameters ,Control_States})
matlabFunction(D,'File','D','Vars',{Design_Parameters ,Control_States})
matlabFunction(G,'File','G','Vars',{Design_Parameters ,Control_States})

matlabFunction(pFR,'File','visualize','Vars',{q,dimensionparameters})
%% simulation
qRn = 0.2*ones(1,n); % Joints of Right array
dqRn = 00*ones(1,n); %Angular rates
Mn = 0.1*ones(1,n); %mass of each link
Ln = ones(1,n); %length of each link
Jn = 0.1*ones(1,n); %inertia of each link
Kn = 0.5*ones(1,n); %inertia of each link
MR_COM_xn = ones(1,n); % COM position for each link
MR_COM_yn = ones(1,n);

MB_COM_xn =1;
MB_COM_yn =1;
qCentern =0; 
dqCentern =0;

% % For Center Body
LCentern = 1;
JCentern = 10;
MCentern = 10;

dimensionparametersn = [Ln,LCentern];
Design_Parametersn = [Mn,Ln,Jn,Kn,MR_COM_xn,MR_COM_yn,LCentern,JCentern,MCentern,MB_COM_xn,MB_COM_yn];
init = transpose([qCentern,qRn,dqCentern,dqRn]);

[t,xt] = ode113(@(t,x)dynamics_half(t,x,n,Design_Parametersn),[0 1000],init);
dt = diff(t)

h = animatedline('Marker','o');
axis([-2 6 -5 5])
axis square

for (i=1:length(xt))
    p = visualize(xt(i,1:n+1),dimensionparametersn);
    x = [0 p(1,:)];
    y = [0 p(2,:)];
    clearpoints(h);
    addpoints(h,x,y);
     pause(0.001);
end