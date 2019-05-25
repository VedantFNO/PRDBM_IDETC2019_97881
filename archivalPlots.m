%% Results
clear all

% Uniform
% rectangular = true;
% K_list = [2.7441 2.8032 2.8830 2.9815 3.0960 3.2246 3.3652 3.5164 3.6776];
% 
% L_list =[ 0.7447 0.0466 0.0725 0.1362];

% K_list = [3.2301 3.2812 3.3504 3.4356 3.5349 3.6462];
% 
% L_list =[ 0.7555 0.0369 0.1323 0.0281 0.0472];

% % Random
rectangular = false;
% K_list = [6.5099 6.7109 6.9627 7.2464 7.5471 7.8553];
% 
% L_list =[0.3553 0.0707 0.1242 0.1819 0.2678];

% K_list = [5.0093 5.2555 5.5703 5.9349 6.3344 6.7601];
%     
% L_list =[0.4345 0.0890 0.0447 0.1366 0.2952];
% 
% K_list = [ 2.1305 2.1365 2.1450 2.1557 2.1687];
% 
% L_list =[0.6323 0.1930 0.0667 0.0435 0.0418 0.0228];

% K_list = [5.0093 5.2555 5.5703 5.9349 6.3344 6.7601];
% 
% L_list =[0.4345 0.0890 0.0447 0.1366 0.2952];

%new test--

K_list = 0*[5.0093 5.2555 5.5703 5.9349 6.3344 6.7601];

L_list =ones(1,5)/5;

%% plots!
n = length(L_list)-1

K = [K_list L_list(2:end) L_list(1)]
cd SimMatrixFun
f_names = dir (strcat('J',num2str(n),'.m'));
if(isempty(f_names))
cd ..
gen_arm_dyn_fun(n)
cd SimMatrixFun
end


myxlabel = '$x [m]$'; % x label with latex
myylabel = '$y [m]$'; % y label with latex
myzlabel = '$z [m]$'; % x label with latex
mylegend = ''; % legend with latex

% plot variables
fwidth = 550; % figure width in pixels
fheight = 400; % figure height in pixels
fontlabel = 12; % x,y label font size
fontlegend = 12; % x,y legend font size
fonttick = 12; % x,y rick font size
mycolor1 = [0.8500 0.3250 0.0980]; % custom color 1
mycolor2 = [0.4940 0.1840 0.5560]; % custom color 2
wcolor = [1 1 1]; % white color
bcolor = [0 0 0]; % black color

set(0,'DefaultTextInterpreter','latex'); % change the text interpreter
set(0,'DefaultLegendInterpreter','latex'); % change the legend interpreter
set(0,'DefaultAxesTickLabelInterpreter','latex'); % change the tick interpreter
hf = figure; % create a new figure and save handle
%hf.Color = wcolor; % change the figure background color



B = str2func(strcat('B',num2str(n)));
C = str2func(strcat('C',num2str(n)));
D = str2func(strcat('D',num2str(n)));
G = str2func(strcat('G',num2str(n)));
J = str2func(strcat('J',num2str(n)));
visualize = str2func(strcat('visualize',num2str(n)));
cd ..
random = false;
gradient = true;
GA = false;

choice = 17572%randi(50000)%3904,17572
%% generate random panel
if(random)
sections = 6
        X = [];
        if sections>2
            X = sort(rand(1,sections-2));
        end
        X = [0,X,1];
        if sections==0
            Y = rand()*ones(1,length(X));
        else
            Y = rand(1,length(X));
        end
        t = 1*rand();
        %         vars = gen_SASA_model_params(X,Y,t,n,samples);
        %vars = double(subs(vars_sym,[Xs Ys ts],[X Y t]));
        
%% read panel from lookup
elseif(rectangular)
    
X = [0,1]
% Y = 0.5*rand();
Y = 0.1812;
Y = [Y,Y];
t = 0.040532
else
data = csvread('lookup.csv',choice,0,[choice 0 choice 27]);
leng = (data(1));
X = data(3:3+leng-1);
Y = data(13:13+leng-1)*0.5;
t = data(23)/20;
[q,w] = sort(Y);
end

%% define desired props
qRn = 0*ones(1,n); % Joints of Right array
dqRn = -0.0*ones(1,n); %Angular rates

% % For Center Body
b = 2*Y(1);
factor = sqrt((205e9*(b*t^3/12))/(7850*b*t))
Ln = ones(1,n); %length of each link
LCentern = 1;
L_total = [LCentern,Ln];
L_total = L_total/sum(L_total);

base_freq = 2;
desired_omega = [(1.875^2);(4.694^2);(7.885^2)]*factor;%;

if(~rectangular)
candidate = X(w);
candidate= candidate(X(w)<1 & X(w)>0);
L_total = [0 sort(candidate(1:n)) 1];
L_total = diff(L_total);
desired_omega = [24674;189760;551810;];%3898500;5469100
end
LCentern = L_total(1);
Ln = L_total(2:end);
figure()
[mass,  COM_x,  COM_y, ~, Iyy]  = gen_SASA_model_params_poly(X,Y,t,L_total);

Mn = mass; %mass of each link
Jn = Iyy; %inertia of each link
Kn = 0*ones(1,n); %inertia of each link
MR_COM_xn = COM_x; % COM position for each link
MR_COM_yn = COM_y;
K = 0;

Design_Parametersn = [Mn,Ln,Jn,Kn,MR_COM_xn,MR_COM_yn,LCentern];
init = transpose([qRn,dqRn]);




Input.omega = desired_omega;
Input.DParams = Design_Parametersn;
Input.x = init;
Input.n = n;

%% Start optimizer
clc


Loss_list = [];
Init = [];
for i = 0.1:0.05:(1*pi/(2*n))
    qRn = linspace(0,i,n+1);
    %qRn = qRn(2:end);
    qRn = i*ones(1,n);
    init = transpose([qRn,dqRn]);
%     Input.x = init;
    Init = [Init init];
end
Input.x = Init;
Input.X = X;
Input.Y = Y;
Input.t = t;

design_init = [1e-3*ones(1,size(Input.x,2)), Ln, LCentern];
A = [];
b = [];
A_eq = zeros(length(design_init));
A_eq(1,end - n:end) = 1;
b_eq = zeros(1,length(A_eq));
b_eq(1) = 1;
space = 20;
lb = [0*ones(1,size(Input.x,2)) ,1/space*ones(1,(length(design_init)-size(Input.x,2)))];
ub = [inf*ones(1,size(Input.x,2)) ,((space-1)/space)*ones(1,(length(design_init)-size(Input.x,2)))];
design_init_list = repmat([ones(1,size(Input.x,2)), Ln, LCentern],40,1);
for qw=2:(256*4)
        
        lengths = diff(sort(randperm(space-1,n+1)))/space;
        lengths = [lengths 1-sum(lengths)];
        design_init_list(qw,:)= [rand(1,size(Input.x,2)), lengths];
    end

L_list =ones(1,5)/5;

K = [K_list L_list(2:end) L_list(1)]

%% all plots

Ln = K(end-Input.n:end-1); %length of each link
Lcenter = K(end);
L = 0;
clear L_tot
L_tot(1) = Lcenter;
L_tot= [L_tot Ln];
figure(40)
subplot(2,1,1)
[mass,  COM_x,  COM_y, Ixx, Iyy]  = gen_SASA_model_params_poly(Input.X,Input.Y,Input.t,L_tot);

Mn = mass(2:end); %mass of each link
Jn = Iyy(2:end); %inertia of each link
MR_COM_xn = COM_x(2:end); % COM position for each link
MR_COM_yn = COM_y(2:end);

for i=1:size(Input.x,2)
K_diag = ones(1,(Input.n));
for k = 1:(Input.n)
    K_diag(k) = sqrt(200e9*Jn(k)/Mn(k));
end
K_list_unscaled(i,:) = K(i)*K_diag;
end

K_list_unscaled;
% K_list
L_list(end,:)
% Loss_list;

dimensionparametersn = [L_list(end,2:end),L_list(end,1)];
%Kn =    1.0 * [   20.3989,17.9315]
figure()
set(gca,'FontSize',32)
yyaxis left
plot(rad2deg(Init(n,:)),K_list,'bd')
xlabel('Deflection angle($\theta_i$) [degrees]','FontSize',22,'interpreter','latex')%
ylabel('Stiffness correction factor ($k_i$)','FontSize',22,'interpreter','latex')%

yyaxis right
plot(rad2deg(Init(n,:)),(K_list-mean(K_list))*100/mean(K_list),'rd')
ylabel('Percent change from mean','FontSize',22,'interpreter','latex')%




%% Simulate oscillation 
% 
% qRn = linspace(0,1.0,n+1);
%     qRn = qRn(2:end);
%   
    qRn = 0.5*ones(1,n); % Joints of Right array
dqRn = -.0*ones(1,n); %Angular rates
    init = transpose([qRn,dqRn]);
    
    Design_Parametersn = [Mn,L_list(end,1:end-1),Jn,K_list_unscaled(1,:),MR_COM_xn,MR_COM_yn,L_list(end,end)];
% Design_Parametersn = [Mn,Ln,Jn,Kn,MR_COM_xn,MR_COM_yn,LCentern];

[t,xt] = ode113(@(t,x)dynamicsf(t,x,n,Design_Parametersn,B,C,D,G,J,0,0),[0 100],init);
figure()
h = animatedline('Marker','o');
axis([0 1 -1 1])
axis square
dt = diff(t);
for (i=1:length(xt))
    p = visualize(xt(i,1:n),dimensionparametersn);
    x = [0 p(1,:)];
    y = [0 p(2,:)];
        clearpoints(h);
    addpoints(h,x,y);
    x = [p(1,1:end-1)];
    y = [p(2,1:end-1)];
    addpoints(h,x,y);

    if(i ~= length(xt))
        pause(dt(i));
    end
end
figure()
plot(rad2deg(Init(n,:)),K_list,'.b')
figure(4)
plot(t,xt(:,end))
xlim([0 1])

%% Simulate deflection   
  
    qRn = 0*ones(1,n); % Joints of Right array
dqRn = 0*ones(1,n); %Angular rates
    init = transpose([qRn,dqRn]);
    
    Design_Parametersn = [Mn,L_list(end,1:end-1),Jn,Kn,MR_COM_xn,MR_COM_yn,L_list(end,end)];
% Design_Parametersn = [Mn,Ln,Jn,Kn,MR_COM_xn,MR_COM_yn,LCentern];

[t,xt] = ode113(@(t,x)dynamicsf(t,x,n,Design_Parametersn,B,C,D,G,J,0,-1),[0 5000],init);


figure()
h = animatedline('Marker','o');
axis([0 1 -1 1])
axis square
dt = diff(t);
for (i=(length(xt)-50):length(xt))
% i = length(xt)
    p = visualize(xt(i,1:n),dimensionparametersn);
    x = [0 p(1,:)];
    y = [0 p(2,:)];
    clearpoints(h);
    addpoints(h,x,y);
    if(i ~= length(xt))
        pause(dt(i));
    end
end

figure()
plot(t,xt(:,1:3))
xlim([t(end-50) t(end)])

K_list
L_list
K_list_unscaled