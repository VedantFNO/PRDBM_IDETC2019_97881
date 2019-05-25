%% simulation
tic
clear all
clc
close all
n = 4;

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
rectangular = false;
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
t = rand()/20
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
desired_omega = [48.342;359.4;1013.5;5711.4;7816];%
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

K_list = [] ;
Loss_list = [];
L_list = [];
Init = [];
for i = 0.1:0.01:(1*pi/(2*n))
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
lb = [0*ones(1,size(Input.x,2)) ,1/(2*space)*ones(1,(length(design_init)-size(Input.x,2)))];
ub = [inf*ones(1,size(Input.x,2)) ,((space-1)/space)*ones(1,(length(design_init)-size(Input.x,2)))];
design_init_list = repmat([ones(1,size(Input.x,2)), Ln, LCentern],40,1);
for qw=2:(256*4)
        
        lengths = diff(sort(randperm(space-1,n+1)))/space;
        lengths = [lengths 1-sum(lengths)];
        design_init_list(qw,:)= [rand(1,size(Input.x,2)), lengths];
    end
if (gradient)
    pctRunOnAll warning('off')
    opts = optimoptions('fmincon','Display','iter','ConstraintTolerance', 1e-6,'StepTolerance', 1e-6,'MaxIterations',10000,...
    'OptimalityTolerance',1e-6,'MaxFunctionEvaluations' , 5000 );%,'Algorithm','sqp'
    problem = createOptimProblem('fmincon','x0',design_init,'objective',@(x) findKL(x, Input,B,C,D,G),...
        'Aeq',[A_eq],'beq',[b_eq],'lb',lb,'ub',ub,'options',opts);
    ms = MultiStart('StartPointsToRun','bounds','UseParallel',true);
    tpoints = CustomStartPointSet(design_init_list);
[K_gard,Loss_gard] = run(ms,problem,{tpoints});


% [K_gard,Loss_gard] = fmincon(problem);
% [K,Loss] = fmincon(@(x) findKL(x, Input,B,C,D,G), design_init,[],[],[A_eq],[b_eq],[1e-3*ones(1,size(Input.x,2)) (1/(5*n))*ones(1,(length(design_init)-size(Input.x,2)))],[inf*ones(1,size(Input.x,2)) (5/(n))*ones(1,(length(design_init)-size(Input.x,2)))],[],opts);
end
if(GA)
    pctRunOnAll warning('off')
    hybridopts = optimoptions('fmincon','Display','iter','Algorithm','sqp','ConstraintTolerance', 1e-8,'StepTolerance', 1e-8,'MaxIterations',10000,...
                            'OptimalityTolerance',1e-8,'MaxFunctionEvaluations' , 10000 );%
    opts = optimoptions('ga','Display','iter','ConstraintTolerance', 1e-8,'HybridFcn',{@fmincon,hybridopts},'MaxStallGenerations',15,...
                    'InitialPopulationMatrix',design_init_list);%'PopulationSize',250'PlotFcn', @gaplotbestf,
%     problem = createOptimProblem('ga','objective',@(x) findKL(x, Input,B,C,D,G),'Aeq',[A_eq],'beq',[b_eq],'lb',lb,'ub',ub,...
%     'options',opts);
[K_ga,Loss_ga,exitflag,output,population,scores] = ga(@(x) findKL(x, Input,B,C,D,G),length(design_init),A,b,A_eq,b_eq,lb,ub,[],opts);
end
toc
%% Analysis of results

% [K_list_grad,L_list_grid] = Analysis_of_results(K_gard,Loss_gard,Input)
% [K_list_ga,L_list_ga] = Analysis_of_results(K_ga,Loss_ga,Input)

K = K_gard

Loss = Loss_gard

K;
% Loss_list = [Loss_list Loss];
K_list = [K_list; K(1:end-n-1)];
L_list = [L_list; [K(end) K(end-n:end-1)]];

Ln = K(end-Input.n:end-1); %length of each link
Lcenter = K(end);
L = 0;

L_tot(1) = Lcenter;
L_tot= [L_tot Ln];
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

%% Simulate oscillation 
% 
% qRn = linspace(0,1.0,n+1);
%     qRn = qRn(2:end);
%   
    qRn = 0*ones(1,n); % Joints of Right array
dqRn = -.5*ones(1,n); %Angular rates
    init = transpose([qRn,dqRn]);
    
    Design_Parametersn = [Mn,L_list(end,1:end-1),Jn,K_list_unscaled(1,:),MR_COM_xn,MR_COM_yn,L_list(end,end)];
% Design_Parametersn = [Mn,Ln,Jn,Kn,MR_COM_xn,MR_COM_yn,LCentern];

[t,xt] = ode113(@(t,x)dynamicsf(t,x,n,Design_Parametersn,B,C,D,G,J,0,0),[0 2],init);
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

