function [mass,  COM_x,  COM_y, Ixx, Iyy]  = gen_SASA_model_params_poly(X,Y,t,Ln)
% clear all
% clc
colorblind =[0         0    1.0000;...
    1.0000         0         0;...
    1.0000    1.0000         0;...
    0.6602    0.6602    0.6602;...
         0         0         0;...
    1.0000    0.6445         0;...
    1.0000         0    1.0000;...
         0    0.5000    0.5000;...
         0         0    0.5430;...
         0    0.3906         0;...
         0    1.0000    1.0000;...
    0.5977    0.1953    0.7969];
    
close all
%% Analytical (Polygonal only)
% f = figure('visible','off');
% set(f,'color','none')
% patch([0 X 1],[0 Y 0], [1,1-t,t], 'LineStyle','None' );
% hold on;
% patch([0 X 1],[0 -Y 0], [1,1-t,t], 'LineStyle','None' );
% xlim([-0.1 1.1]);
% filename = strcat('fig_',num2str(n,'%08.f'),'_',num2str(samples,'%08.f'),'.png')
% saveas(f,[pwd '\Figures\' filename])
nodes(1) = 0;
mu = 7850;
for k=1:length(Ln)
    nodes(k+1) = nodes(k)+Ln(k);
end
clear xy;
xy = [X' Y'];
xy = [0 0;xy];
xy = [xy;[(fliplr(X))' (fliplr(-Y))']];
xy = [xy;[0 0]];


% xlimit = [0.25 0.65];
% ylimit = [-2  2];
% xbox = xlimit([1 1 2 2 1]);
% ybox = ylimit([1 2 2 1 1]);
% mapshow(xbox,ybox,'DisplayType','polygon','LineStyle','none')


% mapshow(xy(:,1),xy(:,2),'Marker','+')
% 
% [xi,yi] = polyxpoly(xy(:,1),xy(:,2),xbox,ybox);
% mapshow(xi,yi,'DisplayType','point','Marker','o')
P = polyshape(xy(:,1),xy(:,2));
% x = P.Vertices(P.Vertices(:,2)>0,1)
% y = P.Vertices(P.Vertices(:,2)>0,2)
% area(x,y,'LineStyle','none','FaceColor',(0.5*[0 0 1]+0.5*ones(1,3)))%,'LineWidth',2
% hold on
% area(x,-y,'LineStyle','none','FaceColor',(0.5*[0 0 1]+0.5*ones(1,3)))%,'LineWidth',2
% set(gca,'FontSize',24)
% xlabel('Length [m]')%,'FontSize',20
% ylabel('Length [m]')%,'FontSize',20
% axis([-0.01 1.01 -0.5 0.5])
% line([-10,10],[0,0],'Color','k')
for i = 1:(length(nodes)-1)

% if i==1
% pgon = polyshape([0 0 nodes(i) nodes(i)],[3 -3 -3 3]);
% lengths = nodes(i)-0;
% elseif i == (length(nodes)+1)
% pgon = polyshape([nodes(i-1) nodes(i-1) 1 1],[3 -3 -3 3]);
% lengths = 1-nodes(i-1);
% else
pgon = polyshape([nodes(i) nodes(i) nodes(i+1) nodes(i+1)],[2 -2 -2 2]);
lengths = nodes(i+1)-nodes(i);
a = nodes(i);
b = nodes(i+1);
% end
% ------------
% plot(pgon)
hold on
polyout = intersect(P,pgon);
% plot(polyout,'LineWidth',1.5)

    
% plot(nodes(i),0,'*k')
% fill(polyout.Vertices(:,1),polyout.Vertices(:,2),(0.5*colorblind(i,:)+0.5*ones(1,3)))%,'LineWidth',2
x = polyout.Vertices(polyout.Vertices(:,2)>0,1);
y = polyout.Vertices(polyout.Vertices(:,2)>0,2);
<<<<<<< HEAD
area(x,y,'LineStyle','none','FaceColor',(0.5*colorblind(i,:)+0.5*ones(1,3)))%,'LineWidth',2
area(x,-y,'LineStyle','none','FaceColor',(0.5*colorblind(i,:)+0.5*ones(1,3)))%,'LineWidth',2
set(gca,'FontSize',20)
xlabel('Length [m]','FontSize',20)%
ylabel('Width [m]','FontSize',20)%
=======
area(x,y,'LineStyle','none','FaceColor',(0.5*colorblind(i,:)+0.5*ones(1,3)));%,'LineWidth',2
area(x,-y,'LineStyle','none','FaceColor',(0.5*colorblind(i,:)+0.5*ones(1,3)));%,'LineWidth',2
set(gca,'FontSize',24);
xlabel('Length [m]');%,'FontSize',20
ylabel('Length [m]');%,'FontSize',20
>>>>>>> e37170208c9459dd85b67786de0bdeb26227442b


% xlim([0 1.2])
% ylim([-1 1])

axis([-0.05 1.05 -1.0 0.5])
% axis equal
% ----------
polyout.Vertices(:,1) = polyout.Vertices(:,1)-min(polyout.Vertices(:,1));
xy_part = polyout.Vertices;
props = PolygonMoments(xy_part,[],0);
lengths_list(i) = lengths;
mass(i) =  abs(props.Area*t)*1*mu;
COM_y(i) = props.MAy/props.Area;
COM_x(i) = props.MAx/props.Area;
Ixx(i) =   abs(props.Ixx*t*mu);
Iyy(i) =   abs(props.Iyy*t*mu);
str = {("Section"+" "+dec2rom(i)),strcat("Mass: ",num2str(mass(end))),  strcat("$COM_x$: ",num2str(COM_x(end))),  strcat("$COM_y$: ",num2str(COM_y(end))), strcat("$I_{xx}$: ",num2str(Ixx(end))), strcat("$I_{yy}$: ",num2str(Iyy(end)))};
text(a+(COM_x(end)),-0.6,str,'FontSize',14,'HorizontalAlignment','center');%,'FontWeight','bold'


% SASA = [mass,  COM_y,  COM_x, Ixx, Iyy];
% gen_SASA_model_params_poly(X,Y,0.1,[0.1 0.2 0.5 0.7])
end
end
function ans = dec2rom(z)
d = [ 1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1];
c =  {'M', 'CM', 'D', 'CD', 'C', 'XC', 'L', 'XL', 'X', 'IX', 'V', 'IV', 'I'};
[];
for ii = 1:numel(d)
    if z >= d(ii)    
        ans = [ans,repmat(c{ii},1,fix(z/d(ii)))];
        z = rem(z,d(ii));
    end
end
end