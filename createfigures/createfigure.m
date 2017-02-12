function createfigure(X1, YMatrix1, X2, YMatrix2, X3, YMatrix3)
%CREATEFIGURE(X1,YMATRIX1,X2,YMATRIX2,X3,YMATRIX3)
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data
%  X2:  vector of x data
%  YMATRIX2:  matrix of y data
%  X3:  vector of x data
%  YMATRIX3:  matrix of y data

%  Auto-generated by MATLAB on 06-May-2016 16:33:23

% Create figure
figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1,'YScale','log','YMinorTick','on',...
    'YMinorGrid','on',...
    'XScale','log',...
    'XMinorTick','on',...
    'XMinorGrid','on',...
    'FontSize',12,...
    'FontName','Arial',...
    'xlim',[10^5 10^6],...
    'ylim',[10^-6 10^1]);
% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[100000 1000000]);
% Uncomment the following line to preserve the Y-limits of the axes
% ylim(axes1,[1e-06 10]);
box(axes1,'on');
grid(axes1,'on');
hold(axes1,'all');

% Create multiple lines using matrix input to loglog
loglog1 = loglog(X1,YMatrix1,'Parent',axes1,'MarkerSize',12,'LineWidth',1.5,...
    'Color',[0 0 1]);
set(loglog1(1),'Marker','o','DisplayName','Singlerate (i=V)');
set(loglog1(2),'Marker','square','DisplayName','Singlerate (i=K\_A)');
set(loglog1(3),'Marker','diamond','DisplayName','Singlerate (i=Ca)');

% Create multiple lines using matrix input to loglog
loglog2 = loglog(X2,YMatrix2,'Parent',axes1,'MarkerFaceColor',[0 0 0],...
    'MarkerSize',12,...
    'LineWidth',1.5,...
    'Color',[0 0 1]);
set(loglog2(1),'Marker','o','DisplayName','Multirate (i=V)');
set(loglog2(2),'Marker','square','DisplayName','Multirate (i=K\_A)');
set(loglog2(3),'Marker','diamond','DisplayName','Multirate (i=Ca)');

% Create xlabel
xlabel('Number of ODE calls');

% Create ylabel
ylabel('\epsilon_{i} [%]');

% Create legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.645089285714282 0.66815476190476 0.290178571428571 0.333630952380952]);
% Create multiple lines using matrix input to loglog
loglog(X3,YMatrix3,'Parent',axes1,'LineWidth',1.5,'Color',[0 0 1],'LineStyle','--','Color','b');
