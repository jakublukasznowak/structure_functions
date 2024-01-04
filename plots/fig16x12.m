
function [f,ax,co,ls,mk] = fig16x12(scale,mgrid,varargin)

% FIG16x12 prepares a single figure optimized for printing purposes.
%
% [F,AX] = fig16x12(SCALE, MGRID) provides the handles to the figure F
% and to the axes AX with position and margins optimized for printing in
% 16 x 12 cm size. 
% SCALE controls the scale of the horizontal and vertical axes and can be
% one of the following: 'linlin','loglin','linlog','loglog', where the
% abrreviations correspond to linear and logarithmic scale.
% MGRID controls the grid appearance and needs to be a two element logical
% vector where the first element enables/disables the major grid and the
% second element enables/disables the minor grid.
%
% [...,CO,LS,MK] = fig16x12 additionally outputs color order CO, line style
% order LS, marker style order MK which are useful to control the styles
% while plotting multiple data series in a single figure.
%
% See also FIGURE, AXES.


width=16;
height=12;
mx=0.14;
my=0.15;
fontsize=14;

axpos=[mx my 1-1.5*mx 1-1.5*my];


co=repmat([0  0.4470    0.7410;
    0.8500    0.3250    0.0980;
    0.9290    0.6940    0.1250;
    0.4940    0.1840    0.5560;
    0.4660    0.6740    0.1880;
    0.3010    0.7450    0.9330;
    0.6350    0.0780    0.1840;
    1 0 0;
    0 1 0;
    0 0 1;
    1 1 0;
    1 0 1;
    0 1 1;
    0 0 0;
    0.7 0.7 0.7],10,1);
ls=repmat({'-','--','-.',':'},1,5);
mk=repmat({'o','^','d','s'},1,7);


if nargin<1, scale=''; end
if nargin<2 || isempty(mgrid), mgrid=[0 0]; end

if strcmp(scale,'loglog')
    xscale='log'; yscale='log';
elseif strcmp(scale,'linlog')
    xscale='linear'; yscale='log';
elseif strcmp(scale,'loglin')
    xscale='log'; yscale='linear';
else
    xscale='linear'; yscale='linear';
end

if mgrid(1), xmgrid='on'; else, xmgrid='off'; end
if mgrid(2), ymgrid='on'; else, ymgrid='off'; end


f = figure('Color','white','PaperUnits','centimeters',...
    'PaperSize',[width height],'PaperPosition',[0 0 width height]);

ax = axes('Parent',f,'Position',axpos,...
    'Color','none','FontSize',fontsize,'Box','on',...
    'XGrid','on','YGrid','on','GridAlpha',0.2,...
    'XMinorGrid',xmgrid,'YMinorGrid',ymgrid,'MinorGridAlpha',0.5,...
    'XScale',xscale,'YScale',yscale,varargin{:});
hold on                  


end