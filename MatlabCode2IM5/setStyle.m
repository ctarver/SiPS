function [] = setStyle( fh, style, bScaleLine )
%SETSTYLE Changes a figure to match a certain style
%   Inputs:
%       fh - figure handle, e.g. gcf for current figure
%       style - path to style file
%       bScaleLine - 0: don't scale line width of plots, 1: scale line
%           width of plots (default)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Authors:      Janis Werner, Jaakko Marttila
% Organization: Tampere University of Technology (TUT)
% Contact:      janis.werner@tut.fi, jaakko.marttila@tut.fi
% Year:         2012
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin == 2
    bScaleLine = 1;
end

%Undock figure in order to make position setting possible
undock(fh);

%Load style
eval(['run ''' style '''']);

%Axis handle
ah = findobj(fh,'type','axes','-not','Tag','legend','-not','Tag','Colorbar');

%Legend handle
lh =  findobj(fh,'Type','axes','Tag','legend');

%Xlabel handle
xh = get(ah, 'xlabel');

%Ylabel handle
yh = get(ah, 'ylabel');

%Set units
set(fh, 'Units', unitFigure);
set(ah, 'Units', unitAxis);

%Get figure position
figPos = get(fh, 'Position');
figX = figPos(1);
figY = figPos(2);

%Scale all variables
sFigureWidth = scaleFactor*figureWidth;
sFigureHeight = scaleFactor*figureHeight;
sAxisX = scaleFactor*axisX;
sAxisY = scaleFactor*axisY;
sAxisWidth = scaleFactor*axisWidth;
sAxisHeight = scaleFactor*axisHeight;
sAxisFontSize = scaleFactor*axisFontSize;
sLabelFontSize = scaleFactor*labelFontSize;

%Change pixels per inch on screen if desired
if exist('screenPixelsPerInch', 'var')
    sScreenPixelsPerInch = screenPixelsPerInch/scaleFactor;
    set(0, 'ScreenPixelsPerInch', sScreenPixelsPerInch)
end

%Adjust size and position
set(fh, 'Position', [figX, figY, sFigureWidth, sFigureHeight]);
set(ah, 'Position', [sAxisX, sAxisY, sAxisWidth, sAxisHeight]);

%Set text interpreter
set(lh, 'Interpreter', legendInterpreter);
set(xh, 'Interpreter', labelInterpreter);
set(yh, 'Interpreter', labelInterpreter);

%Set the font for the x and y tick
set(ah, 'FontName', tickFont);

%Adjust text size
set(ah, 'FontSize', sAxisFontSize);
set(xh, 'FontSize', sLabelFontSize);
set(yh, 'FontSize', sLabelFontSize);

%Scale line width of the axis
% aLineWidth = get(ah, 'LineWidth');
aLineWidth = 0.5;   %standard
set(ah, 'LineWidth', scaleFactor*aLineWidth);

if bScaleLine
    %Apply the scaling to the linewidth of the plots
    ih = findobj(ah, 'Type', 'line');
    for ii = 1:length(ih)
        lineWidth = get(ih(ii), 'LineWidth');
        set(ih(ii), 'LineWidth', scaleFactor*lineWidth);
    end
end

box(ah, 'on');

end

