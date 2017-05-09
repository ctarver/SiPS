%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Stylefile for IEEE paper with two columns, to be used with 'setStyle'.
% For square figures.
%
% Authors:      Janis Werner, Jaakko Marttila
% Organization: Tampere University of Technology (TUT)
% Contact:      janis.werner@tut.fi, jaakko.marttila@tut.fi
% Year:         2012
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Control how many pixels are displayed per inch
% - Comment if no change desired
% - Affects only how the figure is displayed
% screenPixelsPerInch = 120;

%Unfortuantely Matlab does not scale plots properly, e.g. the legends
%become to big for small figures. This factor can be used to to create a
%scaleFactor times bigger figure with respective properties. In latex
%figures should be important with a scaling 1/scaleFactor.
scaleFactor = 2;

%Unit used for figure
unitFigure = 'inches'; 

%Unit used for axis
unitAxis = 'inches';

%Width and height of figure
figureWidth = 3.5;
figureHeight = 1.5;

%Position of axis
axisX = 0.5; 
axisY = 0.4;

%Width of axis
axisWidth = 2.7;
axisHeight = 1;

%Font size of x and y labels
labelFontSize = 9;

%Interpreter of x and y labels 
labelInterpreter = 'latex';

%Font size of axis and legend
axisFontSize = 7;

%Interpreter of legend 
legendInterpreter = 'latex';

%Tick font (It is really hard to set the tick interpreter to latex)
tickFont = 'Arial';

%Interpreter for ticks (at the moment only 'latex' results in changes)
% tickInterpreter = 'latex';
