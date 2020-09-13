function [fh, y_fine, y_parareal, y_exact] = ExpGridPlot(params)
%EXPGRIDPLOT Summary of this function goes here
%   Detailed explanation goes here
[ y_fine, y_parareal, y_exact, se_parareal ] = ExpGridPlotData(params);
[ fh ] = ExpGridPlotPlotter(params, y_fine, y_parareal, y_exact, se_parareal);
end