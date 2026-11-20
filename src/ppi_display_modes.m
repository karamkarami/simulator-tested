function [modes, names] = ppi_display_modes(sys)
%PPI_DISPLAY_MODES  Display modes of the PPI and their names in the panel.
%   [modes, names] = PPI_DISPLAY_MODES(sys): 'plots', 'detections', 'both', 'video', 'dbf' (only
%   with the DBF: sys.dbf.enabled), 'slc' (only with the sidelobe canceller: sys.slc.enable) and
%   'mf'. The C receiver lists the same (ppi_qt.cpp).
modes = {'plots', 'detections', 'both', 'video'};
names = {'Plots (after extractor)', 'Detections (before extractor)', 'Plots and detections', ...
         'Video (DDC, elements averaged)'};
if sys.dbf.enabled
    modes{end+1} = 'dbf';  names{end+1} = 'DBF out (selected beam)';
end
if isfield(sys, 'slc') && sys.slc.enable
    modes{end+1} = 'slc';  names{end+1} = 'SLC out (selected beam)';
end
modes{end+1} = 'mf';  names{end+1} = 'MF out (selected beam)';
end
