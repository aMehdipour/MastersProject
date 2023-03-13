function save_figures(filepath)
% Save all open figures as .fig and .png files to the specified filepath
%
% Input:
% - filepath: a string representing the file path where the figures will be saved
%
% Example usage:
% save_figures('/home/username/my_figures')

% Convert to Linux-style file path
filepath = strrep(filepath, '\', '/');

% Create directory if it doesn't exist
if ~exist(filepath, 'dir')
    mkdir(filepath);
end

% Get handles to all open figures
fig_handles = findobj('Type', 'figure');

% Loop through all open figures
for i = 1:length(fig_handles)
    % Get the figure name
    fig_name = get(fig_handles(i), 'Name');
    
    % Set the figure as the current figure
    set(groot, 'CurrentFigure', fig_handles(i));
    
    % Construct the file names
    fig_filename = sprintf('%s.fig', fig_name);
    png_filename = sprintf('%s.png', fig_name);
    
    % Construct the file paths
    fig_filepath = fullfile(filepath, fig_filename);
    png_filepath = fullfile(filepath, png_filename);
    
    % Save the figure as a .fig file
    savefig(fig_filepath);
    
    % Save the figure as a .png file
    saveas(gcf, png_filepath, 'png');
end
