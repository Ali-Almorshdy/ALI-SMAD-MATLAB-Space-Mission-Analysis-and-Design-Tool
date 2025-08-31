function setFigureAlwaysOnTop(h)
    % setFigureAlwaysOnTop - Set a MATLAB figure window to always stay on top.
    % 
    % Syntax: setFigureAlwaysOnTop(h)
    % 
    % Input:
    %    h - Handle to the MATLAB figure. If not provided, the current figure (gcf) will be used.
    % 
    % Example:
    %    h = figure;
    %    setFigureAlwaysOnTop(h);

    if nargin < 1
        h = gcf; % Use the current figure if no figure handle is provided
    end
    
    % Ensure the figure is fully rendered
    drawnow;
    
    try
        % Get the JavaFrame (MATLAB R2014b and later)
        jFrame = get(h, 'JavaFrame');
        jWindow = jFrame.getFigurePanelContainer.getTopLevelAncestor;

        % Set the window to always be on top
        jWindow.setAlwaysOnTop(true);
    catch
        % Handle errors and provide feedback
        warning('This method might not work on your version of MATLAB.');
    end
end
