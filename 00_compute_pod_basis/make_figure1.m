%% FEM APPROXIMATION VISUALIZATION (t0, tp/3, 2tp/3)
% Brusselator reaction-diffusion equations witgh parameters :
% alpha = 1, beta = 3, nu = 0.002
% This code generates PDF plots for FEM approximation at three different
% times for components (u) and (v) of the system.

%% 1. Load the data.
% Tri    : connectivity matrix for quadratic elements (80 subdivisions)
% z      : nodes (x,y)-coordinates
% tiempos: time vector (1x2049)
% Zh     : coupled FEM solution vector [u; v]
load snapshots_alpha1_beta3_nu0p002.mat

% Define mesh variables for plotting
TriLin = Tri(:,1:3);  % Extract linear connectivity
x = z(:,1);               % x-coordinates
y = z(:,2);               % y-coordinates
nn = length(z);           % Number of nodes
nTimeSteps = length(tiempos)-1;

%% Compute time at maximum H1 norm
NormZhH1 = sqrt((abs(sum(Zh .* (S2 * Zh), 1))));
% S2 = kron(speye(2),Sh)
[maxNormZhH1, idxNormZh1] = max(NormZhH1);
timeMaxNormZhH1 = tiempos(idxNormZh1);

%% 2. Take FEM solution vector at selected times
% Select time indices for visualization: initial, 1/3 period, and 2/3 period
Indices = [floor([1, nTimeSteps/3, nTimeSteps*2/3]), idxNormZh1];
nIndices = length(Indices); 

Uh = zeros(nn, nIndices); 
Vh = zeros(nn, nIndices);

for k = 1:nIndices
    idx = Indices(k); % Take time corresponding time index
    Uh(:,k) = Zh(1:nn, idx);     
    Vh(:,k) = Zh(nn+1:end, idx); 
end

%% 3. Define plotting variables
fontSizeLabel = 20;      
fontSizeAxis  = 18;      
fontSizeTitle = 20;      

% Vertical offset for Z-axis label to improve visibility
zPos = [1.2236, 1.0484, 2]; 

%% 4. Plotting component u
for k = 1:nIndices
    figU = figure('Visible', 'off'); % Figure Window does not appear
    ax = gca;
    
    trimesh(TriLin, x, y, Uh(:,k));
    view(244, 12);
    axis square;  
    
    % Axis style
    ax.FontSize = fontSizeAxis;
    ax.TickLabelInterpreter = 'latex';
    
    xlabel('$x$', 'Interpreter','latex','FontSize',fontSizeLabel)
    ylabel('$y$', 'Interpreter','latex','FontSize',fontSizeLabel)
    zlabel('$u$', 'Interpreter','latex','FontSize',fontSizeLabel, ...
          'Position', zPos + [0.2, 0, 0], 'Rotation', 0)

    xticks([0, 0.5, 1]); yticks([0, 0.5, 1]);
    zlim([0 4]); zticks(0:4);

    title(['$t = ', num2str(tiempos(Indices(k))), ' \ s $'], ...
          'Interpreter','latex', 'FontSize',fontSizeTitle)

    % Export to PDF
    exportgraphics(figU, sprintf('fig_u_t%d.pdf', k), ...
                   'ContentType', 'vector', 'BackgroundColor', 'none');
    close(figU); 
    % Content type 'vector' saves the image as a vector/geometry way to
    % pixels, make transparent backgrounds
end

%% 5. Plotting component v
for k = 1:nIndices
    figV = figure('Visible', 'off'); % Figure Window does not appear
    ax = gca;
    
    trimesh(TriLin, x, y, Vh(:,k));
    view(244, 12);
    axis square;  
    
    % Axis style
    ax.FontSize = fontSizeAxis;
    ax.TickLabelInterpreter = 'latex';
    
    xlabel('$x$', 'Interpreter','latex','FontSize',fontSizeLabel)
    ylabel('$y$', 'Interpreter','latex','FontSize',fontSizeLabel)
    zlabel('$v$', 'Interpreter','latex','FontSize',fontSizeLabel, ...
           'Position', zPos + [0.2, 0, 0], 'Rotation', 0)

    xticks([0, 0.5, 1]); yticks([0, 0.5, 1]);
    zlim([0 5]); zticks(0:5);

    title(['$t = ', num2str(tiempos(Indices(k))), ' \ s $'], ...
          'Interpreter','latex', 'FontSize',fontSizeTitle)

    % Export to PDF
    exportgraphics(figV, sprintf('fig_v_t%d.pdf', k), ...
                   'ContentType', 'vector', 'BackgroundColor', 'none');
    close(figV); 
    % Content type 'vector' saves the image as a vector/geometry way instead of
    % pixels
end