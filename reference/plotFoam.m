%% Plotting openFOAM fields in matlab

clear all;
close all;

% Load the data
path = './stokesFoam/postProcessing/surfaces/';
num = dir(path);
num = num(3).name; % Fetch the last iteration
p_path = [path num '/p_normal.raw'];
u_path = [path num '/U_normal.raw'];
p_data = readmatrix(p_path, "FileType","text");
u_data = readmatrix(u_path, "FileType","text");

Nx = length(unique(p_data(:,1)));
Ny = length(unique(p_data(:,2)));

U = zeros([Nx, Ny]);
V = zeros([Nx, Ny]);
p = zeros([Nx, Ny]);
[X, Y] = meshgrid(unique(p_data(:,1)), unique(p_data(:,2)));
for i = 1:Nx
    for j = 1:Ny
        idx = p_data(:,1) == X(i,j) & p_data(:,2) == Y(i,j);
        U(i,j) = u_data(idx, 4);
        V(i,j) = u_data(idx, 5);
        p(i,j) = p_data(idx, 4);
    end
end

% Make pressure zero mean
p = p - mean(p,"all");

% Plot p
figure()
surface(X,Y,p, 'EdgeColor','none');
colorbar()
xlabel('x');
ylabel('y');
saveas(gcf, 'p_field.png');

% Plot U
c = 3;      % Plot every third value to make the plot readable
figure()
hold on
quiver(X(1:c:Nx,1:c:Nx),Y(1:c:Nx,1:c:Nx),U(1:c:Nx,1:c:Nx),V(1:c:Nx,1:c:Nx), 2, 'LineWidth',1.05)
xlabel('x');
ylabel('y');
axis tight
%contour(X,Y,sqrt(U.^2+V.^2))
saveas(gcf, 'U_field.png');

% Plot streamlines (assumes [0,1]^2 domain)
d = 0.2;    % Spacing for diagonal lines
dd = 0.01;  % Spacing for corner lines
figure()
startX = 0:d:1;
startY = 0:d:1;
startX = [startX, startX, 0.5*ones(size(startX)), 0:dd:0.2, 0.8:dd:1];
startY = [startY, flip(startY), startY, 0:dd:0.2, flip(0:dd:0.2)];
streamline(X,Y,U,V,startX,startY, [0.4, 1e4]) 
xlabel('x');
ylabel('y');
saveas(gcf, 'streamlines.png');