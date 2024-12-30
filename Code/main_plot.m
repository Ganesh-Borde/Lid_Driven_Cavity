tic
clear all
close all

% Reynolds numbers to process
Re_values = [100, 400];

% Initialize variables
u_coll_all = {};
v_coll_all = {};
p_coll_all = {};
omega_all = {};
shi_all = {};
y_all = {};
x_all = {};

for i = 1:length(Re_values)
    Rey = Re_values(i); % Current Reynolds number
    U = 1;
    imax = 129;                             % Number of points in x direction 
    jmax = 129;                             % Number of points in y direction 
    L_x = 1;
    L_y = 1;
    h = L_x / (imax - 1);
    dx = h;
    dy = h;
    x = 0:h:L_x;                            % X domain span
    y = 0:h:L_y;                            % Y domain span
    nuu = 1 / Rey;
    % Relaxation factors
    al_v = 0.8;
    al_p = 0.8;

    % Residual error
    residual_err = 1e-6;

    % Calling the SIMPLE algorithm
    [iterations, u, v, p] = simple(imax, jmax, nuu, h, Rey, residual_err);

    % Transforming to collocated grid
    [u_coll, v_coll, p_coll] = collactedgrid(u, v, p, imax, jmax);

    % Calculating vorticity and stream function
    [omega] = vorticityfunction(imax, jmax, u, v, dx, dy);
    [shi] = streamfunction(imax, jmax, dx, v);

    % Store results for later plotting
    u_coll_all{i} = u_coll;
    v_coll_all{i} = v_coll;
    p_coll_all{i} = p_coll;
    omega_all{i} = omega;
    shi_all{i} = shi;
    y_all{i} = 1 - y; % Reverse y-axis direction for plotting
    x_all{i} = x;
end

% Create mesh grid for contour plots
x_ = ((1:imax)-1) .* h;
y_ = 1 - ((1:jmax)-1) .* h;
[X, Y] = meshgrid(x_, y_);

figureWidth = 800; % Width in pixels
figureHeight = 400; % Height in pixels

%% Plot and Save u-Velocity Contours Side by Side
figure('Name', 'u-Velocity Contours', 'NumberTitle', 'off', 'Position', [100, 100, figureWidth, figureHeight]);
for i = 1:length(Re_values)
    subplot(1, length(Re_values), i);
    contourf(X, Y, u_coll_all{i}, 22, 'LineStyle', 'none');
    colorbar;
    colormap('jet');
    xlabel('$x$', 'FontSize', 12);
    ylabel('$y$', 'FontSize', 12);
    title([' Re = ', num2str(Re_values(i))], 'FontSize', 12);
    axis equal; % Ensures square aspect ratio
end
%saveas(gcf, 'u_velocity_contours.pdf');

%% Plot and Save v-Velocity Contours Side by Side
figure('Name', 'v-Velocity Contours', 'NumberTitle', 'off', 'Position', [100, 100, figureWidth, figureHeight]);
for i = 1:length(Re_values)
    subplot(1, length(Re_values), i);
    contourf(X, Y, v_coll_all{i}, 22, 'LineStyle', 'none');
    colorbar;
    colormap('jet');
    xlabel('$x$', 'FontSize', 12);
    ylabel('$y$', 'FontSize', 12);
    title(['Re = ', num2str(Re_values(i))], 'FontSize', 12);
    axis equal; % Ensures square aspect ratio
end
%saveas(gcf, 'v_velocity_contours.pdf');

%% Plot and Save Vorticity Contours Side by Side
figure('Name', 'Vorticity Contours', 'NumberTitle', 'off', 'Position', [100, 100, figureWidth, figureHeight]);
for i = 1:length(Re_values)
    subplot(1, length(Re_values), i);
    contour(X, Y, omega_all{i}, 200, '--');
    colorbar;
    colormap('jet');
    xlabel('$x$', 'FontSize', 12);
    ylabel('$y$', 'FontSize', 12);
    title(['Re = ', num2str(Re_values(i))], 'FontSize', 12);
    axis equal; % Ensures square aspect ratio
end
%saveas(gcf, 'vorticity_contours.pdf');

%% Plot and Save Stream Function Contours Side by Side
figure('Name', 'Stream Function Contours', 'NumberTitle', 'off', 'Position', [100, 100, figureWidth, figureHeight]);
for i = 1:length(Re_values)
    subplot(1, length(Re_values), i);
    contour(X, Y, shi_all{i}, 1000, '--');
    colorbar;
    colormap('jet');
    xlabel('$x$', 'FontSize', 12);
    ylabel('$y$', 'FontSize', 12);
    title([' Re = ', num2str(Re_values(i))], 'FontSize', 12);
    axis equal; % Ensures square aspect ratio
end
%saveas(gcf, 'stream_function_contours.pdf');

% Plot v-velocity side by side for Re = 100 and Re = 400
figure('Name', 'v-Velocity Comparison', 'NumberTitle', 'off', 'Position', [100, 100, figureWidth, figureHeight]);
for i = 1:length(Re_values)
    subplot(1, length(Re_values), i);
    Rey = Re_values(i);
    [M_x, RE_100_v, RE_400_v] = ghiavalues_v();
    if Rey == 100
        plot(M_x, RE_100_v, 'o', 'DisplayName', 'Ghia');
    elseif Rey == 400
        plot(M_x, RE_400_v, 'o', 'DisplayName', 'Ghia');
    end
    hold on;
    plot(x_all{i}, v_coll_all{i}((imax+1)/2, :), 'LineWidth', 1.5, 'DisplayName', 'Numerical');
    title(['Re = ', num2str(Rey)]);
    xlabel('$x$');
    ylabel('$v$');
    legend('Location', 'southeast');
    hold off;
end


figure('Name', 'u-Velocity Comparison', 'NumberTitle', 'off', 'Position', [100, 100, figureWidth, figureHeight]);
for i = 1:length(Re_values)
    subplot(1, length(Re_values), i);
    Rey = Re_values(i);
    [M_y, RE_100_u, RE_400_u] = ghiavalues_u();
    if Rey == 100
        plot(RE_100_u, M_y, 'o', 'DisplayName', 'Ghia');
    elseif Rey == 400
        plot(RE_400_u, M_y, 'o', 'DisplayName', 'Ghia');
    end
    hold on;
    plot(u_coll_all{i}(:, (jmax+1)/2), y_all{i}, 'LineWidth', 1.5, 'DisplayName', 'Numerical');
    title(['Re = ', num2str(Rey)]);
    xlabel('$u$');
    ylabel('$y$');
    legend('Location', 'southeast');
    hold off;
end

figure('Name', 'streamline', 'NumberTitle', 'off', 'Position', [100, 100, figureWidth, figureHeight]);
for i = 1:length(Re_values)
    subplot(1, length(Re_values), i);
    streamslice(X, Y, u_coll_all{i}, v_coll_all{i}, 4);
    axis([0, 1, 0, 1]);
    colorbar;
    colormap('jet');
    xlabel('$x$', 'FontSize', 12);
    ylabel('$y$', 'FontSize', 12);
    title([' Re = ', num2str(Re_values(i))], 'FontSize', 12);
    axis equal; % Ensures square aspect ratio
end
% figure(8);
% streamslice(X, Y, u_coll, v_coll, 4);
% axis([0, 1, 0, 1]);
% title(['Streamline for Re' num2str(Rey)], 'FontSize', 12);
% xlabel('x/L', 'FontSize', 12);
% ylabel('y/L', 'FontSize', 12);

toc
