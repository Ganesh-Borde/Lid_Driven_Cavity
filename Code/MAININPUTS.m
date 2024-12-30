tic
clear all
close all


%%Defining the problem domain
Rey = 1000; %Reynolds number
U=1;
imax= 513;                             %%no.of points in x direction
jmax =513;                               % no_of_points in y direction
L_x= 1;
L_y=1;
h = L_x/(imax-1);
dx=h;
dy=h;
x = 0:h:L_x;                                 %X domain span
y = 0:h:L_y;                                  %Y domain span
nuu = 1/Rey;
% Relaxation factors
al_v = 0.8;
al_p = 0.8;


%Residual error
residual_err= 1e-6;

%calling the whole simple algorithm function
[iterations,u,v,p]=simple(imax,jmax,nuu,h,Rey,residual_err);

%calling the collacted grid values from the staggered grid values
[u_coll,v_coll,p_coll]=collactedgrid(u,v,p,imax,jmax);

%calling the vorticity function for the vorticity value found 
[omega]=vorticityfunction(imax,jmax,u,v,dx,dy);

%calling the stream function value
[shi]=streamfunction(imax,jmax,dx,v);

%Centreline u velocity 
%plotting the u and v velocity with ghai paper,which data has been stored in a function
%calling the ghia data for u velocity

[M_y, RE_100, RE_400] = ghiavalues_u();
figure(2);

if Rey == 100
    plot(RE_100, M_y, 'o');
    hold on;
elseif Rey == 400
    plot(RE_400, M_y, 'o');
    hold on;
else
    % Skip comparison if Re is not 100 or 400
    disp('No comparison data available for this Reynolds number.');
end

%% Plotting the u velocity along the centerline
plot(u_coll(:, (jmax+1)/2), 1-y, 'LineWidth', 1);
xlabel('u');
ylabel('y/L');
legendEntries = {'Numerical'};
if Rey == 100 || Rey == 400
    legendEntries = [{'Ghia'}, legendEntries];
end
legend(legendEntries, 'Location', 'southeast');

%% Plotting v velocity along the x-direction
figure(3);
[M_x, RE_100, RE_400] = ghiavalues_v();
if Rey == 100
    plot(M_x, RE_100, 'o');
    hold on;
elseif Rey == 400
    plot(M_x, RE_400, 'o');
    hold on;
else
    disp('No comparison data available for this Reynolds number.');
end

plot(x, v_coll((imax+1)/2, :), 'LineWidth', 1);
xlabel('x/L');
ylabel('v');
legendEntries = {'Numerical'};
if Rey == 100 || Rey == 400
    legendEntries = [{'Ghia'}, legendEntries];
end
legend(legendEntries);

if Rey == 100
    title('v velocity along the x at y=0.5 RE100');
elseif Rey == 400
    title('v velocity along the x at y=0.5 RE400');
else
    title(['v velocity along the x at y=0.5 RE' num2str(Rey)]);
end
hold off;

%% Creating mesh grid for plotting the contour
x_ = ((1:imax)-1) .* h;
y_ = 1 - ((1:jmax)-1) .* h;
[X, Y] = meshgrid(x_, y_);

%% Plotting the u velocity contour
figure(4);
contourf(X, Y, u_coll, 22, 'LineStyle', 'none');
colorbar;
colormap('jet');
xlabel('x/L', 'FontSize', 12);
ylabel('y/L', 'FontSize', 12);
title(['u velocity contour Re' num2str(Rey)], 'FontSize', 12);

%% Plotting the v velocity contour
figure(5);
contourf(X, Y, v_coll, 22, 'LineStyle', 'none');
colorbar;
colormap('jet');
xlabel('x/L', 'FontSize', 12);
ylabel('y/L', 'FontSize', 12);
title(['v velocity contour Re' num2str(Rey)], 'FontSize', 12);

%% Plotting the vorticity contour
figure(6);
contour(X, Y, omega, 200, '--');
colorbar;
colormap('jet');
xlabel('x/L', 'FontSize', 12);
ylabel('y/L', 'FontSize', 12);
title(['vorticity contour Re' num2str(Rey)], 'FontSize', 12);

%% Plotting the stream function contour
figure(7);
contour(X, Y, shi, 1000, '--');
colorbar;
colormap('jet');
xlabel('x/L', 'FontSize', 12);
ylabel('y/L', 'FontSize', 12);
title(['stream contour Re' num2str(Rey)], 'FontSize', 12);

%% Printing the total number of iterations
fprintf('The total number of iterations: %d\n', iterations);

%% Printing CPU time
fprintf('The total time taken: %d seconds\n', cputime);

%% Using the streamline function
figure(8);
streamslice(X, Y, u_coll, v_coll, 4);
axis([0, 1, 0, 1]);
title(['Streamline for Re' num2str(Rey)], 'FontSize', 12);
xlabel('x/L', 'FontSize', 12);
ylabel('y/L', 'FontSize', 12);

%%using tic and tac for execution time
toc








