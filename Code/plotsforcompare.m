function []=plotsforcompare(u_coll,v_coll,y,x,imax,jmax) 
%%plotting the u and v velocity with ghai paper
figure(11);
plot(u_coll(:,(jmax+1)/2),1-y, 'LineWidth', 1)


xlabel('u')
ylabel('y')
legend('Numerical', 'Benchmark', 'location', 'southeast')
figure(12)
plot(x,v_coll((imax+1)/2,:),'LineWidth',1)
xlabel('x')
ylabel('v')
legend('Numerical', 'Benchmark', 'location', 'southeast')