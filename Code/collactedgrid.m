function[u_coll,v_coll,p_coll]=collactedgrid(u,v,p,imax,jmax)
%Final collocated variables
u_coll(imax,jmax)=0;
v_coll(imax,jmax)=0;
p_coll(imax,jmax)=1;
u_coll(1,:) = 1;

% After the converged solution, we map the staggered variables to
% collocated variables
for i = 1:imax
    for j = 1:jmax
        u_coll(i,j) = 0.5*(u(i,j) + u(i+1,j));
        v_coll(i,j) = 0.5*(v(i,j) + v(i,j+1));
        p_coll(i,j) = 0.25*(p(i,j) + p(i,j+1) + p(i+1,j) + p(i+1,j+1));
    end
end
return
end