function[omega]=vorticityfunction(imax,jmax,u,v,dx,dy)
%initiazlize the omega function
omega=zeros(imax,jmax);
for i=2:imax-1
    for j=2:jmax-1
        omega(i,j)=((v(i,j+1)-v(i,j))/dx)-((u(i+1,j)-u(i,j))/dy);
    end
end
%boundary condition for the omega function
omega(1,1)=0.5*(omega(1,2)+omega(2,1));
omega(1,jmax)=0.5*(omega(2,jmax)+omega(1,jmax-1));
omega(imax,1)=0.5*(omega(imax-1,1)+omega(imax,2));
omega(imax,jmax)=0.5*(omega(imax-1,jmax)+omega(imax,jmax-1));

return
end
