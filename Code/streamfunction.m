function[shi]=streamfunction(imax,jmax,dx,u)
shi=zeros(imax,jmax);
for i=2:jmax-1
    shi(2:imax-1,i)=shi(2:imax-1,i-1)+dx*u(2:imax-1,i);
end

return
end
