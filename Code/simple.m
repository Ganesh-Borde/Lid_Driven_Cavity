function[iterations,u,v,p]=simple(imax,jmax,nuu,h,Rey,residual_err)
%Staggered variables
%u,v=the initial guess
%star variables are the intermediate results found from the momentum
%equation
%final varibales are for the final values computed after correction
%coll variables are for the values of variable form stagggered grid to
%collacted grid
%correction variable for the variable at corrected step
%intizaling the u velocity variables for u u_star u_final u_coll

u=zeros(imax+1,jmax);
u_star=zeros(imax+1,jmax);
u(1,:)=2;
%intizaling the v velocity variables for v v_star v_final v_coll
v=zeros(imax,jmax+1);
v_star=zeros(imax,jmax+1);

%intizaling the p velocity variables for p   p_final p_correction
p=ones(imax+1,jmax+1);
%p_final=ones(imax+1,jmax+1);
p_corr=zeros(imax+1,jmax+1);

%the coefficient for the pressure term
b=ones(imax+1,jmax+1);
D_h=zeros(imax+1,jmax+1);
D_v=zeros(imax+1,jmax+1);

%the continuity for star velcoiteis
continuity_star=zeros(imax+1,jmax+1);

% Under-relaxation factors
al_v = 0.8;
al_p = 0.8;

a=nuu;
err = 1;
iterations = 0;

%%using the while loop for the iterative method
while err > residual_err
    %%solving the u_star for finding the pressure Poisson equation
    for i = 2:imax
        for j = 2:jmax - 1
            u_Ea = 0.5*(u(i,j) + u(i,j+1)); a_E = -0.5*u_Ea*h + a;
            u_We = 0.5*(u(i,j) + u(i,j-1)); a_W = 0.5*u_We*h + a;
            v_No = 0.5*(v(i-1,j) + v(i-1,j+1));a_N = -0.5*v_No*h + a;
            v_So = 0.5*(v(i,j) + v(i,j+1)); a_S = 0.5*v_So*h + a;      
            a_e = 0.5*u_Ea*h - 0.5*u_We*h + 0.5*v_No*h - 0.5*v_So*h + 4*a;
            A_e = -h;
            D_h(i,j) = A_e/a_e;
        
       %%calculating the u_star            
            u_star(i,j) = (a_E*u(i,j+1) + a_W*u(i,j-1) + a_N*u(i-1,j) + a_S*u(i+1,j))/a_e + D_h(i,j)*(p(i,j+1) - p(i,j));
            
        end
    end    
    %Boundary condition for u-momentum i.e u_star terms
    u_star(1,:) = 2 - u_star(2,:);
    u_star(imax + 1,:) = -u_star(imax,:);
    u_star(2:imax,1) = 0;
    u_star(2:imax,jmax) = 0;
    
    %solving the v_star for finding the pressure Poisson equation 
    for i = 2:imax - 1
        for j = 2:jmax
            u_Ea = 0.5*(u(i,j) + u(i+1,j));a_E = -0.5*u_Ea*h + a;
            u_We = 0.5*(u(i,j-1) + u(i+1,j-1));a_W = 0.5*u_We*h + a;
            v_No = 0.5*(v(i-1,j) + v(i,j));a_N = -0.5*v_No*h + a;
            v_So = 0.5*(v(i,j) + v(i+1,j));a_S = 0.5*v_So*h + a;                   
            a_n = 0.5*u_Ea*h - 0.5*u_We*h + 0.5*v_No*h - 0.5*v_So*h + 4*a;
            A_n = -h;
            D_v(i,j) = A_n/a_n;
            
            %D_v(i,j) = A_n/a_n;            
            v_star(i,j) = (a_E*v(i,j+1) + a_W*v(i,j-1) + a_N*v(i-1,j) + a_S*v(i+1,j))/a_n + D_v(i,j)*(p(i,j) - p(i+1,j));
        end
    end    
    %Boundary condition for v-momentum i.e v_star terms
    v_star(:,1) =-v_star(:,2);
    v_star(:,jmax + 1) =-v_star(:,jmax);
    v_star(1,2:jmax) = 0;
    v_star(imax,2:jmax) = 0;    
    %intizaling the pressure correction to zero
    p_corr(1:imax+1,1:jmax+1)=0;
    %solving the pressure poisson equation
    for i = 2:imax
        for j = 2:jmax
            %initiating the coefficient for pressure terms
            a_Ea = -D_h(i,j)*h;a_We = -D_h(i,j-1)*h;
            a_No = -D_v(i-1,j)*h;a_So = -D_v(i,j)*h;
            a_P = a_Ea+a_We+a_No+a_So;
            %the RHS of pressure poisson equation
            b(i,j) = -(u_star(i,j)-u_star(i,j-1))*h +(v_star(i,j)-v_star(i-1,j))*h;
            continuity_star(i,j)=b(i,j);
            %the pressure Poisson equation            
            p_corr(i,j) = (a_Ea*p_corr(i,j+1) + a_We*p_corr(i,j-1) + a_No*p_corr(i-1,j) + a_So*p_corr(i+1,j) + b(i,j))/a_P;
            %calculting the Final pressure by adding the pressure correction
            
            p(i,j) = p(i,j) + al_p*p_corr(i,j);

        end
    end    
    
    %Boundary condition i.e neumann type
    p(1,:) = p(2,:);
    p(imax + 1,:) = p(jmax,:);
    p(:,1) = p(:,2);
    p(:,jmax + 1) = p(:,jmax);
    
    %Updating the  u -velocities from calculated correction
    for i = 2:imax
        for j = 2:jmax - 1
            %correction term
            u(i,j)=D_h(i,j)*(p_corr(i,j+1) - p_corr(i,j));
            u(i,j) = u_star(i,j) + al_v*u(i,j);
            
        end
    end    
    

    %Updating the  v-velocities from calculated correction
    
    for i = 2:imax - 1
        for j = 2:jmax
            %correction term
            v(i,j)=D_v(i,j)*(p_corr(i,j) - p_corr(i+1,j));
            %final term
            v(i,j) = v_star(i,j) + al_v*v(i,j);        
        end
    end
    %Boundary condition for u-momentum i.e u terms
    u(1,:) = 2 - u(2,:);
    u(imax + 1,:) = -u(imax,:);
    u(2:imax,1) = 0;
    u(2:imax,jmax) = 0;

    %Boundary condition for v-momentum i.e v_star terms
    v(:,1) = -v(:,2);
    v(:,jmax + 1) = -v(:,jmax);
    v(1,2:jmax) = 0;
    v(imax,2:jmax) = 0;  
    
    
    %the residual calcualtions form the RHS of pressure poisson equation
    %which is the 
    err = 0;
    for i = 2:imax
        for j = 2:jmax
            err = err + abs(continuity_star(i,j));
        end
    end     
    
    %calcuating and monitoring the iterations vs residuals for the compuations
    if(rem(iterations, 200)) == 0  
       figure(1);
       semilogy(iterations, err,'ok')  %%using the semiology function which plots 
       %the iterations vs residuals for every step with saving in memory for new
       %variables which save the computation time ,since we don’t require the error 
       %value at each iterations
       hold on
       
       xlabel('Iterations')
       ylabel('Residual Error')
       if Rey==100
         title('Residual plot for Re 100');
       else
         title('Residual plot for Re 400');
       end
       
       
    end
    
    %intializing the initial variable for next iterations with the final %%variables found in the loop
 
   

    
    %continuing the iterations by adding it up
    iterations = iterations + 1;
    
end

end



