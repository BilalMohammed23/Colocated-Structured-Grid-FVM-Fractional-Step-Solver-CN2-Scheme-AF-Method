clear
close all
clc

Re = 7500;
gamma = 1/Re;

Nx = 128;
x = linspace(0,1,Nx);
del_x = x(2) - x(1);
y = linspace(0,1,Nx);
del_y = y(2) - y(1);

u_cell_cntr = 0*ones(Nx+2,Nx+2);
v_cell_cntr = 0*ones(Nx+2,Nx+2);

u_cell_face_sides_intp = 0*ones(Nx,Nx+3);
v_cell_face_sides_intp = 0*ones(Nx,Nx+3);

u_cell_face_tb_intp = 0*ones(Nx+3,Nx);
v_cell_face_tb_intp = 0*ones(Nx+3,Nx);

P = 0*ones(Nx+2,Nx+2);

%% B.C
[u_cell_cntr,v_cell_cntr,u_cell_face_sides_intp,v_cell_face_sides_intp,u_cell_face_tb_intp,v_cell_face_tb_intp] = bound_cond(u_cell_cntr,v_cell_cntr,u_cell_face_sides_intp,v_cell_face_sides_intp,u_cell_face_tb_intp,v_cell_face_tb_intp);

%% interpolated values
[u_cell_face_sides_intp,v_cell_face_sides_intp,u_cell_face_tb_intp,v_cell_face_tb_intp] = interpl_faces(Nx,u_cell_cntr,v_cell_cntr,u_cell_face_sides_intp,v_cell_face_sides_intp,u_cell_face_tb_intp,v_cell_face_tb_intp);

%% Time-advancement

T = 50;
t = 0;
count = 0;

while t<T
    %% Step-1:
    % tridiag system
    
    delta_t_max = dt_lid_driven_cavity(u_cell_cntr,v_cell_cntr,del_x,del_y,Nx,gamma);

    if t + delta_t_max > T
        delta_t_max = T - t;
    end

    % Have used Approximate-factorization method instead of Direct Solve -
    % Reduces to 4 tridiagonal systems
    
    [del_ustar_bar] = delustarbar_fnc_tridiag(Nx,gamma,del_x,del_y,delta_t_max,u_cell_face_sides_intp,v_cell_face_tb_intp,u_cell_cntr,P);
    
    [del_vstar_bar] = delvstarbar_fnc_tridiag(Nx,gamma,del_x,del_y,delta_t_max,u_cell_face_sides_intp,v_cell_face_tb_intp,v_cell_cntr,P,del_ustar_bar);
    
    [del_vstar] = delvstar_fnc_tridiag(Nx,gamma,del_y,delta_t_max,v_cell_face_tb_intp,v_cell_cntr,del_vstar_bar);
    
    [del_ustar] = delustar_fnc_tridiag(Nx,gamma,del_y,delta_t_max,v_cell_face_tb_intp,u_cell_cntr,del_vstar,del_ustar_bar);
    
    %% 
    
    % yield predicted cell center velocities and face normal velocities
    
    u_cell_cntr_pred1 = u_cell_cntr + del_ustar;
    
    v_cell_cntr_pred1 = v_cell_cntr + del_vstar;
    
    % u_v_face_predicated1 % actually the pred1_faces not required, while we
    % only need pred1_cell_cntr b.c to satisfy because it moves to pred2
    % but we use a common b.c function so get 
    
    u_cell_face_sides_intp_pred1 = 0*ones(Nx,Nx+3);
    v_cell_face_sides_intp_pred1 = 0*ones(Nx,Nx+3);
    
    u_cell_face_tb_intp_pred1 = 0*ones(Nx+3,Nx);
    v_cell_face_tb_intp_pred1 = 0*ones(Nx+3,Nx);
    
    % B.C
    [u_cell_cntr_pred1,v_cell_cntr_pred1,u_cell_face_sides_intp_pred1,v_cell_face_sides_intp_pred1,u_cell_face_tb_intp_pred1,v_cell_face_tb_intp_pred1] = bound_cond(u_cell_cntr_pred1,v_cell_cntr_pred1,u_cell_face_sides_intp_pred1,v_cell_face_sides_intp_pred1,u_cell_face_tb_intp_pred1,v_cell_face_tb_intp_pred1);
    
    % interpolated values of predicted velocities
    [u_cell_face_sides_intp_pred1,v_cell_face_sides_intp_pred1,u_cell_face_tb_intp_pred1,v_cell_face_tb_intp_pred1] = interpl_faces(Nx,u_cell_cntr_pred1,v_cell_cntr_pred1,u_cell_face_sides_intp_pred1,v_cell_face_sides_intp_pred1,u_cell_face_tb_intp_pred1,v_cell_face_tb_intp_pred1);
    
    %% Step-2:
    
    % finding u**
    
    u_cell_cntr_pred2 = u_cell_cntr_pred1;
    v_cell_cntr_pred2 = v_cell_cntr_pred1;
    
    for k = 2:Nx+1
        for l = 2:Nx+1
    
            dq1dx = delta_t_max * (P(k,l+1) - P(k,l-1))/(2*del_x);
            dq1dy = delta_t_max * (P(k-1,l) - P(k+1,l))/(2*del_y);
    
            u_cell_cntr_pred2(k,l) = u_cell_cntr_pred1(k,l) + ((1/2) * dq1dx);
            v_cell_cntr_pred2(k,l) = v_cell_cntr_pred1(k,l) + ((1/2) * dq1dy);
    
        end
    end
    
    %% Step-3:
    
    % finding u(n+1)
    
    % u_v_face_predicated2
    
    u_cell_face_sides_intp_pred2 = 0*ones(Nx,Nx+3);
    v_cell_face_sides_intp_pred2 = 0*ones(Nx,Nx+3);
    
    u_cell_face_tb_intp_pred2 = 0*ones(Nx+3,Nx);
    v_cell_face_tb_intp_pred2 = 0*ones(Nx+3,Nx);
    
    % B.C
    [u_cell_cntr_pred2,v_cell_cntr_pred2,u_cell_face_sides_intp_pred2,v_cell_face_sides_intp_pred2,u_cell_face_tb_intp_pred2,v_cell_face_tb_intp_pred2] = bound_cond(u_cell_cntr_pred2,v_cell_cntr_pred2,u_cell_face_sides_intp_pred2,v_cell_face_sides_intp_pred2,u_cell_face_tb_intp_pred2,v_cell_face_tb_intp_pred2);
    
    % interpolated values of predicted velocities
    [u_cell_face_sides_intp_pred2,v_cell_face_sides_intp_pred2,u_cell_face_tb_intp_pred2,v_cell_face_tb_intp_pred2] = interpl_faces(Nx,u_cell_cntr_pred2,v_cell_cntr_pred2,u_cell_face_sides_intp_pred2,v_cell_face_sides_intp_pred2,u_cell_face_tb_intp_pred2,v_cell_face_tb_intp_pred2);
    
    % Pressure-Poisson Equation
    RHS_P_pred2 = zeros(Nx,Nx);
    
    for i = 1 : Nx
        for j = 1 : Nx
            P_R1 = (u_cell_face_sides_intp_pred2(i,j+2) - u_cell_face_sides_intp_pred2(i,j+1)) * del_y;
            P_R2 = (v_cell_face_tb_intp_pred2(i+1,j) - v_cell_face_tb_intp_pred2(i+2,j)) * del_x;
            RHS_P_pred2(i,j) = P_R1 + P_R2;
        end
    end
    
    delt_P_cell_cntr = pressure_poisson_GS(Nx,del_x,del_y,RHS_P_pred2);
    
    [delt_P_sides,delt_P_tb] = delt_P_faces(Nx,del_x,del_y,delt_P_cell_cntr);
    
    u_cell_face_sides_corr = u_cell_face_sides_intp_pred2;
    v_cell_face_tb_corr    = v_cell_face_tb_intp_pred2;
    
    
    
    % Correct normal face velocities
    u_cell_face_sides_corr(:,2:Nx+2) = u_cell_face_sides_intp_pred2(:,2:Nx+2) - (1/2)*delt_P_sides;
    
    v_cell_face_tb_corr(2:Nx+2,:) = v_cell_face_tb_intp_pred2(2:Nx+2,:) - (1/2)*delt_P_tb;
    
    % B.C
    u_cell_face_sides_corr(:,2)    = 0;  
    u_cell_face_sides_corr(:,Nx+2) = 0;  
    
    v_cell_face_tb_corr(2,:)       = 0;  
    v_cell_face_tb_corr(Nx+2,:)    = 0;  
    
    % Correct cell center velocities
    u_cell_cntr_corr = u_cell_cntr_pred2;
    v_cell_cntr_corr = v_cell_cntr_pred2;
    
    for k = 2:Nx+1
        for l = 2:Nx+1
    
            dqdx = (delt_P_cell_cntr(k,l+1) - delt_P_cell_cntr(k,l-1))/(2*del_x);
            dqdy = (delt_P_cell_cntr(k-1,l) - delt_P_cell_cntr(k+1,l))/(2*del_y);
    
            u_cell_cntr_corr(k,l) = u_cell_cntr_pred2(k,l) - (1/2)*dqdx;
            v_cell_cntr_corr(k,l) = v_cell_cntr_pred2(k,l) - (1/2)*dqdy;
    
        end
    end
    
    P_new = P + delt_P_cell_cntr / delta_t_max;
    
    % B.C
    [u_cell_cntr_corr,v_cell_cntr_corr,u_cell_face_sides_corr,v_cell_face_sides_intp,u_cell_face_tb_intp,v_cell_face_tb_corr] = bound_cond(u_cell_cntr_corr,v_cell_cntr_corr,u_cell_face_sides_corr,v_cell_face_sides_intp,u_cell_face_tb_intp,v_cell_face_tb_corr);
    
    u_cell_cntr = u_cell_cntr_corr;
    v_cell_cntr = v_cell_cntr_corr;
    
    u_cell_face_sides_intp = u_cell_face_sides_corr;
    v_cell_face_tb_intp    = v_cell_face_tb_corr;
    
    P = P_new;

    t = t + delta_t_max;
    count = count + 1;
    fprintf('%4d: t = %.4f, dt = %.3e\n',count,t,delta_t_max)

end

%% Lid-driven cavity contour

u_plot = u_cell_cntr(2:Nx+1,2:Nx+1);
v_plot = v_cell_cntr(2:Nx+1,2:Nx+1);

u_plot = flipud(u_plot);
v_plot = flipud(v_plot);

vel_mag = sqrt(u_plot.^2 + v_plot.^2);

[X,Y] = meshgrid(x,y);

figure
contourf(X,Y,vel_mag,30,'LineColor','none')
colorbar
xlabel('x')
ylabel('y')
title('Velocity Magnitude Contour')
axis equal tight

%% Streamlines
figure
streamslice(X,Y,u_plot,v_plot)
xlabel('x')
ylabel('y')
title('Streamlines')
axis equal tight

%% u versus y at x-center
mid_x = round(Nx/2);

u_center = u_plot(:,mid_x);

figure
plot(y,u_center,'-o','LineWidth',1.5)
xlabel('y')
ylabel('u')
title('u velocity at x-center')
grid on

%% v versus x at y-center
mid_y = round(Nx/2);

v_center = v_plot(mid_y,:);

figure
plot(x,v_center,'-o','LineWidth',1.5)
xlabel('x')
ylabel('v')
title('v velocity at y-center')
grid on