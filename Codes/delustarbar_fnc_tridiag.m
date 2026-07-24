function [del_ustar_bar] = delustarbar_fnc_tridiag(Nx,gamma,del_x,del_y,delta_t,u_cell_face_sides_intp,v_cell_face_tb_intp,u_cell_cntr,P)
m1 = zeros(Nx,Nx);
m2 = zeros(Nx,1);
del_ustar_bar = zeros(Nx+2,Nx+2);

% tri-diagonal system of del_ustar_bar
for i = 1 : Nx
    for j = 1 : Nx
        m2(j) = RHS_u_mom(i+1,j+1,gamma,del_x,del_y,delta_t,u_cell_face_sides_intp,v_cell_face_tb_intp,u_cell_cntr,P);

        [ap,aw,ae] = coeff_sides_u(i+1,j+1,gamma,del_x,delta_t,u_cell_face_sides_intp,u_cell_cntr);

        if j == 1
            m1(j,j) = (1 + ap - aw);
            m1(j,j+1) = ae;
        elseif j == Nx
            m1(j,j) = (1 + ap - ae);
            m1(j,j-1) = aw;
        else
            m1(j,j) = (1 + ap);
            m1(j,j+1) = ae;
            m1(j,j-1) = aw;
        end
    end
    m = (m1\m2).';
    for k = 1 : Nx
        del_ustar_bar(i+1,k+1) = m(k);
    end
    m1 = zeros(Nx,Nx);
end
del_ustar_bar(:,1)    = -del_ustar_bar(:,2);
del_ustar_bar(:,Nx+2) = -del_ustar_bar(:,Nx+1);
end