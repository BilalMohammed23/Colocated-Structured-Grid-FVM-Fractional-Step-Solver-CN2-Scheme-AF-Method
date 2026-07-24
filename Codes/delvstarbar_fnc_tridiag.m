function [del_vstar_bar] = delvstarbar_fnc_tridiag(Nx,gamma,del_x,del_y,delta_t,u_cell_face_sides_intp,v_cell_face_tb_intp,v_cell_cntr,P,del_ustar_bar)
m1 = zeros(Nx,Nx);
m2 = zeros(Nx,1);
del_vstar_bar = zeros(Nx+2,Nx+2);

% tri-diagonal system of del_ustar_bar
for i = 1 : Nx
    for j = 1 : Nx
        [adash_p,adash_w,adash_e,hp,hw,he] = coeff_sides_v(i+1,j+1,del_x,delta_t,v_cell_cntr,gamma,u_cell_face_sides_intp);

        m2(j) = RHS_v_mom(i+1,j+1,gamma,del_x,del_y,delta_t,u_cell_face_sides_intp,v_cell_face_tb_intp,v_cell_cntr,P,del_ustar_bar,hp,hw,he);

        if j == 1
            m1(j,j) = (1 + adash_p - adash_w);
            m1(j,j+1) = adash_e;
        elseif j == Nx
            m1(j,j) = (1 + adash_p - adash_e);
            m1(j,j-1) = adash_w;
        else
            m1(j,j) = (1 + adash_p);
            m1(j,j+1) = adash_e;
            m1(j,j-1) = adash_w;
        end
    end
    m = (m1\m2).';
    for k = 1 : Nx
        del_vstar_bar(i+1,k+1) = m(k);
    end
    m1 = zeros(Nx,Nx);
end
del_vstar_bar(:,1)    = -del_vstar_bar(:,2);
del_vstar_bar(:,Nx+2) = -del_vstar_bar(:,Nx+1);
end