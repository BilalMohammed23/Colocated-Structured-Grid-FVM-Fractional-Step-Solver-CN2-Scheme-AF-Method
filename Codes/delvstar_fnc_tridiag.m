function [del_vstar] = delvstar_fnc_tridiag(Nx,gamma,del_y,delta_t,v_cell_face_tb_intp,v_cell_cntr,del_vstar_bar)
m1 = zeros(Nx,Nx);
m2 = zeros(Nx,1);
del_vstar = zeros(Nx+2,Nx+2);

% tri-diagonal system of del_vstar
for i = 1 : Nx
    for j = 1 : Nx
        [bdash_p,bdash_n,bdash_s] = coeff_tb_v(j+1,i+1,gamma,delta_t,del_y,v_cell_cntr,v_cell_face_tb_intp);

        m2(j) = del_vstar_bar(j+1,i+1);

        if j == 1
            m1(j,j) = (1 + bdash_p - bdash_n);
            m1(j,j+1) = bdash_s;
        elseif j == Nx
            m1(j,j) = (1 + bdash_p - bdash_s);
            m1(j,j-1) = bdash_n;
        else
            m1(j,j) = (1 + bdash_p);
            m1(j,j+1) = bdash_s;
            m1(j,j-1) = bdash_n;
        end
    end
    m = (m1\m2).';
    for k = 1 : Nx
        del_vstar(k+1,i+1) = m(k);
    end
    m1 = zeros(Nx,Nx);
end
del_vstar(1,:) = -del_vstar(2,:);
del_vstar(Nx+2,:) = -del_vstar(Nx+1,:);
end