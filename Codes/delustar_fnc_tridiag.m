function [del_ustar] = delustar_fnc_tridiag(Nx,gamma,del_y,delta_t,v_cell_face_tb_intp,u_cell_cntr,del_vstar,del_ustar_bar)
m1 = zeros(Nx,Nx);
m2 = zeros(Nx,1);
del_ustar = zeros(Nx+2,Nx+2);

% tri-diagonal system of del_ustar
for i = 1 : Nx
    for j = 1 : Nx
        [bp,bs,bn,gp,gs,gn] = coeff_tb_u(j+1,i+1,del_y,delta_t,gamma,v_cell_face_tb_intp,u_cell_cntr);

        m2(j) = RHS_final(i+1,j+1,del_vstar,del_ustar_bar,gp,gs,gn);

        if j == 1
            m1(j,j) = (1 + bp - bn);
            m1(j,j+1) = bs;
        elseif j == Nx
            m1(j,j) = (1 + bp - bs);
            m1(j,j-1) = bn;
        else
            m1(j,j) = (1 + bp);
            m1(j,j+1) = bs;
            m1(j,j-1) = bn;
        end
    end
    m = (m1\m2).';
    for k = 1 : Nx
        del_ustar(k+1,i+1) = m(k);
    end
    m1 = zeros(Nx,Nx);        
end

end