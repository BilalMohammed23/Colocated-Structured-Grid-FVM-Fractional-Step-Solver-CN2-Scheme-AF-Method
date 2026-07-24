function [del_t_P_sides,del_t_P_tb] = delt_P_faces(Nx,del_x,del_y,del_t_P_cell_cntr)

del_t_P_sides = zeros(Nx,Nx+1);
del_t_P_tb = zeros(Nx+1,Nx);

for i = 2 : Nx+1
    for j = 1 : Nx+1
        del_t_P_sides(i-1,j) = (del_t_P_cell_cntr(i,j+1) - del_t_P_cell_cntr(i,j))/del_x ; 
    end
end

for k = 1 : Nx+1
    for l = 2 : Nx+1
        del_t_P_tb(k,l-1) = (del_t_P_cell_cntr(k,l) - del_t_P_cell_cntr(k+1,l))/del_y;
    end
end

end