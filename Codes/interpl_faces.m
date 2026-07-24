function [u_cell_face_sides_intp,v_cell_face_sides_intp,u_cell_face_tb_intp,v_cell_face_tb_intp] = interpl_faces(Nx,u_cell_cntr,v_cell_cntr,u_cell_face_sides_intp,v_cell_face_sides_intp,u_cell_face_tb_intp,v_cell_face_tb_intp)
    for i = 1 : Nx
        for j = 3 : Nx+1
            u_cell_face_sides_intp(i,j) = (u_cell_cntr(i+1,j) + u_cell_cntr(i+1,j-1))/2;
            v_cell_face_sides_intp(i,j) = (v_cell_cntr(i+1,j) + v_cell_cntr(i+1,j-1))/2;
        end
    end
    
    for k = 3 : Nx+1
        for l = 1 : Nx
            u_cell_face_tb_intp(k,l) = (u_cell_cntr(k,l+1) + u_cell_cntr(k-1,l+1))/2;
            v_cell_face_tb_intp(k,l) = (v_cell_cntr(k,l+1) + v_cell_cntr(k-1,l+1))/2;
        end
    end
end