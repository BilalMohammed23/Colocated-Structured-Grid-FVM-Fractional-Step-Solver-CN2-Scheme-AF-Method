function [bdash_p,bdash_n,bdash_s] = coeff_tb_v(i,j,gamma,delta_t,del_y,v_cell_cntr,v_cell_face_tb_intp)

    bdash_p = delta_t * (1/2) * ((2*gamma/(del_y^2)) + (v_cell_cntr(i-1,j)/(4*del_y)) - (v_cell_cntr(i+1,j)/(4*del_y)) + (v_cell_face_tb_intp(i,j-1)/(2*del_y)) - (v_cell_face_tb_intp(i+1,j-1)/(2*del_y)));
    
    bdash_s = delta_t * (1/2) * (-(gamma/(del_y^2)) - (v_cell_cntr(i,j)/(4*del_y)) - (v_cell_cntr(i+1,j)/(4*del_y)) - (v_cell_face_tb_intp(i+1,j-1)/(2*del_y)));
    
    bdash_n = delta_t * (1/2) * (-(gamma/(del_y^2)) + (v_cell_cntr(i,j)/(4*del_y)) + (v_cell_cntr(i-1,j)/(4*del_y)) + (v_cell_face_tb_intp(i,j-1)/(2*del_y)));
end