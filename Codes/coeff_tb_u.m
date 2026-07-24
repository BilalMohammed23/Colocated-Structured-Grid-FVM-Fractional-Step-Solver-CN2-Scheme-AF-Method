function [bp,bs,bn,gp,gs,gn] = coeff_tb_u(i,j,del_y,delta_t,gamma,v_cell_face_tb_intp,u_cell_cntr)

bp = delta_t * (1/2) * ((2*gamma/(del_y^2)) + (v_cell_face_tb_intp(i,j-1) / (2*del_y)) - (v_cell_face_tb_intp(i+1,j-1) / (2*del_y)));

bs = delta_t * (1/2) * ((-gamma/(del_y^2)) - (v_cell_face_tb_intp(i+1,j-1) / (2*del_y)));

bn = delta_t * (1/2) * ((-gamma/(del_y^2)) + (v_cell_face_tb_intp(i,j-1) / (2*del_y)));

gp = delta_t * (1/2) * ((u_cell_cntr(i-1,j) / (4*del_y)) - (u_cell_cntr(i+1,j) / (4*del_y)));

gs = delta_t * (1/2) * (-(u_cell_cntr(i,j) / (4*del_y)) - (u_cell_cntr(i+1,j) / (4*del_y)));

gn = delta_t * (1/2) * ((u_cell_cntr(i,j) / (4*del_y)) + (u_cell_cntr(i-1,j) / (4*del_y)));

end