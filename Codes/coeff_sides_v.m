function [adash_p,adash_w,adash_e,hp,hw,he] = coeff_sides_v(i,j,del_x,delta_t,v_cell_cntr,gamma,u_cell_face_sides_intp)

adash_p = delta_t * (1/2) * ((2*gamma/(del_x^2)) + (u_cell_face_sides_intp(i-1,j+1) / (2 * del_x)) - (u_cell_face_sides_intp(i-1,j) / (2 * del_x)));

adash_w = delta_t * (1/2) * (-(gamma/(del_x^2)) - (u_cell_face_sides_intp(i-1,j) / (2 * del_x)));

adash_e = delta_t * (1/2) * (-(gamma/(del_x^2)) + (u_cell_face_sides_intp(i-1,j+1) / (2 * del_x)));

hp = delta_t * (1/2) * ((v_cell_cntr(i,j+1) / (4*del_x)) - (v_cell_cntr(i,j-1) / (4*del_x)));

hw = - delta_t * (1/2) * ((v_cell_cntr(i,j) / (4*del_x)) + (v_cell_cntr(i,j-1) / (4*del_x)));

he = delta_t * (1/2) * ((v_cell_cntr(i,j) / (4*del_x)) + (v_cell_cntr(i,j+1) / (4*del_x)));

end