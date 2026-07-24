function [ap,aw,ae] = coeff_sides_u(i,j,gamma,del_x,delta_t,u_cell_face_sides_intp,u_cell_cntr)

ap = delta_t * (1/2) * ((2*gamma/(del_x^2)) + (u_cell_cntr(i,j+1)/(4*del_x)) - (u_cell_cntr(i,j-1)/(4*del_x)) + (u_cell_face_sides_intp(i-1,j+1)/(2*del_x)) - (u_cell_face_sides_intp(i-1,j)/(2*del_x))); 

aw = delta_t * (1/2) * (-(gamma/(del_x^2)) - (u_cell_cntr(i,j)/(4*del_x)) - (u_cell_cntr(i,j-1)/(4*del_x)) - (u_cell_face_sides_intp(i-1,j)/(2*del_x)));

ae = delta_t * (1/2) * (-(gamma/(del_x^2)) + (u_cell_cntr(i,j)/(4*del_x)) + (u_cell_cntr(i,j+1)/(4*del_x)) + (u_cell_face_sides_intp(i-1,j+1)/(2*del_x)));

end