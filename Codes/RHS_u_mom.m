function [RHS_u] = RHS_u_mom(i,j,gamma,del_x,del_y,delta_t,u_cell_face_sides_intp,v_cell_face_tb_intp,u_cell_cntr,P)
    
    % 3 terms : Conv + Pressure + Diff
    RHS_u_11 = (((u_cell_cntr(i,j) + u_cell_cntr(i,j+1))/2)*u_cell_face_sides_intp(i-1,j+1)*del_y);
    RHS_u_12 = (((u_cell_cntr(i,j) + u_cell_cntr(i,j-1))/2)*u_cell_face_sides_intp(i-1,j)*del_y);
    RHS_u_13 = (((u_cell_cntr(i,j) + u_cell_cntr(i-1,j))/2)*v_cell_face_tb_intp(i,j-1)*del_x);
    RHS_u_14 = (((u_cell_cntr(i,j) + u_cell_cntr(i+1,j))/2)*v_cell_face_tb_intp(i+1,j-1)*del_x);

    RHS_u_1 = -(RHS_u_11 - RHS_u_12 + RHS_u_13 - RHS_u_14);

    RHS_u_2 = - del_y * (P(i,j+1) - P(i,j-1))/2;

    RHS_u_31 = (((u_cell_cntr(i,j+1) - u_cell_cntr(i,j))/del_x)*del_y);
    RHS_u_32 = (((u_cell_cntr(i,j-1) - u_cell_cntr(i,j))/del_x)*del_y);
    RHS_u_33 = (((u_cell_cntr(i-1,j) - u_cell_cntr(i,j))/del_y)*del_x);
    RHS_u_34 = (((u_cell_cntr(i+1,j) - u_cell_cntr(i,j))/del_y)*del_x);

    RHS_u_3 = gamma * (RHS_u_31 + RHS_u_32 + RHS_u_33 + RHS_u_34);

    RHS_u = delta_t * (RHS_u_1 + RHS_u_2 + RHS_u_3) / (del_x*del_y);
end