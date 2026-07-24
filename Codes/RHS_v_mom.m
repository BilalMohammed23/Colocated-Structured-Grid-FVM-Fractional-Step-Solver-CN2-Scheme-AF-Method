function [RHS_v] = RHS_v_mom(i,j,gamma,del_x,del_y,delta_t,u_cell_face_sides_intp,v_cell_face_tb_intp,v_cell_cntr,P,del_ustar_bar,hp,hw,he)

    % 3 terms : Conv + Pressure + Diff
    RHS_v_11 = (((v_cell_cntr(i,j) + v_cell_cntr(i,j+1))/2)*u_cell_face_sides_intp(i-1,j+1)*del_y);
    RHS_v_12 = (((v_cell_cntr(i,j) + v_cell_cntr(i,j-1))/2)*u_cell_face_sides_intp(i-1,j)*del_y);
    RHS_v_13 = (((v_cell_cntr(i,j) + v_cell_cntr(i-1,j))/2)*v_cell_face_tb_intp(i,j-1)*del_x);
    RHS_v_14 = (((v_cell_cntr(i,j) + v_cell_cntr(i+1,j))/2)*v_cell_face_tb_intp(i+1,j-1)*del_x);
    
    RHS_v_1 = -(RHS_v_11 - RHS_v_12 + RHS_v_13 - RHS_v_14);
    
    RHS_v_2 = - del_x * (P(i-1,j) - P(i+1,j))/2;
    
    RHS_v_31 = (((v_cell_cntr(i,j+1) - v_cell_cntr(i,j))/del_x)*del_y);
    RHS_v_32 = (((v_cell_cntr(i,j-1) - v_cell_cntr(i,j))/del_x)*del_y);
    RHS_v_33 = (((v_cell_cntr(i-1,j) - v_cell_cntr(i,j))/del_y)*del_x);
    RHS_v_34 = (((v_cell_cntr(i+1,j) - v_cell_cntr(i,j))/del_y)*del_x);
    
    RHS_v_3 = gamma * (RHS_v_31 + RHS_v_32 + RHS_v_33 + RHS_v_34);
    
    RHS_v_4 = - (hp * del_ustar_bar(i,j)) - (hw * del_ustar_bar(i,j-1)) - (he * del_ustar_bar(i,j+1));
    
    RHS_v = (delta_t * (RHS_v_1 + RHS_v_2 + RHS_v_3) / (del_x*del_y))  + RHS_v_4;

end