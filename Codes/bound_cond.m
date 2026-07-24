function [u_cell_cntr,v_cell_cntr,u_cell_face_sides_intp,v_cell_face_sides_intp,u_cell_face_tb_intp,v_cell_face_tb_intp] = bound_cond(u_cell_cntr,v_cell_cntr,u_cell_face_sides_intp,v_cell_face_sides_intp,u_cell_face_tb_intp,v_cell_face_tb_intp);

    U_lid = 1;
    
    u_cell_cntr(1,:) = 2*U_lid - u_cell_cntr(2,:);
    u_cell_cntr(:,1) = -u_cell_cntr(:,2);
    u_cell_cntr(end,:) =  -u_cell_cntr(end-1,:);
    u_cell_cntr(:,end) = -u_cell_cntr(:,end-1);
    
    v_cell_cntr(1,:) = -v_cell_cntr(2,:);
    v_cell_cntr(:,1) = -v_cell_cntr(:,2);
    v_cell_cntr(end,:) = -v_cell_cntr(end-1,:);
    v_cell_cntr(:,end) = -v_cell_cntr(:,end-1);
    
    u_cell_face_sides_intp(:,1) = 0;
    u_cell_face_sides_intp(:,2) = 0;
    u_cell_face_sides_intp(:,end) = 0;
    u_cell_face_sides_intp(:,end-1) = 0;
    
    v_cell_face_sides_intp(:,1) = 0;
    v_cell_face_sides_intp(:,2) = 0;
    v_cell_face_sides_intp(:,end) = 0;
    v_cell_face_sides_intp(:,end-1) = 0;
    
    u_cell_face_tb_intp(1,:) = 0;
    u_cell_face_tb_intp(2,:) = U_lid;
    u_cell_face_tb_intp(end,:) = 0;
    u_cell_face_tb_intp(end-1,:) = 0;
    
    v_cell_face_tb_intp(1,:) = 0;
    v_cell_face_tb_intp(2,:) = 0;
    v_cell_face_tb_intp(end,:) = 0;
    v_cell_face_tb_intp(end-1,:) = 0;

end