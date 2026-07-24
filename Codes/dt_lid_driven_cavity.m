function [delta_t_max] = dt_lid_driven_cavity(u_cell_cntr,v_cell_cntr,del_x,del_y,Nx,gamma)
    CFL = 0.5;      
    Fourier_max = 5;      
                      
    umax = max(max(abs(u_cell_cntr(2:Nx+1,2:Nx+1))));
    vmax = max(max(abs(v_cell_cntr(2:Nx+1,2:Nx+1))));

    if umax < 1e-12
        umax = 1e-12;
    end
    if vmax < 1e-12
        vmax = 1e-12;
    end

    dt_conv    = 1 / (umax/del_x + vmax/del_y);
    dt_diff    = Fourier_max / (2*gamma*(1/del_x^2 + 1/del_y^2));

    delta_t_max = CFL * min(dt_conv, dt_diff);
end