function [RHS_f] = RHS_final(i,j,del_vstar,del_ustar_bar,gp,gs,gn)
    
    RHS_f_1 = del_ustar_bar(j,i);

    RHS_f_2 = - (gp * del_vstar(j,i)) - (gs * del_vstar(j+1,i)) - (gn * del_vstar(j-1,i));

    RHS_f = RHS_f_1 + RHS_f_2;
end