function q = pressure_poisson_GS(Nx,del_x,del_y,RHS_P_pred)

q = zeros(Nx+2,Nx+2);  

Ax = del_y/(2*del_x);
Ay = del_x/(2*del_y);

% Ax = del_y/(del_x);
% Ay = del_x/(del_y);

denom = 2*Ax + 2*Ay;

tol = 1e-8;
max_iter = 10000;

for iter = 1:max_iter

    q_old = q;

    % neumann pressure BC: dP/dn = 0
    q(1,:)     = q(2,:);       
    q(end,:)   = q(end-1,:);   
    q(:,1)     = q(:,2);       
    q(:,end)   = q(:,end-1);   

    % Gauss-Seidel update for interior cells
    for i = 2:Nx+1
        for j = 2:Nx+1
            q(i,j) = (Ax*(q(i,j+1) + q(i,j-1)) + Ay*(q(i-1,j) + q(i+1,j)) - RHS_P_pred(i-1,j-1) ) / denom;
        end
    end

    % fixing pressure singularity
    q(2,2) = 0;

    % convergence check
    err = max(abs(q(:) - q_old(:)));

    if err < tol
        fprintf('Pressure GS converged in %d iterations\n',iter);
        break
    end

end

end