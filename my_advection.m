function drift_term = my_advection(B,chi,dx) %input: matrix;model domain; resolution

% % function drift_term = my_advection(B,chi,dx,dt,DB) %input: matrix;model domain; resolution

lnB = log(B + 1e-12);

% Periodic central differences in x-direction
dlnB_dx = (circshift(lnB, 1) - circshift(lnB, -1))./ (2 * dx);
dchi_dx = (circshift(chi, 1) - circshift(chi, -1))./ (2 * dx);

drift_term = dlnB_dx .* dchi_dx;




% % % % 
% % % % lnB = log(B + 1e-12);
% % % % 
% % % % q = B .* chi;                % trait mass density
% % % % Nx = length(B);
% % % % % Periodic neighbors
% % % % ip1 = mod((1:Nx), Nx) + 1;         % i+1
% % % % im1 = mod((1:Nx) - 2, Nx) + 1;     % i-1
% % % % 
% % % % % Step 1: Compute face-centered advection velocities v_{i+1/2}
% % % % dlnB = lnB(ip1) - lnB;            % lnB difference across each face
% % % % v = ((2 * DB).* dlnB)./ dx;           % velocity at each face
% % % % 
% % % % % Step 2: Compute upwind fluxes
% % % % F = zeros(size(B));
% % % % ind_pos = find(v >= 0);           % flow from cell i to i+1
% % % % ind_neg = find(v < 0);            % flow from cell i+1 to i
% % % % 
% % % % F(ind_pos) = v(ind_pos) .* q(ind_pos);
% % % % F(ind_neg) = v(ind_neg) .* q(ip1(ind_neg));
% % % % 
% % % % % Step 3: Update trait mass via conservative flux
% % % % q_new = q - (dt / dx).* (F - F(im1));
% % % % 
% % % % % Step 4: Recover new trait means
% % % % chi_new = q_new ./ B; %WARNING: dividing by B might cause problems as B can be really small! 
% % % % 
% % % % drift_term = chi_new - chi;

end

