function drift_term = my_advection(B,chi,dx) %input: matrix;model domain; resolution

% % function drift_term = my_advection(B,chi,dx,dt,DB) %input: matrix;model domain; resolution

lnB = log(B + 1e-12);

% Periodic central differences in x-direction
dlnB_dx = (circshift(lnB, 1) - circshift(lnB, -1))./ (2 * dx);
dchi_dx = (circshift(chi, 1) - circshift(chi, -1))./ (2 * dx);

drift_term = dlnB_dx .* dchi_dx;



end

