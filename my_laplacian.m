function u = my_laplacian(M,n,dx) %input: matrix;model domain; resolution
    
    i = 2:n-1; % interior rows
 
    u = zeros(n,1);
    
    %internal nodes
    u(i,1) = (M(i+1,1) + M(i-1,1) - 2 * M(i,1))./(dx^2);

    % upper and lower b.c.
    u([1,n],:) = (M([2,1],:) + M([n,n-1],:) - 2.*M([1,n],:))./dx^2;
    
end

