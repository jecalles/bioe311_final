function y = diffm(matrix, alpha, alpha0, p, dx, t, dt, noise, C);
    matrix(:,:,t+1) = matrix(:,:,t) + dt *( -matrix(:,:,t) + alpha ./ (1 + p(:,:,t).^2) + alpha0 + noise );
    % remove negative conc.
    matrix(:,:,t+1) = max(0, matrix(:,:,t+1));
    y = matrix;
end