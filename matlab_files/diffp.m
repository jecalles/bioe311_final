function y = diffp(matrix, diff, beta, m, dx, dt, numStepsT, C)
    matrix(:,:,i+1) = matrix(:,:,i) + dt *( beta*(m(:,:,1)-matrix(:,:,i))...
        + ((diff)/(dx)^2)*(C*matrix(:,:,i)+ matrix(:,:,i)*C) );
y = matrix;