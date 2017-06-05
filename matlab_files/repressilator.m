classdef repressilator < handle
    properties
        length
        width
        initCond
        D_pi
        param
        noiseParam
        C_x
        C_y
    end
    
    methods
        %Constructor
        function obj = repressilator()
        
        end
       
        function simulate()
    
        end
        
        % construct diffusive coefficient matrix
        function setCoeffMat(periodic_x, periodic_y)
            % periodic flags should be used as follows:
            %       flag = 0 --> not periodic
            %       flag = 1 --> periodic
            
            % unpack length and width
            L = obj.length;
            W = obj.width;
            %define coefficient matrix for y direction
            c_y = eye(L);
            for i=1:L
                for j=1:L
                    if i == j
                        c_y(i,j) = -2;
                    elseif i == j+1
                        c_y(i,j) = 1;
                    elseif i == j-1
                        c_y(i,j) = 1;
                    end
                end
            end
            
            %define coefficient matrix for x direction
            c_x = eye(W);
            for i=1:W
                for j=1:W
                    if i == j
                        c_x(i,j) = -2;
                    elseif i == j+1
                        c_x(i,j) = 1;
                    elseif i == j-1
                        c_x(i,j) = 1;
                    end
                end
            end
            % apply BC
            
            % if not periodic in y
            if periodic_y == 0
                c_y(1,1) = -1;
                c_y(end,end) = -1;
            % else if periodic in y
            else
                % allow transer across edges
                c_y(1,end) = 1;
                c_y(end,1) = 1;            
            end
            
             % if not periodic in x
            if periodic_x == 0
                c_x(1,1) = -1;
                c_x(end,end) = -1;
            % else if periodic in x
            else
                % allow transer across edges
                c_x(1,end) = 1;
                c_x(end,1) = 1;            
            end
            
            
            %apply sparse function
            obj.C_y = sparse(c_y);
            obj.C_x = sparse(c_x);
        end
        
    
    % nonlinearity static methods
    methods(static)
        % diff eq for mRNAs
        function dm = diffm(m, alpha, alpha0, p, noise, C, dt)
            % calcualte differential
            dm = dt*(-m + alpha ./ (1 + p.^2) + alpha0 + noise);
        end
        
        % diff eq for proteins
        function dp = diffp(p, diff, beta, m, dx, dt, numStepsT, C)
            % calculate differential
            dp = dt*( beta*(m- p) + ((diff)/(dx)^2)*(C*p+ p*C) );
        end
        
        % noise matrix generator
        function noise = noisegen(matrix, amplitude)
            % generate a random matrix with mean 0
            % and same shape as matrix 
            seed = (rand(size(matrix)) - 0.5) *2
            % return matrix with max amplitude equal
            % to amplitude
            noise = seed*amplitude
        end
    
    end
end