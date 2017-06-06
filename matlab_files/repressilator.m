classdef repressilator
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
        function obj = repressilator(L, W, initCond, D_pi, param, noiseParam, periodicity)
            % assign properties
            obj.length = L;
            obj.width = W;
            obj.initCond = initCond;
            obj.D_pi = D_pi;
            obj.param = param;
            obj.noiseParam = noiseParam;
            % initialize C_x and C_y as empty arrays
            obj.C_x = obj.CoeffMat(W, periodicity(1));
            obj.C_y = obj.CoeffMat(L, periodicity(2));

        end
        % main simulat method
        function output = simulate(obj, t, dt)
            %% unpack parameters

            % unpack geometry
            L = obj.length;
            W = obj.width;

            % assert dxy = 1
            dxy = 1;

            % unpack diffusion constants
            Dp1 = obj.D_pi(1);
            Dp2 = obj.D_pi(2);
            Dp3 = obj.D_pi(3);

            % unpack diff eq parameters
            alpha = obj.param(1);
            alpha0 = obj.param(2);
            beta = obj.param(3);

            % unpack noise parameters
            mnoise = obj.noiseParam(1);
            pnoise =obj.noiseParam(2);

            % unpack diffusion matrices
            C_x = obj.C_x;
            C_y = obj.C_y;

            % calculate time length of simulation
            numStepsT = ceil(t/dt);

            % initialize x and y arrays
            x = 1:W;
            y = 1:L;

            % declare concentration arrays

            m1 = zeros(L,W,numStepsT);
            m2 = zeros(L,W,numStepsT);
            m3 = zeros(L,W,numStepsT);

            p1 = zeros(L,W,numStepsT);
            p2 = zeros(L,W,numStepsT);
            p3 = zeros(L,W,numStepsT);

            % unpack initial concentrations
            m1(:,:,1) = obj.initCond(:,:,1);
            m2(:,:,1) = obj.initCond(:,:,2);
            m3(:,:,1) = obj.initCond(:,:,3);

            p1(:,:,1) = obj.initCond(:,:,4);
            p2(:,:,1) = obj.initCond(:,:,5);
            p3(:,:,1) = obj.initCond(:,:,6);

            % calculate noise random population
            m1noise = obj.noisegen(m1(:,:,1), mnoise);
            m2noise = obj.noisegen(m2(:,:,1), mnoise);
            m3noise = obj.noisegen(m3(:,:,1), mnoise);
            p1noise = obj.noisegen(p1(:,:,1), pnoise);
            p2noise = obj.noisegen(p2(:,:,1), pnoise);
            p3noise =obj.noisegen(p3(:,:,1), pnoise);

            %% evaluate PDEsusing diffp and diffm functions
            for t = 1:(numStepsT - 1)
                m1(:,:,t+1) = max(0, m1(:,:,t) + obj.diffm(m1(:,:,t), alpha, alpha0, p2(:,:,t), m1noise, dt));
                m2(:,:,t+1) = max(0, m2(:,:,t) + obj.diffm(m2(:,:,t), alpha, alpha0, p3(:,:,t), m2noise, dt));
                m3(:,:,t+1) = max(0, m3(:,:,t) + obj.diffm(m3(:,:,t), alpha, alpha0, p1(:,:,t), m3noise, dt));

                p1(:,:,t+1) = max(0, p1(:,:,t) + obj.diffp(p1(:,:,t), Dp1, beta, m1(:,:,t), dxy, dt, numStepsT, C_x, C_y, p1noise));
                p2(:,:,t+1) = max(0, p2(:,:,t) + obj.diffp(p2(:,:,t), Dp2, beta, m2(:,:,t), dxy, dt, numStepsT, C_x, C_y, p2noise));
                p3(:,:,t+1) = max(0, p3(:,:,t) + obj.diffp(p3(:,:,t), Dp3, beta, m3(:,:,t), dxy, dt, numStepsT, C_x, C_y, p3noise));
            end
            % return values
            output = struct('m1', m1, 'm2', m2, 'm3', m3, ...
                'p1', p1, 'p2', p2, 'p3', p3);
        end

        % construct diffusive coefficient matrix
        function C = CoeffMat(obj, X, periodic)
            % periodic flags should be used as follows:
            %       flag = 0 --> not periodic
            %       flag = 1 --> periodic

            %define coefficient matrix for y direction
            c = eye(X);
            for i=1:X
                for j=1:X
                    if i == j
                        c(i,j) = -2;
                    elseif i == j+1
                        c(i,j) = 1;
                    elseif i == j-1
                        c(i,j) = 1;
                    end
                end
            end

            %% apply BC

            % if not periodic in y
            if periodic == 0
                c(1,1) = -1;
                c(end,end) = -1;
            % else if periodic in y
            else
                % allow transer across edges
                c(1,end) = 1;
                c(end,1) = 1;
            end

            %apply sparse function
            C = sparse(c);
        end

        % diff eq for mRNAs
        function dm = diffm(obj, m, alpha, alpha0, Km, Kp, p, noise, dt)
            % calcualte differential
            dm = dt*(alpha ./ (1 + p.^2.1) - m  + alpha0 + noise);
        end

        % diff eq for proteins
        function dp = diffp(obj, p, diff, beta, Km, Kp, m, dx, dt, numStepsT, C_x, C_y, noise)
            % calculate differential
            dp = dt*(beta*(m - p) + ((diff)/(dx)^2)*(C_y*p+ p*C_x) + noise);
        end

        % noise matrix generator
        function noise = noisegen(obj, matrix, amplitude)
            % generate a random matrix with mean 0
            % and same shape as matrix
            seed = (rand(size(matrix)) - 0.5) *2;
            % return matrix with max amplitude equal
            % to amplitude
            noise = seed*amplitude;
        end
    end
end
