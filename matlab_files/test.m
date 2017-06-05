function test(periodic_x, periodic_y)
    % periodic flags should be used as follows:
    %       flag = 0 --> not periodic
    %       flag = 1 --> periodic
    
    % unpack length and width
    L = 5;
    W = 6;
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
    c_x
    c_y
 end