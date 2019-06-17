function stacked = stack(input)

    [m,n] = size(input);
    
    stacked = reshape(input,m*n,1);