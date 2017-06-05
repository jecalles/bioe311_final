function initialRandom = randomize2(matrix, amplitude, L)
    seed = amplitude*rand(L,L);
    initialRandom = seed - mean(seed(:));
    matrix = initialRandom;
end