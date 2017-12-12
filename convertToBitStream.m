function [stream, breakLength] = convertToBitStream(array_3d,N,m,n,numBlocks)
% Takes in a 3D array of an image and then returns the bitstream equivalent
% N is the group size of 8 by 8 blocks
    
    [m, n, numBlocks] = size(array_3d); % Rows, Columns, Number of Blocks
    assert(mod(numBlocks,N)==0);  % make sure we can make N groups of blocks
    streamDummy = reshape(array_3d,[m*n*numBlocks 1]); % converts to vector of ints
    streamDummyBi = de2bi(streamDummy); % Converts int to binary
    bits = size(streamDummyBi,2); % Gets number of bits
    stream = reshape(streamDummyBi',[1 bits*m*n*numBlocks]); % Turns to vec of binary
    breakLength = length(stream)/N; % Calculates break length
    
end