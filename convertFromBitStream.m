function data = convertFromBitStream(bitStream,m,n,numBlocks,bits)
% Bitstream should be a column vector
% Will return a 3D array of integers
    
    bitStream = uint8(bitStream);
    dataDumBi = reshape(bitStream, [bits m*n*numBlocks])'; % gets binary formatting
    dataDum = bi2de(dataDumBi); % turns into vec of ints
    data = reshape(dataDum,[m n numBlocks]); % turns to desired block shape
end