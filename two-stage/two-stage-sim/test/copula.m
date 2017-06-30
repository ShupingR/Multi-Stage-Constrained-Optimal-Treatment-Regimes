function [retSample] =  copula(x, y , z, M)
    % three dimension copula 
    %xTrans = norminv(tiedrank(x) ./ (length(x)+1)); 
    %yTrans = norminv(tiedrank(y) ./ (length(y)+1)); 

    xTrans = norminv(ceil(tiedrank(x)) ./ (length(x)+1)); 
    yTrans = norminv(ceil(tiedrank(y)) ./ (length(y)+1)); 
    zTrans = norminv(ceil(tiedrank(z)) ./ (length(z)+1)); 
    xyzTrans = [ xTrans, yTrans, zTrans];
    copulaCov = cov( xyzTrans );
    xyzSample = mvnrnd( [0 0 0], copulaCov, M );
    retSample = xyzSample;
    for m  = 1:M
        [ ~ , iX] = min(abs(xyzSample(m, 1) - xTrans));
        retSample(m, 1) = x( iX );
        [ ~ , iY] = min(abs(xyzSample(m, 2) - yTrans));
        retSample(m, 2) = y( iY );
        [ ~ , iZ] = min(abs(xyzSample(m, 3) - zTrans));
        retSample(m, 3) = z( iZ );
    end
end