function [h, pval, Jstat, Dev, Md, nu] = SynchTest_static(theta_star, X_star, odds, ordKidx, sigma, alpha)
    T = size(X_star, 1);
    M = numel(theta_star);
    Mr = length( ordKidx );
    Md = M - Mr;
    
    XR_star = X_star(:, ordKidx);
    thetaR_star = zeros(Mr, 1);
    
    SrCompidx = 1:M; SrCompidx(ordKidx) = [];
    for gditer=1:2000
        grad = sum(XR_star)' - T*( ( exp(thetaR_star) )/( 1+sum(exp(thetaR_star))+sum(odds(SrCompidx)) ) );
        thetaR_star = thetaR_star + alpha*grad;
    end
    
    thetaR_star_tmp = thetaR_star;

    thetaR_star = log(odds);
    thetaR_star(ordKidx) = thetaR_star_tmp;
    
    ll_diff = (sum(X_star)*theta_star - T*log( 1+sum(exp(theta_star)) )) - (sum(X_star)*thetaR_star - T*log( 1+sum(exp(thetaR_star(:))) ));

    Dev = 2*ll_diff;    
    nu = max(Dev - Md, 0);
    
    h = and(Md>0, (1 - sigma) < chi2cdf(Dev, Md));
    pval = (1 - chi2cdf(Dev, Md));
    Jstat = h * ( 1 - sigma - ncx2cdf( chi2inv( 1-sigma, Md), Md, nu ) );    
    
end