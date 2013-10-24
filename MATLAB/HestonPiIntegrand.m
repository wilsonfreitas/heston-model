% Computation of the integrand of P1tio and P0tio.
% The notation is the same like in "The Volatility Surface", Jim Gatheral.
%   lambda is the speed of the reversion mean of the volatility.
%   eta is the volatility of volatility.
%   rho is the correlation between the Brownian motion of St and vt.
%   vbar is the long-term mean of the volatility.
%   K is the strike price of the option.
%   vo is the initial variance.
%   s0 is the initial price of asset.
%   r is the risk-free rate.
%   T is the maturity time.
%   t is the actual time.
%   type = 1 is P1 otherwise P0.
function ret = HestonPiIntegrand(u, x, tau, typo, parms)

    lambda = parms(1);
    eta = parms(2);
    rho = parms(3);
    vbar = parms(4);
    v0 = parms(5);
    
    alpha = -u.^2./2 - 1i.*u./2 + 1i.*typo.*u;
    beta = lambda - rho.*eta.*typo - rho.*eta.*1i.*u;
    gamma = eta.^2./2;
    d = sqrt(beta.^2 - 4.*alpha.*gamma);
    r_plus = (beta + d)./eta.^2;
    r_minus = (beta - d)./eta.^2; 
    g = r_minus./r_plus;

    D = r_minus.*( ( 1 - exp(-d.*tau) )./( 1 - g.*exp(-d.*tau) ) );
    C = lambda.*(r_minus.*tau - (2./eta.^2).*log( ( 1 - g.*exp(-d.*tau))./(1 - g)));

    ret = real( exp( C.*vbar + D.*v0 + 1i.*u.*x ) ./ (1i.*u) );