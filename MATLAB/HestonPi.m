% Computation of the of P1 and P0, numerical integral.
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

function ret = HestonPi(x, tau, typo, parms)
    Q = @(u) HestonPiIntegrand(u, x, tau, typo, parms);
    ret = 1./2 + (1./pi).*quadl(Q, 0, 100);