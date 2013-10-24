% Value of the European call option by Heston.
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
%   F is the forward price of the underlying object

function ret = CallHeston(F, K, T, t, parms)

    % parms = [ lambda eta rho vbar v0 ]
    tau = T - t;
    x = log(F/K);
    P1 = HestonPi(x, tau, 1, parms);
    P0 = HestonPi(x, tau, 0, parms);
    ret = K .* ( exp(x).*P1 - P0 );   