% Parameters calibration in the Heston Model
%   OptionData is the data from the market with risk-free rate r,
%   time to maturity T, s0 is the spot price at t = 0, the strike 
%   price K, Option value, bid and offer prices. 
%   lambda is the speed of the reversion mean of the volatility.
%   eta is the volatility of volatility.
%   rho is the correlation between the Brownian motion of St and vt.
%   vbar is the long-term mean of the volatility.
%   K is the strike price of the option.
%   vo is the initial variance.
%   s0 is the initial price of asset.
%   r is the risk-free rate.
%   T is the maturity time.

function hdl = HestonDifferences(PricingFunction, OptionData, DebugMode)

    NoOfOptions = size(OptionData, 1);
    NoOfIterations = 0;
    PriceDifference = zeros(NoOfOptions, 1);
    
    hdl = @HestonDifferences_handler;

    function ret = HestonDifferences_handler(input)

        for i = 1:NoOfOptions;
            call = PricingFunction( OptionData(i,1), OptionData(i,3), ...
                OptionData(i,2), 0, input );
            PriceDifference(i) = (OptionData(i,4) - call);  
            PriceDifference(i) = (OptionData(i,4) - call)/OptionData(i,5);  
        end

        % counts the number of iterations run to calibrate model.
        NoOfIterations = NoOfIterations + 1;

        if DebugMode
            msqrt = sqrt(PriceDifference'*PriceDifference)/NoOfOptions;
%            msqrt = PriceDifference'*PriceDifference;
            fprintf('[%5d] <%.10f> \n', NoOfIterations, msqrt);
        end

        ret = PriceDifference;
    end
end