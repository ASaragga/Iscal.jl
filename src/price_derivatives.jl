#=
Price Derivative Instruments(12/12): 
OK  binprice    - Binomial put and call American option pricing using Cox-Ross-Rubinstein model
OK  blkimpv	    - Implied volatility for futures options from Black model
OK  blkprice	- Black model for pricing futures options
OK  blsdelta	- Black-Scholes sensitivity to underlying price change
OK  blsgamma	- Black-Scholes sensitivity to underlying delta change
OK  blsimpv	    - Black-Scholes implied volatility
OK  blslambda	- Black-Scholes elasticity
OK  blsprice	- Black-Scholes put and call option pricing
OK  blsrho	    - Black-Scholes sensitivity to interest-rate change
OK  blstheta	- Black-Scholes sensitivity to time-until-maturity change
OK  blsvega	    - Black-Scholes sensitivity to underlying price volatility
OK  opprofit	- Option profit
===============================
OK  mcprice     - Monte Carlo simulation put and call option pricing. At present only implemented the Longstaff–Schwartz least-squares Monte Carlo simulation method
OK  binprice_eu - Prices European options using a binomial tree and returns the difference to Black-Scholes price
OK  blspsi      - Black-Scholes sensitivity to underlying dividend yield change
OK  blsvanna    - Black-Scholes sensitivity once to the underlying spot price and once to volatility
OK  bingreeks   - Binomial Tree option price and the greeks: delta, rho, vega, theta, psi
OK  mcgreeks    - Monte Carlo simulation option price and the greeks: delta, rho, vega, theta, psi
=#

#= 
opprofit returns the profit of an option.

    profit = opprofit(assetprice, strike, cost, position, optiontype)

    Input arguments:
        assetprice — Asset price
        strike — Strike or exercise price
        cost — Cost of option
        position — Option position. Using the values :long or :short
        optiontyle - Option type. Using the values :call or :put
=#
function opprofit(assetprice, strike, cost, position, optiontype)
    if optiontype == :call
        if position == :long
            return  max(assetprice - strike, 0.0) - cost
        elseif position == :short
            return -max(assetprice - strike, 0.0) + cost
        end
    elseif optiontype == :put
        if position == :long
            return  max(strike - assetprice, 0.0) - cost
        elseif position == :short
            return -max(strike - assetprice, 0.0) + cost
        end
    end
end

#= 
blsprice computes European put and call option prices using a Black-Scholes model. 

    (call,put) = blsprice(S0, K, r, T, sigma, d = 0.0)

    Input arguments:
        S0      — Current price of underlying asset
        K       — Exercise price of the option
        r       — Annualized continuously compounded risk-free rate of return over life of the option
        T       — Time to expiration of option
        sigma   — Annualized asset price volatility
        d       — Annualized continuously compounded yield of underlying asset over life of the option    Default: 0.0.   
    Output parameters:
        call - Price (i.e., value) of a European call option.
        put  - Price (i.e., value) of a European put option.
=#
function blsprice(S0, K, r, T, sigma, d = 0.0)
	d1 = (log(S0/K)+(r - d + sigma^2/2)*T) / (sigma * sqrt(T))
    d2 = d1 - sigma * sqrt(T)
	call = S0 * exp(-d*T) * cdf(Normal(0,1),d1) - K * exp(-r*T) * cdf(Normal(0,1),d2);
	put = -S0 * exp(-d*T) * cdf(Normal(0,1),-d1) + K * exp(-r*T) * cdf(Normal(0,1),-d2);
    return call, put
end

#=
blkprice computes European put and call futures option prices using Black's model.

    blkprice(F0, K, r, T, sigma)

=#
function blkprice(F0, K, r, T, sigma)
	return blsprice(F0, K, r, T, sigma, r)
end

#=
blsdelta returns delta, the sensitivity in option value to change in the underlying asset price. Delta is also known as the hedge ratio.

        (calldelta, putdelta) = blsdelta(S0, K, r, T, sigma, d = 0.0) 

blsdelta can handle other types of underlies like Futures and Currencies. When pricing Futures (Black model), enter the input argument d as:
        d = r
When pricing currencies (Garman-Kohlhagen model), enter the input argument d as:
        d = ForeignRate
where ForeignRate is the continuously compounded, annualized risk-free interest rate in the foreign country.
=#
function blsdelta(S0, K, r, T, sigma, d = 0.0) 
	d1 = (log(S0/K) + (r - d + sigma^2/2) * T)/(sigma * sqrt(T))
	calldelta = exp(-d*T) * cdf(Normal(0,1),d1)
	putdelta = -exp(-d*T) * cdf(Normal(0,1),-d1)
    return calldelta, putdelta
end

#=
blsgamma returns gamma, the sensitivity of delta to change in the underlying asset price. 

        gamma = blsgamma(S0, K, r, T, sigma, d = 0.0)

=#
function blsgamma(S0, K, r, T, sigma, d = 0.0) 
	d1 = (log(S0/K) + (r - d + sigma^2/2) * T)/(sigma * sqrt(T))
	gamma = exp(-d*T) * pdf(Normal(0,1), d1)/(S0 * sigma * sqrt(T))
    return gamma
end


#=
blsvega returns vega, the rate of change of the option value with respect to the volatility of the underlying asset. 

        blsvega(S0, K, r, T, sigma, d = 0.0)

=#
function blsvega(S0, K, r, T, sigma, d = 0.0)
	d1 = (log(S0/K) + (r - d + sigma^2/2)*T)/(sigma * sqrt(T));
	vega = S0 * exp(-d*T) * pdf(Normal(0,1),d1) * sqrt(T);
    return vega
end


#=
blsrho returns the call option rho - callrho -, and the put option rho - putrho -. Rho is the rate of change in value of derivative securities with respect to interest rates. 


=#
function blsrho(S0, K, r, T, sigma, d = 0.0)
	d2=(log(S0/K) + (r - d - sigma^2/2)*T)/(sigma * sqrt(T))
	callrho = K * exp(-r*T) * cdf(Normal(0,1),d2) * T
	putrho = -K * exp(-r*T) * cdf(Normal(0,1),-d2) * T
    return callrho, putrho
end


#=
blstheta returns the call option theta - calltheta -, and the put option theta - puttheta -. Theta is the sensitivity in option value with respect to time and is measured in years. CallTheta or PutTheta can be divided by 365 to get Theta per calendar day or by 252 to get Theta by trading day.

        (calltheta, puttheta) = blstheta(S0, K, r, T, sigma, d = 0.0)

=#
function blstheta(S0, K, r, T, sigma, d = 0.0)
	d1 = (log(S0/K) + (r - d + sigma^2/2) * T)/(sigma * sqrt(T))
	shift = -exp(-d * T) * S0 * pdf(Normal(0,1),d1) * sigma/2 / sqrt(T)
	t1 = r * K * exp(-r*T)
	t2 = d * S0 * exp(-d*T)
    coef = cdf(Normal(0,1), d1 - sigma * sqrt(T))
	calltheta = shift - t1 * coef + t2 * cdf(Normal(0,1),d1)
	puttheta = shift + t1 * (1 - coef) + t2 * (cdf(Normal(0,1),d1) - 1)
    return calltheta, puttheta
end


#=
blslambda returns the elasticity of an option. calllambda is the call option elasticity or leverage factor, and putlambda is the put option elasticity or leverage factor. Elasticity (the leverage of an option position) measures the percent change in an option price per 1 percent change in the underlying asset price. 

        (calllambda, putlambda) = blslambda(S0, K, r, T, sigma, d = 0.0)

=#
function blslambda(S0, K, r, T, sigma, d = 0.0)
	(call, put) = blsprice(S0, K, r, T, sigma, d)
	(calldelta, putdelta) = blsdelta(S0, K, r, T, sigma, d)
	calllambda = calldelta * S0/call
    putlambda = putdelta * S0/put 
    return calllambda, putlambda
end


#=
blsimpv using a Black-Scholes model computes the implied volatility of an underlying asset from the market value of European options.

callimplsigma or callimplsigma = blsimpv(S0, K, r, T, Price, d, class, root_finding_method)

=#
function blsimpv(S0, K, r, T, Price, d, class, rm::Robust)
    if class == :call
        fcall(x)= (blsprice(S0, K, r, T, x, d)[1]-Price)
        callimplsigma = find_zero(fcall, rm.x0, Order0())
        return callimplsigma
    elseif class == :put
        fput(x) = (blsprice(S0, K, r, T, x, d)[2]-Price)
        putimplsigma = find_zero(fput, rm.x0, Order0())     
        return putimplsigma
    end
end

function blsimpv(S0, K, r, T, Price, d, class, rm::Secant)
    if class == :call 
        fcall(x)= (blsprice(S0, K, r, T, x, d)[1]-Price)
        callimplsigma = find_zero(fcall, rm.x0, Order1())
        return callimplsigma
    elseif class == :put 
        fput(x) = (blsprice(S0, K, r, T, x, d)[2]-Price)
        putimplsigma = find_zero(fput, rm.x0, Order1())     
        return putimplsigma
    end
end

function blsimpv(S0, K, r, T, Price, d, class, rm::Steffensen)
    if class == :call 
        fcall(x)= (blsprice(S0, K, r, T, x, d)[1]-Price)
        callimplsigma = find_zero(fcall, rm.x0, Order2())  # Steffensen
        return callimplsigma
    elseif class == :put    
        fput(x) = (blsprice(S0, K, r, T, x, d)[2]-Price)
        putimplsigma = find_zero(fput, rm.x0, Order2())    # Steffensen   
        return putimplsigma
    end
end

function blsimpv(S0, K, r, T, Price, d, class, rm::Brent) 
    if class == :call
        fcall(x)= (blsprice(S0, K, r, T, x, d)[1]-Price)
        callimplsigma = find_zero(fcall, (rm.inflimit,rm.suplimit), Roots.Brent()) 
        return callimplsigma
    elseif class == :put 
        fput(x) = (blsprice(S0, K, r, T, x, d)[2]-Price)
        putimplsigma = find_zero(fput, (0.001, 1.2), Roots.Brent())     
        return putimplsigma
    end
end

function blsimpv(S0, K, r, T, Price, d, class)
    if class == :call
        fcall(x)= (blsprice(S0, K, r, T, x, d)[1]-Price)
        callimplsigma = find_zero(fcall, 0.20, Order1())
        return callimplsigma
    elseif class == :put
        fput(x) = (blsprice(S0, K, r, T, x, d)[2]-Price)
        putimplsigma = find_zero(fput, 0.20, Order1())     
        return putimplsigma
    end
end

function blsimpv(S0, K, r, T, Price, class) 
    if class == :call
        fcall(x)= (blsprice(S0, K, r, T, x, 0.0)[1]-Price)
        callimplsigma = find_zero(fcall, 0.20, Order1())
        return callimplsigma
    elseif class == :put
        fput(x) = (blsprice(S0, K, r, T, x, 0.0)[2]-Price)
        putimplsigma = find_zero(fput, 0.20, Order1())     
        return putimplsigma
    end
end

#=
blkipmv computes the implied volatility of a futures price from the market value of European futures options using Black's model

        blkimpv(F0, K, r, T, Price, class, root_finding_method)

=#
function blkimpv(F0, K, r, T, Price, class, method)
	return blsimpv(F0, K, r, T, Price, r, class, method)
end
blkimpv(F0, K, r, T, Price, class) = blkimpv(F0, K, r, T, Price, class, Secant(0.20))

#=
blspsi returns Psi (also know as epsilon) that is the percentage change in option value per percentage change in the underlying dividend yield,

        function blspsi(S0, K, r, T, sigma, d = 0.0)

=#
function blspsi(S0, K, r, T, sigma, d = 0.0)
	d1 = (log(S0/K) + (r - d + sigma^2/2)*T)/(sigma * sqrt(T))
	callpsi = -S0 * exp(-d*T) * cdf(Normal(0,1),d1) * T
	putpsi = S0 * exp(-d*T) * cdf(Normal(0,1),-d1) * T
    return callpsi, putpsi
end

#=
blsvanna is a second order derivative of the option value, once to the underlying spot price and once to volatility (source: wikipedia)

        vanna = blsvanna(S0, K, r, T, sigma, d = 0.0)

=#
function blsvanna(S0, K, r, T, sigma, d = 0.0)
	d1 = (log(S0/K) + (r - d + sigma^2/2)*T)/(sigma*sqrt(T))
	d2 = d1 - sigma * sqrt(T)
	vanna = - exp(-d*T) * pdf(Normal(0,1),d1) * d2/sigma
    return vanna
end

#=
binprice computes the prices of an American option using the Cox-Ross-Rubinstein binomial pricing model. An American option can be exercised any time until its expiration date.

    optionvalue = binprice(S, K, r, T, sigma, q, class::Symbol, n::Int64)

Input arguments
    S       - Price
    K       - Strike
    r       - Interest rate
    T       - Time
    sigma   - Volatility
    q       - Dividend rate, specified as a scalar decimal.
    class   - Either :call, :put  
Output arguments:
    optionvalue - Option value, returned as a vector that represents each node of the Cox-Ross-Rubinstein (CRR) binary tree.
=#
        
# Returns the price for EITHER the :call and :put American option
binprice(S, K, r, sigma, T, dyield, class, method) = binprice2.(S, K, r, sigma, T, dyield, class, Ref(method))
function binprice2(S, K, r, sigma, t, dyield, class::Symbol, method::CRR) 
    N = method.steps
    deltat = t/N
    act = exp(-r * deltat)
    
    U = exp(sigma * sqrt(deltat))
    D = 1/U
    
    R = exp((r-dyield) * deltat)
    p = (R-D)/(U-D)
    q = (U-R)/(U-D)
    
    if class == :put
        Z = [max(0, K - S*exp((2*i-N) * sigma * sqrt(deltat))) for i = 0:N]
        for n = N-1:-1:0
            for i = 0:n
                x = K - S * exp((2*i-n) * sigma * sqrt(deltat))
                y = act * (q * Z[i+1] + p * Z[i+2])
                Z[i+1] = max(x,y)
            end
        end
    elseif class == :call
        Z = [max(0, S*exp((2*i-N) * sigma * sqrt(deltat)) - K) for i = 0:N]
        for n = N-1:-1:0
            for i = 0:n
                x = S * exp((2*i-n) * sigma * sqrt(deltat)) - K
                y = act * (q * Z[i+1] + p * Z[i+2])
                Z[i+1] = max(x,y)
            end
        end
    end
    return Z[1]
end # Adapted from: https://wilmott.com/automatic-for-the-greeks/
binprice(S, K, r, sigma, T, dyield, class::Symbol) = binprice(S, K, r, sigma, T, dyield, class, CRR(1000))



# Returns the price for BOTH :call and :put American options
binprice(S, K, r, sigma, T, dyield, method) = binprice2.(S, K, r, sigma, T, dyield, Ref(method))
binprice(S, K, r, sigma, T, dyield) = binprice2.(S, K, r, sigma, T, dyield, Ref(CRR(1000))) 
function binprice2(S, K, r, sigma, t, dyield, method::CRR) 
    N = method.steps
    deltat = t/N
    act = exp(-r * deltat)

    U = exp(sigma * sqrt(deltat))
    D = 1/U
    R = exp((r-dyield) * deltat)
    p = (R-D)/(U-D)
    q = (U-R)/(U-D)
    
    Zc = [max(0, S*exp((2*i-N) * sigma * sqrt(deltat)) - K) for i = 0:N]
    Zp = [max(0, K - S*exp((2*i-N) * sigma * sqrt(deltat))) for i = 0:N]

    for n = N-1:-1:0
        for i = 0:n
            xc = S * exp((2*i-n) * sigma * sqrt(deltat)) - K
            yc = act * (q * Zc[i+1] + p * Zc[i+2])
            Zc[i+1] = max(xc,yc)
            xp = -xc
            yp = act * (q * Zp[i+1] + p * Zp[i+2])
            Zp[i+1] = max(xp,yp)
        end
    end
    return Zc[1], Zp[1]
end


#= binprice_eu prices European options using a binomial tree method.

    (price, diff) = binprice_eu(S, K, r, sigma, T, q, class::Symbol, method)

Input arguments:
    q - carry-cost, for equities the dividend yield
    class - Either :call or :put
    method - Binomial tree method. E.g. CRR(100) for CRR method with 100 steps

Output arguments:
    price - Price of the Call or Put
    diff  - Difference between the price of a call or put computed by a binomial tree method and the price computed by the Black-Scholes formula.
=#
binprice_eu(S, K, r, sigma, T, dyield, class::Symbol, method) = binprice_eu2.(S, K, r, sigma, T, dyield, class::Symbol, Ref(method))
binprice_eu(S, K, r, sigma, T, dyield, class::Symbol) = binprice_eu2.(S, K, r, sigma, T, dyield, class::Symbol, Ref(CRR(1000)))
# Returns the price for EITHER the :call or :put European option
function binprice_eu2(S, K, r, sigma, t, dyield, class::Symbol, method::CRR) 
    N = method.steps
    deltat = t/N
    act = exp(-r * deltat)
    
    U = exp(sigma * sqrt(deltat))
    D = 1/U
    
    R = exp((r-dyield) * deltat)
    p = (R-D)/(U-D)
    q = (U-R)/(U-D)
    
    if class == :put
        Z = [max(0, K - S*exp((2*i-N) * sigma * sqrt(deltat))) for i = 0:N]
        for n = N-1:-1:0
            for i = 0:n
                x = K - S * exp((2*i-n) * sigma * sqrt(deltat))
                y = act * (q * Z[i+1] + p * Z[i+2])
                Z[i+1] = y
            end
        end
        diff = Z[1] - blsprice(S, K, r, t, sigma, dyield)[2]
    elseif class == :call
        Z = [max(0, S*exp((2*i-N) * sigma * sqrt(deltat)) - K) for i = 0:N]
        for n = N-1:-1:0
            for i = 0:n
                x = S * exp((2*i-n) * sigma * sqrt(deltat)) - K
                y = act * (q * Z[i+1] + p * Z[i+2])
                Z[i+1] = y
            end
        end
        diff = Z[1] - blsprice(S, K, r, t, sigma, dyield)[1]
    end
    return Z[1], diff
end

#=
bingreeks computes American option price and the greeks: delta, rho, vega, theta, psi using a binomial tree

    (spotprice, delta, rho, vega, theta, psi) = bingreeks(S, K, r, sigma, t, dyield, class, method)

Uses automatic differentiation (ForwardDiff.jl package) to compute the greeks. The gamma term is exactly zero. This turns out to be a property of the binomial model: as S only ever occurs as a linear term in the code, delta is a discontinuous stepwise function, so gamma (where it is defined) is zero.
=#
function bingreeks(S, K, r, sigma, t, dyield, class, method)
    return binprice(Dual(S, 1,0,0,0,0), K, Dual(r, 0,1,0,0,0), Dual(sigma, 0,0,1,0,0), Dual(t, 0,0,0,1,0), Dual(dyield, 0,0,0,0,1), class, method)
end

#=
mcprice prices Americal put and call options using Monte Carlo simulation methods 

    optionprice = mcprice(S, K, r, sigma, t, dyield, class, method) 

Input arguments:
    class - Option class, possible values are :call and :put
    method - Algorithm used in silulation with method(N,P), where
        N - Number of steps
        P - Number of paths
    E.g. LS(1000,10000) chooses the Longstaff–Schwartz simulation method with 1000 steps and 10000 paths
=#

# Adapted from the source: https://wilmott.com/automatic-for-the-greeks/
function mcprice(S, K, r, sigma, t, dyield, class, method::LS)
    N = method.steps
    P = method.paths
    deltat = t/N
    act = exp(-r * deltat)

    R = exp((r-dyield) * deltat)
    T = typeof(S * act  * exp(-sigma^2 * deltat/2 + sigma * sqrt(deltat)))
    X = Array{T}(undef, N+1,P)
    for p = 1:P
        X[1, p] = x = S
        for n = 1:N
            x *= R * exp(-sigma^2 * deltat / 2 + sigma * sqrt(deltat) * randn())
            X[n+1, p] = x
        end
    end
    if class == :put
        V = [act * max(K - x, 0) for x in X[N+1, :]]
        for n = N-1:-1:1
            I = V .!= 0  # in the money options
            A = [x^d for d = 0:3, x in X[n+1, :]] # design matrix
            β = A[:, I]' \ V[I] # least-square regression
            cV = A' * β         # estimated continuation values
            for p = 1:P
                ev = max(K - X[n+1, p], 0)
                if I[p] && cV[p] < ev
                    V[p] = act * ev
                else
                    V[p] = act * V[p]
                end
            end
        end
        return max(mean(V), K - S)
    elseif class == :call
        V = [act * max(x - K, 0) for x in X[N+1, :]]
        for n = N-1:-1:1
            I = V .!= 0  # in the money options
            A = [x^d for d = 0:3, x in X[n+1, :]] # design matrix
            β = A[:, I]' \ V[I] # least-square regression
            cV = A' * β         # estimated continuation values
            for p = 1:P
                ev = max(X[n+1, p] - K, 0)
                if I[p] && cV[p] < ev
                    V[p] = act * ev
                else
                    V[p] = act * V[p]
                end
            end
        end
        return max(mean(V), S - K)
    end
end

#=
mcgreeks computes American option price and the greeks: delta, rho, vega, theta, psi using Monte Carlo simulation

    (spotprice, delta, rho, vega, theta, psi) = mcgreeks(S, K, r, sigma, t, dyield, class, method)

Uses automatic differentiation (ForwardDiff.jl package) to compute the greeks.
=#
function mcgreeks(S, K, r, sigma, t, dyield, class, method)
    return mcprice(Dual(S, 1,0,0,0,0), K, Dual(r, 0,1,0,0,0), Dual(sigma, 0,0,1,0,0), Dual(t, 0,0,0,1,0), Dual(dyield, 0,0,0,0,1), class, method)
end


