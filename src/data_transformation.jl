#=
Data Transformation (9/9)
OK  arith2geom	    - Arithmetic to geometric moments of asset returns
OK  boxcox	        - Box-Cox transformation
OK  corr2cov	    - Convert standard deviation and correlation to covariance
OK  cov2corr	    - Convert covariance to standard deviation and correlation coefficient
OK  geom2arith	    - Geometric to arithmetic moments of asset returns
OK  holdings2weights    - Portfolio holdings into weights
OK  ret2tick	    - Convert return series to price series
OK  tick2ret	    - Convert price series to return series
OK  weights2holdings    - Portfolio values and weights into holdings
=#

#= holdings2weights converts portfolio holdings into portfolio weights. The weights must satisfy a budget constraint such that the weights sum to Budget for each portfolio.

        weights = holdings2weights(holdings,prices,budget)

Input arguments:
    holdings    - Number of portfolios (NPORTS) by number of assets (NASSETS) matrix with the holdings of NPORTS portfolios containing NASSETS assets.
    prices      - NASSETS vector of asset prices.
    budget      - (Optional) Scalar or NPORTS vector of nonzero budget constraints. Default = 1.

Output arguments:
    weights     - Weights is a NPORTS by NASSETS matrix containing the normalized weights of NPORTS portfolios containing NASSETS assets.

Notes:
    Holdings may be negative to indicate a short position, but the overall portfolio weights must satisfy a nonzero budget constraint.

    The weights in each portfolio sum to the Budget value (which is 1 if Budget is unspecified.)
=#
function holdings2weights(holdings, prices, budget)
    numports = size(holdings,1)
    numassets = size(holdings,2)
    weights = zeros(numports, numassets)
    for i = 1:numports
	    norm = sum(prices .* holdings[i,:]')
	    weights[i,:] = budget[i] .* (prices .* holdings[i,:]') ./ norm
    end
    return weights
end
holdings2weights(holdings, prices) = holdings2weights(holdings, prices, ones(size(holdings)[1]))


#=
weights2holdings converts portfolio values and weights into portfolio holdings.

    holdings = weights2holdings(values,weights,prices)

Input arguments:
    values      - Scalar or number of portfolios (NPORTS) vector containing portfolio values.
    weights     - NPORTS by number of assets (NASSETS) matrix with portfolio weights. The weights sum to the value of a Budget constraint, which is usually 1. (See holdings2weights for information about budget constraints.)
    prices      - NASSETS vector of prices.

Output arguments:
    holdings    -  NPORTS-by-NASSETS matrix containing the holdings of NPORTS portfolios that contain NASSETS assets

    Note: This function does not create round-lot positions. Holdings are floating-point values.
=#
function weights2holdings(values,weights,prices)
    return weights .* values ./ prices
end

#= cov2corr computes the correlation matrix from the covariance matrix `C` and a vector of standard
deviations `s` 
=#
cov2corr = cov2cor

#= cor2cov compute the covariance matrix from the correlation matrix `C` and a vector of standard
deviations `s`
=#
corr2cov = cor2cov



#= tick2ret computes asset returns for NUMOBS price observations of NASSETS assets.

    [returns,intervals] = tick2ret(tickseries, ticktime, method = :simple)

Input arguments:
    data    - Data for asset prices, specified as a matrix, table, or timetable. Prices across a given row are assumed to occur at the same time for all columns, and each column is a price series of an individual asset.
    method  - Method to convert asset prices to returns. :simple (default) or :log. If the method is :simple, then simple periodic returns at time t are computed as:
        Returns(t) = Data(t)/Data(t-1) - 1
    If the method is :log, the continuous returns are computed as:
        Returns(t) = log(Data(t)/Data(t-1))
    Likewise package TimeSeries that defines these two return methods as :simple and :log
Output arguments:
    returns — Time series array of asset returns. returned as a NUMOBS-1-by-NASSETS array of asset returns with the same type (matrix, table, or timetable) as the input Data. The first row contains the oldest returns and the last row contains the most recent. Returns across a given row are assumed to occur at the same time for all columns, and each column is a return series of an individual asset.
    intervals - Interval times between successive prices, returned as a NUMOBS-1 length column vector where Intervals(t) = TickTimes(t) - TickTimes(t-1).
=#
function tick2ret(tickseries::Matrix{Float64}, ticktimes::Matrix{Dates.Date}, method = :simple)
    nobs = size(tickseries,1)-1
    nseries = size(tickseries,2)  
    returns = zeros(nobs,nseries)
    intervals = zeros(nobs) 
    if method == :simple
        for i in 1:nobs
            returns[i,:] = tickseries[i+1,:] ./tickseries[i,:] .- 1
            intervals[i] = Dates.value(ticktimes[i+1]-ticktimes[i])
        end
    elseif method == :log
        for i in 1:nobs
            returns[i,:] = log.(tickseries[i+1,:] ./tickseries[i,:])
            intervals[i] = Dates.value(ticktimes[i+1]-ticktimes[i])
        end
    end
    return returns, intervals * 24
end

function tick2ret(tickseries::Matrix{Float64}, method = :simple)
    nobs = size(tickseries,1)-1
    nseries = size(tickseries,2)  
    returns = zeros(nobs,nseries)
    if method == :simple
        for i in 1:nobs
            returns[i,:] = tickseries[i+1,:] ./tickseries[i,:] .- 1
        end
    elseif method == :log
        for i in 1:nobs
            returns[i,:] = log.(tickseries[i+1,:] ./tickseries[i,:])
        end
    end
    return returns
end

#= ret2tick computes prices from the start prices of NASSET assets and NUMOBS return 

[tickseries,ticktimes] = ret2tick(returnseries, startprice, returnintervals, starttime, method = :simple)

Input arguments:
    startprice      — Initial prices for each asset. Default: 1 for all assets
    returnintervals — Return interval between prices. Default: 1 for all assets
    starttime       - Starting time for first observation applied to the prices of all assets. Default: 0 if returnintervals is numeric  
    method          - Method to convert asset returns to prices. :simple (default) or :log. If the method is  :simple, then simple periodic returns at time t are computed as:
        Returns(t) = Data(t)/Data(t-1) - 1
    If the method is :log, the continuous returns are computed as:
        Returns(t) = log(Data(t)/Data(t-1))
    Likewise package TimeSeries that defines these two return methods as :simple and :log
Output arguments:
    tickseries      — Time series array of asset prices. The first row contains the oldest prices and the last row contains the most recent. Prices across a given row are assumed to occur at the same time for all columns, and each column is a price series of an individual asset.
    ticktimes       - Observation times associated with prices in TickSeries returned as a NUMBOBS+1 length column vector of monotonically increasing observation times associated with the prices in TickSeries.
=#
function ret2tick(returnseries, startprice, returnintervals, starttime::Date, method = :simple)
    nobs = size(returnseries,1)+1
    nseries = size(returnseries,2)
    tickseries = zeros(nobs,nseries)
    ticktimes = Array{Date}(undef, nobs)
    tickseries[1,:] = startprice
    ticktimes[1] = starttime
    if method == :simple
        for i = 2:nobs
            tickseries[i,:] = tickseries[i-1,:] .* (1 .+ returnseries[i-1,:])
            ticktimes[i] = ticktimes[i-1] + Day(returnintervals[i-1])            
        end
    elseif method == :log
        for i = 2:nobs
            tickseries[i,:] = tickseries[i-1,:] .* exp.(returnseries[i-1,:])
            ticktimes[i] = ticktimes[i-1] + Day(returnintervals[i-1])            
        end    
    end
    return tickseries, ticktimes
end

function ret2tick(returnseries, startprice, method = :simple)
    nobs = size(returnseries,1)+1
    nseries = size(returnseries,2)
    tickseries = zeros(nobs,nseries)
    tickseries[1,:] = startprice
    if method == :simple    
        for i = 2:nobs
            tickseries[i,:] = tickseries[i-1,:] .* (1 .+ returnseries[i-1,:])          
        end
    elseif method == :log
        for i = 2:nobs
            tickseries[i,:] = tickseries[i-1,:] .* exp.(returnseries[i-1,:])          
        end
    end
    return tickseries
end

#= boxcox computes the Box-Cox transformation

    transformeddata = boxcox(data, lambda)
    (transformeddata, lambda) = boxcox(data)

    boxcox transforms nonnormally distributed data to a set of data that has approximately normal distribution. The Box-Cox transformation is a family of power transformations.

    If λ is not = 0, then
        data(λ) = (data(λ)−1)/λ
    If λ is = 0, then
        data(λ) = log(data)
    The algorithm calls for finding the λ value that maximizes the Log-Likelihood Function. 
=#
function boxcox(lambda,data::Vector{Float64})
    if lambda == 0.0
        return log.(data)
    else
        return data .^lambda ./ lambda
    end
end 

#=
Arithmetic to geometric moments of asset returns

    (mg, Cg) = arith2geom(ma, Ca, t)

Input arguments:
    ma - Arithmetic mean of asset-return data (n-vector)
    Ca - Arithmetic covariance of asset-return data, an n-by-n symmetric, positive semidefinite matrix
    t - Target period of geometric moments in terms of periodicity of arithmetic moments. Default: 1.   Examples (Matlab): 
        Given arithmetic mean m and covariance C of monthly total returns, obtain annual geometric mean mg and covariance Cg. In this case, the output period (1 year) is 12 times the input period (1 month) so that t = 12. 
        Given arithmetic mean m and covariance C of monthly total returns, obtain quarterly continuously compounded return moments. In this case, the output is 3 of the input periods so that t = 3.
=#
function arith2geom(ma, Ca, t)
    nassets = size(ma,2)  
    mg = zeros(nassets)
    Cg = zeros(nassets,nassets)
    for i in 1:nassets
        mg[i] = exp(t*ma[i] + 1/2*t*Ca[i,i]) - 1
        for j in 1:i
            Cg[i,j] = (1 + mg[i])*(1+mg[j]) * (exp(t * Ca[i,j])-1)
            if j != i
                Cg[j,i] = Cg[i,j]
            end
        end
    end
    return mg, Cg
end
arith2geom(ma, Ca) = arith2geom(ma, Ca, 1)

#=
Geometric to arithmetic moments of asset returns

    (ma, Ca) = geom2arith(mg, Cg, t)

=#
function geom2arith(mg, Cg, t)
    nassets = size(mg,2)  
    ma = zeros(nassets)
    Ca = zeros(nassets,nassets)
    for i in 1:nassets
        for j in 1:i
            Ca[i,j] = t * log(1 + Cg[i,j]/((1+mg[i])*(1+mg[j])))
            if j != i
                Ca[j,i] = Ca[i,j]
            end
        end
        ma[i] = t * log(1+mg[i]) - 1/2*Ca[i,i] 
    end
    return ma, Ca
end
geom2arith(mg, Cg) = geom2arith(mg, Cg, 1)
