#= 
Stats (10/10):

OK  movavg  - Compute the moving average (MA) of a financial time series
OK  moving  - Compute the rolling window of a financial time series for a general function
OK  movmin  - Compute the rollwing low of a time series
OK  movmax  - Compute the rollwing high of a time series
OK  movsum  - Compute the rollowing sum of a financial time series
OK  movmean - Compute the rollwing simple mean of a financial time series
OK  movmedian   - Compute the rolling correlation for a time series
OK  movstd  - Compute the rollwing std of a financial time series
OK  movmad  - Compute the rolling median absolute deviation (mad) around the median for a time series
OK  movcor  - Compute the rolling correlation between two financial time series
OK  movslope    - Compute the rolling slope of a regression
=#



# helper functions: simpleMA, weightGenerate and weightMA

function simpleMA(data, windowSize, endpoints)
    # Simple moving average
    # Replicate existing behavior
    weights = weightGenerate(windowSize, 0)
    ma = weightMA(data, windowSize, weights, endpoints)
    return ma
end
    
function weightGenerate(windowSize, alpha)
    i = 1:windowSize
    weights = (windowSize .- i .+ 1) .^ alpha ./ sum((1:windowSize) .^ alpha)
    return weights
end
    
function weightMA(data, windowSize, weights, endpoints)
    # build moving average vectors by filtering asset through weights
    weightedMA = filt(weights, 1, data)         # filt() is from the DSP.jl package
    if endpoints == :zero
        return weightedMA
    elseif endpoints == :fill
        weightedMA = Array{Union{Float64,Missing}}(weightedMA)
        windowSizeL1 = convert(Int64, windowSize-1)
        weightedMA[1:windowSizeL1] .= missing      # replacing missing with NaN gets a +30% speed improvement
        return weightedMA
    elseif endpoints == :fillNaN
        windowSizeL1 = convert(Int64, windowSize-1)
        weightedMA[1:windowSizeL1] .= NaN      # replacing missing with NaN gets a +30% speed improvement
        return weightedMA
    elseif endpoints == :shrink
        # calculate first n points.
        weightsSum = cumsum(weights)
        windowSize = convert(Int64, windowSize)
        weightedMA[1:windowSize] = weightedMA[1:windowSize] ./ weightsSum
        return weightedMA
    end
end


#=
movavg computes the moving average (MA) of a financial time series.

    ma = movavg(data,type,windowSize, weights, endpoints) 

Input arguments:
    type - Type of moving average to compute. Possible values are: 
        :simple 
        :squareroot 
        :linear
        :square 
        :exponential
        :triangular
        :modified
        :custom
    windowSize - Number of observations of the input series to include in moving average
    weights - Custom weights used to compute the moving average. The weights argument applies only to a :custom type of moving average.
    endpoints - Indicates how moving average is calculated at initial points. (before there is enough data to fill the window). Possible values are: 
        :shrink (default) - Initializes the moving average such that the initial points include only observed data 
        :fill - Fills initial points with missing
        :fillNaN    - Fills initial points with NaN
        :zero - Initializes the initial points with 0
=#
function movavg(data, type::Symbol, windowSize::Int64, endpoints::Symbol = :shrink)
    if type == :simple
        ma = simpleMA(data, windowSize, endpoints)
        return ma
    
    elseif type == :squareroot
        weights = weightGenerate(windowSize, 0.5);
        ma = weightMA(data,windowSize, weights, endpoints)
        return ma
    
    elseif type == :linear
        weights = weightGenerate(windowSize,1);
        ma = weightMA(data, windowSize, weights, endpoints)
        return ma
    
    elseif type == :square
        weights = weightGenerate(windowSize,2)
        ma = weightMA(data, windowSize, weights, endpoints)
        return ma

    elseif type == :triangular
        # Moving average of a simple moving average.
        doubleSmooth = ceil((windowSize+1)/2)
        movAvg1 = simpleMA(data, doubleSmooth, endpoints)
        ma = simpleMA(movAvg1, doubleSmooth, endpoints)
        return ma

    elseif type == :exponential
        # Calculate the exponential percentage
        alpha = 2/(windowSize+1)
        # Formula:
        # y_0 = x_0
        # y_i = alpha * x_i + (1 - alpha) * y_(i-1)
        ma = filt(alpha, [1,(alpha-1)], data[2:end,:], (1-alpha) * data[1,:]) # filt() is from DSP
        ma = [data[1,:] ;ma]
        return ma

    elseif type == :modified
        # The first point of the modified moving average is calculated the same way the first point of the simple moving average is calculated. However, all subsequent points are calculated using the modified mov avg formula.
        # Formula:
        # y_0 = x_0
        # y_i = (x_i - y_(i-1))/windowSize + y_(i-1)
        alpha = 1/windowSize
        ma = filt(alpha, [1,(alpha-1)], data[2:end,:], (1-alpha) * data[1,:]) # filt() is from DSP
        ma = [data[1,:]; ma]
        return ma
    end
end

function movavg(data, type::Symbol, weights::Vector{Float64}, endpoints::Symbol = :shrink)
    if type == :custom
        normalizedWeights = weights / sum(weights)
        windowSize = length(weights)
        ma = weightMA(data, windowSize, normalizedWeights, endpoints)
        return ma
    end
end

#=
moving computes a moving window for a general function f with one (f(A)) or two datasets arguments (f(A,B))

    moving(f, data, windowsize, endpoints)  

Input arguments:
    f           - Function for wich to compute the moving window. Eg: minimum, maximum, sum
    data        - Data vector
    windowsize  - windowsize can be one of the following
        [nbefore, nafter] - The calculation includes the element in the current position, 'nbefore' elements backward, and 'nafter' elements forward. Example: moving(minimum, A, [1 1], endpoints) centered moving minimum with window size equal to 3.
        nobs - Number of observations of the input series to include in moving window. Example: moving(minimum, A, 3, endpoints), equivalent to moving(minimum, A, [2 0], endpoints)
    endpoints   - controls how the sum is calculated at the endpoints of X, where there are not enough elements to fill the window. endpoints can be one of the following:
        :shrink - (default) compute the sum over the number of elements of X that are inside the window, effectively reducing the window size to fit X at the endpoints. The window size is automatically truncated at the endpoints when there are not enough elements to fill the window. When the window is truncated, the standard deviation is taken over only the elements that fill the window.
        :fill    - compute the sum over the full window size, filling missing values from X with 'missing'. This is equivalent to padding X with 'missing' at the endpoints.
        :fillNaN    - compute the sum over the full window size, filling missing values from X with 'NaN'. This is equivalent to padding X with 'NaN' at the endpoints.
        :discard   - compute the sum only when the window is filled with elements of X, discarding partial endpoint calculations and their corresponding elements in Y. This truncates the output; for a vector X and window length K, Y has length LENGTH(X)-K+1.
=#
function moving(f::Function, A, windowsize, endpoints::Symbol = :shrink) # f applied to one dataset A
    if length(windowsize) == 1
        nbefore = windowsize-1
        nafter = 0
    elseif length(windowsize) == 2
        nbefore = windowsize[1]
        nafter = windowsize[2]
    end
    if endpoints == :shrink
        len = length(A)
        resp = zeros(len)
        for i in 1:nbefore
            resp[i] = f(A[1:i+nafter])
        end
        for i in (nbefore+1):(len-nafter)
            resp[i] = f(A[i-nbefore:i+nafter])
        end
        for i in (len-nafter+1):len
            resp[i] = f(A[i-nbefore:end])
        end
        return resp    
    elseif endpoints == :fill
        len = length(A)
        resp = allowmissing(zeros(len))
        for i in 1:nbefore
            resp[i] = missing
        end
        for i in (nbefore+1):(len-nafter)
            resp[i] = f(A[i-nbefore:i+nafter])
        end
        for i in (len-nafter+1):len
            resp[i] = missing
        end
        return resp    
    elseif endpoints == :fillNaN
        len = length(A)
        resp = zeros(len)
        for i in 1:nbefore
            resp[i] = NaN
        end
        for i in (nbefore+1):(len-nafter)
            resp[i] = f(A[i-nbefore:i+nafter])
        end
        for i in (len-nafter+1):len
            resp[i] = NaN
        end
        return resp    
    elseif endpoints == :discard
        len = length(A)-nbefore-nafter
        resp = zeros(len)
        for i in 1:len
            resp[i] = f(A[i:i+nbefore+nafter])
        end
        return resp    
    end
end


function moving(f::Function, A, B, windowsize, endpoints::Symbol = :shrink) # # f applied to two datasets A, B
    if length(windowsize) == 1
        nbefore = windowsize-1
        nafter = 0
    elseif length(windowsize) == 2
        nbefore = windowsize[1]
        nafter = windowsize[2]
    end

    if endpoints == :shrink
        len = length(A)
        resp = zeros(len)
        for i in 1:nbefore
            resp[i] = f(A[1:i+nafter], B[1:i+nafter])
        end
        for i in (nbefore+1):(len-nafter)
            resp[i] = f(A[i-nbefore:i+nafter], B[i-nbefore:i+nafter])
        end
        for i in (len-nafter+1):len
            resp[i] = f(A[i-nbefore:end], B[i-nbefore:end])
        end
        return resp    
    elseif endpoints == :fill
        len = length(A)
        resp = allowmissing(zeros(len))
        for i in 1:nbefore
            resp[i] = missing
        end
        for i in (nbefore+1):(len-nafter)
            resp[i] = f(A[i-nbefore:i+nafter], B[i-nbefore:i+nafter])
        end
        for i in (len-nafter+1):len
            resp[i] = missing
        end
        return resp    
    elseif endpoints == :fillNaN
        len = length(A)
        resp = zeros(len)
        for i in 1:nbefore
            resp[i] = NaN
        end
        for i in (nbefore+1):(len-nafter)
            resp[i] = f(A[i-nbefore:i+nafter], B[i-nbefore:i+nafter])
        end
        for i in (len-nafter+1):len
            resp[i] = NaN
        end
        return resp    
    elseif endpoints == :discard
        len = length(A)-nbefore-nafter
        resp = zeros(len)
        for i in 1:len
            resp[i] = f(A[i:i+nbefore+nafter],B[i:i+nbefore+nafter])
        end
        return resp    
    end
end

#= 
movmin compute the rolling minimum of a time series

    movmin(A, windowsize, endpoints::Symbol = :shrink)

=#
movmin(A, windowsize, endpoints::Symbol = :shrink) = moving(minimum, A, windowsize, endpoints)


#= 
movmax compute the rolling maximum of a time series high prices

    movmax(A, windowsize, endpoints::Symbol = :shrink)

=#
movmax(A, windowsize, endpoints::Symbol = :shrink) = moving(maximum, A, windowsize, endpoints)


#= 
movstd compute the rolling std of a time series

    movstd(A, windowsize, endpoints)

Input arguments:
    endpoints - controls how the sum is calculated at the endpoints of X, where there are not enough elements to fill the window. endpoints can be one of the following:
        :shrink - (default) compute the sum over the number of elements of X that are inside the window, effectively reducing the window size to fit X at the endpoints. The window size is automatically truncated at the endpoints when there are not enough elements to fill the window. When the window is truncated, the standard deviation is taken over only the elements that fill the window.  
        :fill    - compute the sum over the full window size, filling missing values from X with NaN. This is equivalent to padding X with NaN at the endpoints.
        :discard   - compute the sum only when the window is filled with elements of X, discarding partial endpoint calculations and their corresponding elements in Y. This truncates the output; for a vector X and window length K, Y has length LENGTH(X)-K+1.
=#
movstd(A, windowsize, endpoints::Symbol = :shrink) = moving(std, A, windowsize, endpoints)


#= 
movsum compute the rolling sum of a time series

    movsum(A, windowsize, endpoints)

Input arguments:
    endpoints - controls how the sum is calculated at the endpoints of X, where there are not enough elements to fill the window. endpoints can be one of the following:
        :shrink - (default) compute the sum over the number of elements of X that are inside the window, effectively reducing the window size to fit X at the endpoints. The window size is automatically truncated at the endpoints when there are not enough elements to fill the window. When the window is truncated, the standard deviation is taken over only the elements that fill the window.
        :fill    - compute the sum over the full window size, filling missing values from X with NaN. This is equivalent to padding X with NaN at the endpoints.
        :discard   - compute the sum only when the window is filled with elements of X, discarding partial endpoint calculations and their corresponding elements in Y. This truncates the output; for a vector X and window length K, Y has length LENGTH(X)-K+1.
=#
movsum(A, windowsize, endpoints::Symbol = :shrink) = moving(sum, A, windowsize, endpoints)


#= 
movmean compute the rolling mean for a time series

    movmean(A, windowsize, endpoints)

Input arguments:
    endpoints - controls how the sum is calculated at the endpoints of X, where there are not enough elements to fill the window. endpoints can be one of the following:
        :shrink - (default) compute the sum over the number of elements of X that are inside the window, effectively reducing the window size to fit X at the endpoints. The window size is automatically truncated at the endpoints when there are not enough elements to fill the window. When the window is truncated, the standard deviation is taken over only the elements that fill the window.
        :fill    - compute the sum over the full window size, filling missing values from X with NaN. This is equivalent to padding X with NaN at the endpoints.
        :discard   - compute the sum only when the window is filled with elements of X, discarding partial endpoint calculations and their corresponding elements in Y. This truncates the output; for a vector X and window length K, Y has length LENGTH(X)-K+1.
=#
movmean(A, windowsize, endpoints::Symbol = :shrink) = moving(mean, A, windowsize, endpoints)

#= 
movmedian compute the rolling median for a time series

    movmedian(A, B, windowsize, endpoints)

=#
movmedian(A, windowsize, endpoints::Symbol = :shrink) = moving(median, A, windowsize, endpoints)


#= 
movmad compute the rolling median absolute deviation (mad) around the median for a time series

    movmad(A, B, windowsize, endpoints)

Note: needs the StatsBase.jl package
=#
movmad(A, windowsize, endpoints::Symbol = :shrink) = moving(mad, A, windowsize, endpoints)



#= 
movcor compute the rolling correlation between two time series

    movcor(A, B, windowsize, endpoints)

=#
movcor(A, B, windowsize, endpoints::Symbol = :shrink) = moving(cor, A, B, windowsize, endpoints)


#= 
movslope compute the rolling slope of a regression of a response time series (resp) on a explanatory time series (expln)

    movbeta(resp, expln, windowsize, endpoints)

=#
slope(resp, expln) = cov(resp,expln)/var(expln) 
movslope(resp, expln, windowsize, endpoints::Symbol = :shrink) = moving(slope, resp, expln, windowsize, endpoints)
    