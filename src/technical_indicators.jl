#=
Technical Indicators (26/26):

Time Series Oscillations
OK  adosc       - Accumulation/Distribution oscillator
OK  chaikosc    - Chaikin oscillator
OK  macd        - Moving Average Convergence/Divergence (MACD)
OK  stochosc    - Stochastic oscillator
OK  tsaccel     - Acceleration between times
OK  tsmom       - Momentum between times

Volatilities
OK  chaikvolat  - Chaikin volatility
OK  willpctr    - Williams %R               

Volumes
OK  negvolidx	- Negative volume index
OK  posvolidx	- Positive volume index
OK  rsindex     - Relative Strength Index (RSI)

Time Series Rate of Change
OK  adline      - Accumulation/Distribution line
OK  bollinger   - Time series Bollinger band
OK  hhigh	    - Highest high
OK  llow        - Lowest low
OK  medprice	- Median price
OK  movavg      - Moving average of a financial time series
OK  onbalvol    - On-Balance Volume (OBV)
OK  prcroc      - Price rate of change
OK  pvtrend     - Price and Volume Trend (PVT)
OK  typprice    - Typical price
OK  volroc      - Volume rate of change
OK  wclose      - Weighted close
OK  willad      - Williams Accumulation/Distribution line

Trend indicators
OK  Exponential moving average (EMA)
OK  Simple moving average (SMA)

Utility Functions
OK  ret2tick	- Convert return series to price series
OK  tick2ret	- Convert price series to return series
=#

#= 
MarketTechnicals.jl

Moving Averages
    Kaufman's Adaptive Moving Average
    Moving Average Envelope
Momentum
    Average Directional Movement Index (ADX)
    Aroon Oscillator
    CCI - Commodity Channel Index
Volatility
    True Range
    Average True Range
    Donchian Channels
    Keltner Bands
    Volume
    Volume Weight-Adjusted Price (VWAP)
=#

#=
Indicators.jl

Moving Averages
    SMA (simple moving average)
    WMA (weighted moving average)
    EMA (exponential moving average)
    TRIMA (triangular moving average)
    KAMA (Kaufman adaptive moving average)
    MAMA (MESA adaptive moving average, developed by John Ehlers)
    HMA (Hull moving average)
    ALMA (Arnaud-Legoux moving average)
    SWMA (sine-weighted moving average)
    DEMA (double exponential moving average)
    TEMA (triple exponential moving average)
    ZLEMA (zero-lag exponential moving average)
    MMA (modified moving average)
    VWMA (volume-weighted moving average)
    MLR (moving linear regression)
        Prediction
        Slope
        Intercept
        Standard error
        Upper & lower bound
        R-squared

Momentum Indicators
    Momentum (n-day price change)
    ROC (rate of change)
    MACD (moving average convergence-divergence)
    RSI (relative strength index)
    ADX (average directional index)
    Parabolic SAR (stop and reverse)
    Fast & slow stochastics
    SMI (stochastic momentum indicator)
    KST (Know Sure Thing)
    Williams %R
    CCI (commodity channel index)
    Donchian channel
    Aroon indicator + oscillator

Volatility Indicators
    Bollinger Bands
    Average True Range
    Keltner Bands
=#

#=
chaikocs calculates the Chaikin oscillator

    chosc = chaikocs(Data)

Input arguments:
    Data with high, low, open, close information, specified as a matrix, table, or timetable. For matrix input, Data is an M-by-4 matrix of high, low, opening, and closing prices. Timetables and tables with M rows must contain variables named 'High', 'Low', 'Close', and 'Volume' (case insensitive).
=#
function chaikocs(data)
    highp = data[:,1]
    lowp = data[:,2]
    closep = data[:,3]
    tvolume = data[:,4]

    # Calculate the ADLINE and its 10-period and 3-period exponential moving averages.
    ADLine  = adline([highp, lowp, closep, tvolume])
    ma10p = movavg(ADLine,:exponential,10)
    ma03p = movavg(ADLine,:exponential,3)

    # Calculate the Chaikin Oscillator.
    oscillator = ma03p - ma10p
    oscillator[1:9] .= NaN
    return oscillator
end



#= adosc calculates the Accumulation/Distribution (A/D) oscillator.

    ado = adosc(Data)

Input arguments:
    Data — Data with high, low, open, close information specified as a matrix, dataframe or timetable. For matrix input, Data is an M-by-4 matrix of high, low, opening, and closing prices. Timetables and dataframes with M rows must contain variables named 'High', 'Low', 'Open', and 'Close' (case insensitive)
=# 
function adosc(data)
    openp = data[:,1]
    highp = data[:,2]
    lowp = data[:,3]
    closep = data[:,4]

    oscillator = 100 * ((highp-openp) + (closep-lowp)) / (2*(highp-lowp))
    return oscilator
end


#=
adline calculates the Accumulation/Distribution Line from a set of high, low, closing prices and volume traded of a security

    ADline = adline(Data)

Input arguments
    Data - Data with high, low, closing prices and volume traded, specified as a matrix, table, or timetable. For matrix input, Data is M-by-4 with high, low, closing prices, and volume traded. Timetables and tables with M rows must contain a variable named 'High', 'Low', 'Close', and 'Volume' (case insensitive)
Output argument:
    ADline - Accumulation/Distribution line, returned with the same number of rows (M) and the same type (matrix, table, or timetable) as the input data.
=#
function adline(data)
    highp = data[:,1]
    lowp = data[:,2]
    closep = data[:,3]
    tvolume = data[:,4]

    tADline = tvolume * ((closep - lowp) - (highp - closep)) / (highp - lowp)
    return tADline
end


#=
tsaccel calculates the acceleration of a data series with time distance of n periods.  It is essentially the difference of the current momentum with the momentum n periods ago. By default, Acceleration is based on 12-period difference.

    acceleration = tsaccel(Data, NumPeriods, Datatype)

Acceleration is defined as the difference of two momentum series separated by N periods.

Input arguments:
    Data        - A vector, matrix, table, or timetable. For vector input, Data is column vector. For matrix input, Data is an M-by-N column-oriented matrix.
    NumPeriods - Integer indicating the period difference for acceleration. The default is 12.
    Datatype   - Integer scalar indicating whether Data contains the data itself or the momentum of the data. Valid choices are: 0 (Data contains the data itself), 1 (Data contains the momentum of the data). The default is 0 (Data contains the data itself).
Output arguments:
    acceleration - Acceleration series with the same number of rows (M), columns (N), and type as the input data.
=#
function tsaccel(data, numPeriods, datatype)
    momentum = data[numPeriods:end,:] - data[1:(end-numPeriods+1),:]
    momentum = [NaN(numPeriods-1,numVar);momentum] # Fill in NaNs in front.
    if datatype == 0
        acceleration = momentum[numPeriods:end,:] - momentum[1:(end-numPeriods+1),:]
        acceleration = [NaN(numPeriods-1,numVar);acceleration] # Fill in NaNs in front.
    elseif datatype == 1
        acceleration = momentum
    end
end


#=
tsmom calculates the momentum of a data series with time distance of n periods. It is essentially the difference of the current data with the data n periods ago. By default, Momentum is based on 12-period difference.

    momentum = tsmom(Data,NumPeriods)

Input arguments:
    Data        - A vector, matrix, table, or timetable. 
    NumPeriods  - Integer indicating the period difference for Momentum. The default is 12.
Output argument:
    momentum - Momentum series with the same number of rows, columns and type as the input data.
=#
function tsmom(data,numPeriods)
    # Calculate the Momentum.
    momentum = data[numPeriods:end,:] - data[1:(end-numPeriods+1),:]
    momentum = [NaN(numPeriods-1,numVar);momentum] # Fill in NaNs in front.
end


#=
willad calculates the Williams Accumulation/Distribution Line from the series of high, low, and closing prices.

    WADLine = willad(Data)

Input arguments:
    Data    - A matrix, table, or timetable. For matrix input, Data is an M-by-3 matrix of high, low, and closing prices stored in the corresponding columns, respectively. Timetables and tables with M rows contain variables named 'High', 'Low', and 'Close' (case insensitive).
=#
function willad(data)
    highp = data(:,1)
    lowp = data(:,2)
    closep = data(:,3)

    TRH = max(highp[2:end],closep[1:end-1])
    TRL = min(lowp[2:end],closep[1:end-1])
    dTRL = closep[2:end] - TRL
    dTRH = closep[2:end] - TRH
    dCLOSE = diff(closep);
    
    # Calculate Today's A/D.
    todaysad = zeros(length(TRH),1)
    if dCLOSE > 0
        todaysad = dTRL
    elseif dCLOSE < 0 
        todaysad = dTRH
    end
    
    # Calculate the Williams A/D Line.
    WADLine = cumsum([closep(1);todaysad])
    return WADLine
end

#=
ema computes an exponential moving average (EMA)

    ema(data, windowSize)

The moving average is initialised such that the initial points include only observed data (em movavg()opção :shrink)
=#
function ema(data, windowSize::Int64)
    movavg(data, :exponential, windowSize)
end


#=
sma computes a simple moving average

    ema(data, windowSize)

TThe moving average is initialised such that the initial points include only observed data (em movavg()opção :shrink)
=#
function sma(data, windowSize::Int64)
    movavg(data, :simple, windowSize)
end
#=
macd computes the Moving Average Convergence/Divergence (MACD) indicator as well as the 9-period exponential moving average from the 'macd' line

    (MACDLine,SignalLine) = macd(sata)


The MACD is calculated by subtracting the 26-period (7.5%) exponential  moving average from the 12-period (15%) moving average.  The 9-period (20%) exponential moving average of the MACD line is used as the "signal" line.

When the two lines are plotted, they can give you indications on when to buy or sell a stock, when overbought or oversold is occurring, and when the end of trend may occur. For example, when the MACD and the 20 moving average line has just crossed and the MACD line becomes below the other line, it is time to sell.
=#
function macd(data)
    # Calculate the 26-period (7.5%)exp mov avg and the 12-period (15%) exp mov avg
    # 26-period (7.5%) exp mov avg
    ema26p = movavg(data, :exponential, 26)
    # 12-period (15%) exp mov avg
    ema12p = movavg(data, :exponential, 12)

    # Calculate the 9-period (20%) exp mov avg of the MACD line.
    MACDLine = ema12p - ema26p;
    SignalLine = movavg(MACDLine, :exponential, 9)
    MACDLine[1:25,:] .= NaN
    SignalLine[1:34,:] .= NaN;
end

#=
onbalvol calculates the On-Balance Volume from the series of closing stock prices and trade volume.

volume = onbalvol(data)

Input arguments:
    data - A matrix, table, or timetable. For matrix input, Data is an M-by-2 matrix of closing prices and trade volume stored in the first and second columns, respectively. Timetables and tables with M rows contain variables named 'Close', 'Volume' (case insensitive).
=#
function onbalvol(data)
    closep = data(:,1)
    tvolume = data(:,2)
    volume = tvolume;
    for didx = 2:length(tvolume)
        if closep[didx] > closep[didx-1]
            volume[didx] = volume[didx-1] + tvolume[didx]
        elseif closep[didx] < closep[didx-1]
            volume[didx] = volume[didx-1] - tvolume[didx]
        elseif closep[didx] == closep[didx-1]
            volume[didx] = volume[didx-1]
        end
    end
    return volume
end

#=
medprice calculates the median prices from the series of high and low prices. The median price is just the average of the high and low prices for each period.

    medprice(data)

=#
function medprice(data)
    highp = data[:,1]
    lowp = data[:,2]
    return (highp + lowp) / 2
end


#=
chaikvolat calculates the Chaikin's volatility from the series of high and low stock prices. By default, Chaikin's volatility values are based on a 10-period exponential moving average and 10-period difference. 

    volatility = chaikvolat(Data, NumPeriods, WindowSize)

Input arguments:
    Data       - A matrix, table, or timetable. For matrix input, Data is an M-by-2 matrix of high and low prices stored in the first and second columns, respectively. Timetables and tables with M rows contain variables named 'High', 'Low' (case insensitive).
    NumPeriods - Positive integer scalar indicating the period difference. The default is 10.
    WindowSize - Positive integer scalar indicating the length of the exponential moving average in periods. The default is 10.
=#
function chaikvolat(data, numperiods = 10, windowSize = 10)
    highp = data[:,1]
    lowp = data[:,2]

    hlEMA = movavg(highp-lowp, :exponential, windowSize)
    for i in 1:windowSize-1
        hlEMA[i] = NaN
    end
    # Calculate the Chakin's volatility based on difference period difference.             
    volatility = (hlEMA[numPeriods:end] - hlEMA[1:end-numPeriods+1]) ./ (hlEMA[1:end-numPeriods+1])*100
end

#=
bollinger calculates the Bollinger bands.

    [middle, upper, lower] = bollinger(data, windowSize, Type, NumStd)

Input arguments:
    WindowSize - Integer indicating the number of observations of the input series to include in the moving average in periods. The default is 10.
    Type - Indicator of the particular type of moving average to compute. Valid choices are :simple, and :linear. The default is :simple.
    NumStd  - Number of standard deviations for the upper and lower bands. The default is 2.
=#
function bollinger(data, windowSize, type, numstd)
    middle = movavg(data, type, windowSize, :fill)
    mstd = movstd(data, windowSize, :fill)
    upper = middle + numstd * mstd
    lower = middle - numstd * mstd
    return middle, upper, lower
end


#=
stochosc calculates the Fast PercentK (F%K), Fast PercentD (F%D), Slow PercentK (S%K), and Slow PercentD (S%D) from the series of high, low, and closing stock prices. By default, Stochastic Oscillator is based on 10-period difference for PercentK and a 3-period exponential moving average for PercentD. 

    percentKnD = stochosc(data, numPeriodsK, numPeriodsD, type)

Input arguments:
    data    - For matrix input, Data is an M-by-3 matrix of high, low, and closing prices stored in the 
            corresponding columns, respectively. Timetables and tables with M rows contain variables named 'High', 'Low', and 'Close' (case insensitive).
    NumPeriodsK - Indicating of the period difference for PercentK. The default is 10.
    NumPeriodsD - Positive integer scalar indicating the length of moving average in periods for PercentD. The default is 3.
    Type        - Indicating of the moving average methods for percentD calculation. Valid options are:
                    :exponential (default)
                    :triangular
Output argument:
   percentKnD   - PercentK and PercentD 
=#
function stochosc(data, numPeriodsK, numPeriodsD, type)
    highp = data[:,1]
    lowp = data[:,2]
    closep = data[:, 3]
    llv = llow(lowp, numPeriodsK)
    hhv = hhigh(highp, numPeriodsK)

    # Calculate the PercentK (%K).
    pctk = 100 * (closep .- llv)./(hhv .- llv)
    # Calculate the PercentD (%D).
    pctd = movavg(pctk, type, numPeriodsD)
    # Calculate the Slow PercentK (%K).
    spctk = pctd;
    spctd = movavg(spctk, type, numPeriodsD)
    # Populate NaNs
    pctd[1:numPeriodsD-1] .= NaN
    spctk[1:numPeriodsD-1] .= NaN
    spctd[1:2*numPeriodsD-1] = NaN;

    return pctk, pctd, spctk, spctd
end


#=
rsindex calculates the Relative Strength Index (RSI) from the series of closing stock prices. By default, RSI values are based on a 14-period window.

    index = rsindex(data, windowSize)

=#
function rsindex(Data,WindowSize)
    closep = data
    numObs = length(closep)
    # Take a diff of the closing prices
    priceChange = diff(closep)

    # Create '+' Delta vectors and '-' Delta vectors
    advances = priceChange
    declines = priceChange

    for i in 1:numObs
        if priceChange[i] < 0
            advances[i] = 0
        elseif priceChange[i] >= 0
            declines[i] = 0
        end
    end
    declines = -declines

    totalGain = movsum(advances,windowSize, :discard)
    totalLoss = movsum(declines,windowSize, :discard)

    # Calculate RSI
    rs = totalGain ./ totalLoss
    index = 100 .- (100 ./ (1 .+ rs))
    for i in 1:windowSize
        index[i] = NaN
    end
end

#=
willpctr calculates the Williams PercentR (%R) values for a data series of
with high, low, and closing prices. By default, Williams PercentR values are based on 14 periods.

    PercentR = willpctr(data, numPeriods)

=#
function  willpctr(data, numPeriods = 14)
    highp = data[:,1]
    lowp = data[:,2]
    closep   = data[:, 3]

    llv = llow(lowp, numPeriods)
    hhv = hhigh(highp, numPeriods)

    # Calculate the Williams PercentR.
    PercentR = -100 * (hhv .- closep)./(hhv .- llv)
end

#=
pvtrend calculates the Price and Volume Trend (PVT) from the series of closing stock prices and trade volume.

    trend = pvtrend(data)

=#
function pvtrend(data)
    closep = data[:,1]
    tvolume = data[:,2]
    trend = [tvolume[1]; tvolume[2:end] .* diff(closep) ./ closep[1:end-1]]
    trend = cumsum(trend)
    return trend
end

#=
typprice calculates the typical prices from the series of high, low, and closing prices. The typical price is just the average of the high, low, and closing prices for each period.

    TypicalPrice = typprice(data)

=#
function typprice(data)
    highp = data[:,1]
    lowp = data[:,2]
    closep = data[:,3]
    return (highp + lowp + closep) / 3
end



#=
posvolidx calculates the positive volume index from the series of closing stock prices and trade volume. By default, the initial values for the negative volume index is set to 100.

    volume = posvolidx(data,initialvalue)

Input arguments:
    InitialValue - Scalar indicating the initial value for positive volume index. The default is 100.
=#
function posvolidx(data,initialvalue)
    closep = data[:,1]
    tvolume = data[:,2]

    volume = initialvalue * ones(length(tvolume),1)

    for didx in 2:length(tvolume)
        if tvolume[didx] > tvolume[didx-1]
            volume[didx] = volume[didx-1] * (1 + ((closep[didx]-closep[didx-1])./closep[didx-1]))
        elseif tvolume[didx] <= tvolume[didx-1]
            volume[didx] = volume[didx-1]
        end
    end
    return volume
end


#=
negvolidx calculates the negative volume index from the series of closing stock prices and trade volume. By default, the initial values for the negative volume index is set to 100.

    volume = negvolidx(data,initialvalue)

Input arguments:
    InitialValue - Scalar indicating the initial value for negative volume index. The default is 100.
=#
function negvolidx(data,initialvalue)
    closep = data[:,1]
    tvolume = data[:,2]

    volume = initialvalue * ones(length(tvolume),1)

    for didx in 2:length(tvolume)
        if tvolume[didx] < tvolume[didx-1]
            volume[didx] = volume[didx-1] * (1 + ((closep[didx]-closep[didx-1])./closep[didx-1]))
        elseif tvolume[didx] >= tvolume[didx-1]
            volume[didx] = volume[didx-1]
        end
    end
    return volume
end



#=
wclose calculates the weighted closing prices from the series of high, low, and closing prices. The weighted closing price is the average of twice the closing price plus the high and low prices.

     WeightedClose = wclose(data)

=#
function wclose(data)
    highp = data[:,1]
    lowp = data[:,2]
    closep = data[:,3]
    return (2*closep + highp + lowp) / 4
end

#=
proc calculates the price rate-of-change, PriceChangeRate, from the series of closing stock prices. By default, the price rate-of-change is calculated between the current closing price and the closing price 
12 periods ago.

    PriceChangeRate = prcroc(data,numPeriods)

=#
function prcroc(data, numPeriods = 12)
    closep = data
    PriceChangeRate = (closep[numPeriods:end] - closep[1:end-numPeriods+1])./closep[1:end-numPeriods+1]*100
    for i in 1:numPeriods-1
        PriceChangeRate[i] = NaN
    end
    return PriceChangeRate
end


#=
volroc calculates the volume rate-of-change from a data series of volume traded. The volume rate-of-change is calculated between the current volume and the volume n periods ago. By default, Volume Rate-of-Change is based on 12-period difference.

    volumeChangeRate = volroc(data,numPeriods)

=#
function volroc(data,numPeriods = 12)
    tvolume = data
    volumeChangeRate = 100 * (tvolume[numPeriods:end] - tvolume[1:(end-numPeriods+1)]) ./ tvolume[1:(end-numPeriods+1)]
    for i in 1:numPeriods-1
        volumeChangeRate[i] = NaN
    end
    return volumeChangeRate[i]
end


#= 
hhigh compute the rolling high of a financial series high prices

    hhigh(A, windowsize = 14, endpoints::Symbol = :shrink)

=#
hhigh(A, windowsize = 14, endpoints::Symbol = :shrink) = movmax(A, windowsize, endpoints)


#= 
llow compute the rolling low of a financial series low prices

    llow(A, windowsize = 14, endpoints::Symbol = :shrink)

=#
llow(A, windowsize = 14, endpoints::Symbol = :shrink) = movmin(A, windowsize, endpoints)


