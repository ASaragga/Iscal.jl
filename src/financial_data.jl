
#=
Setting data providers keys. Example: 
    setkey(:AlphaVantage)
    st = stocks(:IBM)
=#
function setkey(datasource)
    if datasource == :AlphaVantage
        AlphaVantage.global_key!(Key[datasource])
    end
end


# helper functions to get raw data into a DataFrame.
function stocks_to_dataframe(rawData)
    df = DataFrame(rawData[1])
    dfNames = Symbol.(replace.(vcat(rawData[2]...), " " => "_")) # For some intervals we get "adjusted close" and "dividend amount"
    #println(nrow(df),"x",ncol(df), " DataFrame with column names ", dfNames)
    df = DataFrames.rename(df, dfNames)

    df.timestamp = Date.(df.timestamp)
    for x in (:open, :high, :low, :close, :adjusted_close, :dividend_amount)
        df[!, x] = Float64.(df[!, x])
    end 
    df.volume = Int64.(df.volume)

    return df
end

function stocks_intra_to_dataframe(rawData)
    df = DataFrame(rawData[1])
    dfNames = Symbol.(replace.(vcat(rawData[2]...), " " => "_")) # For some intervals we get "adjusted close" and "dividend amount"
    #println(nrow(df),"x",ncol(df), " DataFrame with column names ", dfNames)
    df = DataFrames.rename(df, dfNames)

    df.timestamp = DateTime.(df.timestamp, "yyyy-mm-dd HH:MM:SS")
    for x in (:open, :high, :low, :close)
        df[!, x] = Float64.(df[!, x])
    end 
    df.volume = Int64.(df.volume)
    return df
end


function forex_to_dataframes(rawData)
    df = DataFrame(rawData[1])
    dfNames = Symbol.(replace.(vcat(rawData[2]...), " " => "_")) # For some intervals we get "adjusted close" and "dividend amount"
    #println(nrow(df),"x",ncol(df), " DataFrame with column names ", dfNames)
    df = DataFrames.rename(df, dfNames)

    df.timestamp = Date.(df.timestamp)
    for x in [:open, :high, :low, :close]
        df[!, x] = Float64.(df[!, x])
    end
    return df
end

function forex_intra_to_dataframes(rawData)
    df = DataFrame(rawData[1])
    dfNames = Symbol.(replace.(vcat(rawData[2]...), " " => "_")) # For some intervals we get "adjusted close" and "dividend amount"
    #println(nrow(df),"x",ncol(df), " DataFrame with column names ", dfNames)
    df = DataFrames.rename(df, dfNames)

    df.timestamp = DateTime.(df.timestamp, "yyyy-mm-dd HH:MM:SS")
    for x in (:open, :high, :low, :close)
        df[!, x] = Float64.(df[!, x])
    end 
    return df
end


#=
Import stock prices into a dataframe

    dataframe = alphavantage(symbol,freq, interval)

Input arguments:
    symbol      - Stock symbol
    interval    - Sampling interval for time series. Possible values are:
        "1m"  -  1 minute  (or :m)
        "5m"  -  5 minutes (or :th, twelfth hour) 
        "15" - 15 minutes (or :qh, quarter hour)
        "30m" - 30 minutes (or :hh, half hour)
        "h"   - 1 hour
        "d"   - 1 day
        "wk"   - 1 week
        "mo"  - 1 month
=#
function alphavantage(symbol::Symbol, interval::String = "1d")
    if interval == "1d"
        raw = AlphaVantage.time_series_daily_adjusted(String(symbol), outputsize="full", datatype="csv")
        df = stocks_to_dataframe(raw)
    elseif interval == "1wk"
        raw = AlphaVantage.time_series_weekly_adjusted(String(symbol), outputsize="full", datatype="csv")
        df = stocks_to_dataframe(raw)
    elseif interval == "1mo"
        raw = AlphaVantage.time_series_monthly_adjusted(String(symbol), outputsize="full", datatype="csv")
        df = stocks_to_dataframe(raw)
    elseif interval == "1m" || interval == "5m" || interval == "15m" || interval == "30m" || interval == "1h"
        if interval == "1m"
            sinterval = "1min"
        elseif interval == "5m"
            sinterval = "5min"
        elseif interval == "15m"
            sinterval = "15min"
        elseif interval == "30m"
            sinterval = "30min"
        elseif interval == "1h"
            sinterval = "60min"
        end
         raw = AlphaVantage.time_series_intraday(String(symbol), sinterval, outputsize="full", datatype="csv")
         df = stocks_intra_to_dataframe(raw)
    end
    return df
end



#=
Import foreign exchange rates into a dataframe 

    dataframe = alphavanatge(symbol1, symbol2, freq, interval)

Input arguments:
    symbol1     - Currency symbol
    symbol2     - Currency symbol
    interval    - Sampling interval for time series. Possible values are:
        "1m"  -  1 minute  (or :m)
        "5m"  -  5 minutes (or :th, twelfth hour) 
        "15" - 15 minutes (or :qh, quarter hour)
        "30m" - 30 minutes (or :hh, half hour)
        "h"   - 1 hour
        "d"   - 1 day
        "wk"   - 1 week
        "mo"  - 1 month
=#
function alphavantage(symbol1::Symbol, symbol2::Symbol, interval::String = "1d")
    if interval == "1d"
        raw = AlphaVantage.fx_daily(String(symbol1), String(symbol2), datatype="csv")
        df = forex_to_dataframes(raw)
    elseif interval == "1wk"
        raw = AlphaVantage.fx_weekly(String(symbol1), String(symbol2), datatype="csv")
        df = forex_to_dataframes(raw)
    elseif interval == "1mo"
        raw = AlphaVantage.fx_monthly(String(symbol1), String(symbol2), datatype="csv")
        df = forex_to_dataframes(raw)
    elseif interval == "1m" || interval == "5m" || interval == "15m" || interval == "30m" || interval == "1h"
        if interval == "1m"
            sinterval = "1min"
        elseif interval == "5m"
            sinterval = "5min"
        elseif interval == "15m"
            sinterval = "15min"
        elseif interval == "30m"
            sinterval = "30min"
        elseif interval == "1h"
            sinterval = "60min"
        end
        raw = AlphaVantage.fx_intraday(String(symbol1), String(symbol2), sinterval, outputsize="full", datatype="csv")
        df = forex_intra_to_dataframes(raw)
    end
    return df
end

#=

(CompanyOverview, BalanceSheet, IncomeStatement, CashFlow, Earnings) = fundamentals(symbol, source)

=#

function reports(symbol::Symbol)
    ssymbol = String(symbol)
    rawCO = AlphaVantage.company_overview(ssymbol::String)  # returns a dictionary
    rawBS = AlphaVantage.balance_sheet(ssymbol::String)     # returns a dictionary
    rawIS = AlphaVantage.income_statement(ssymbol::String)  # returns a dictionary
    rawCF = AlphaVantage.cash_flow(ssymbol::String)         # returns a dictionary
    rawE  = AlphaVantage.earnings(ssymbol::String)          # returns a dictionary
    return rawCO, rawBS, rawIS, rawCF, rawE
end
# Reference: http://dm13450.github.io/2021/01/01/Fundamental-AlphaVantage.html#companyoverview


#=

    listingstatus(date, state)

Input arguments:
    date - If no date is set, the API endpoint will return a list of active or delisted symbols as of the latest trading day. If a date is set, the API endpoint will "travel back" in time and return a list of active or delisted symbols on that particular date in history. Any YYYY-MM-DD date later than 2010-01-01 is supported. For example, date=2013-08-03
    state - By default, state=active and the API will return a list of actively traded stocks and ETFs. Set state=delisted to query a list of delisted assets.
=#
listingstatus = AlphaVantage.listing_status

#=

    earningscalendar(horizon)

=#
earningscalendar = AlphaVantage.earnings_calendar



# helper function
cleanup_colname!(ta::TimeArray) = TimeSeries.rename!(s -> replace(s, r"[. -]" => ""), ta, String)
"""
Description
    The fred() method is a wrapper to download financial and economic time series data from the St. Louis Federal Reserve (FRED).
Usage
    DGS = fred("DGS10")
    CPI = fred()
Method Signature(s)
    fred(data::String="CPIAUCNS")
Details
    The fred() method takes a string argument that corresponds to a series code from the St. Louis Federal
    Reserve (FRED) database. It returns the data in the TimeSeries.TimeArray data structure.  When no argument
    is provided, the default data set is the Consumer Price Index for All Urban Consumers: All Items (CPIAUCNS).
References
    https://research.stlouisfed.org/fred2
See Also
    yahoo() which is a wrapper from downloading financial time series for stocks from Yahoo Finance.
"""
function fred(series::Symbol = :CPIAUCNS)
    data = String(series)
    url = "http://research.stlouisfed.org/fred2/series/$data/downloaddata/$data.csv"
    res = HTTP.get(url)
    @assert res.status == 200
    csv = CSV.File(res.body)
    sch = TimeSeries.Tables.schema(csv)
    TimeArray(csv, timestamp = first(sch.names)) |> cleanup_colname!
end


"""
# Download stock market data from Yahoo Finance
    yahoo(symbol, from, to, interval)
## Arguments
    - `symbol`: Market symbol, e.g. :AAPL or :GOOG
    - `from` : Date/DateTime type, e.g. Date(2018,12,26)
    - `to` : Date/DateTime type, e.g. Date(2018,12,20) or DateTime(2018,12,20,8,30,0)
    - 'interval' : Sampling interval, e.g. "1d" (default), "1wk", "1mo" or "1m", "5m", "15m", "30m", "1h"
## Examples
```jldoctest
julia> yahoo("GOOG", Date(2018,12,26), Date(2018,12,20))
julia> yahoo("GOOG", Date(2018,12,26), Date(2018,12,20), "1wk")
julia> yahoo("EURUSD=X",  Date(2018,12,26), Date(2018,12,20))
```
"""
function yahoo(symbol::String, date1 = Date(1900,1,1), date2 = Date(Dates.now()), interval::String = "1d", datatype = :dataframe)
    date1 = DateTime(date1)
    date2 = DateTime(date2)
    if date1 > date2
        date1, date2 = date2, date1
    end
    from = string(round(Int64, datetime2unix(date1)))
    to = string(round(Int64, datetime2unix(date2)))
    host = rand(["query1", "query2"])
    url = "https://$host.finance.yahoo.com/v7/finance/chart/$symbol?&interval=$interval&period1=$from&period2=$to"
    response = HTTP.get(url, cookies = true)
    body = JSON.parse(String(response.body))["chart"]["result"][1]   
    values = body["indicators"]["quote"][1]

    x = DataFrame(
        Time = Dates.Date.(unix2datetime.(body["timestamp"])),
        Open = convert(Vector{Float64},values["open"]),
        High =  convert(Vector{Float64},values["high"]),
        Low = convert(Vector{Float64},values["low"]),
        Close = convert(Vector{Float64},values["close"]),
        AdjClose = convert(Vector{Float64},body["indicators"]["adjclose"][1]["adjclose"]),
        Volume = convert(Vector{Int64},values["volume"])
        )
    deleterows!(x, isnothing.(x).Close)
    if datatype == :dataframe
        return x
    elseif datatype == :timeseries
        return TimeArray(x, timestamp = :Time)
    end
end
yahoo(s::Symbol, date1, date2, interval) = yahoo(String(s), date1, date2, interval)
yahoo(s::Symbol, date1, date2) = yahoo(String(s), date1, date2, "1d")
yahoo(s::Symbol, date1) = yahoo(String(s), date1, Date(Dates.now()), "1d")
yahoo(s::Symbol) = yahoo(String(s), Date(1950,1,1), Date(Dates.now()), "1d")


function stock(symbol, interval = "1d", date1= Date(1950,1,1), date2 = Date(Dates.now()), source = :yahoo)
    if source == :yahoo
        return yahoo(symbol, date1, date2, interval)
    elseif source == :alphavantage
        return alphavantage(symbol, interval)
    end
end

function forex(symbol1, symbol2, interval = "1d", date1 = Date(1950,1,1), date2 = now(), source = :yahoo)
    if source == :yahoo
        fxpair = string(symbol1)* string(symbol2) *"=X"
        return yahoo(fxpair, date1, date2, interval)
    elseif source == :alphavantage
        return alphavantage(symbol1, symbol2, interval)
    end
end


function prices(simbol, startdate, enddate)  
    prices = yahoo(simbol, startdate, enddate)
    return prices[:AdjClose]
end


function returns(simbol, startdate, enddate, returntype) 
    p = prices(simbol, startdate, enddate)
    nobs = size(p,1)       
    returns = zeros(nobs-1)     
    if returntype == :log
        for i = 1:nobs-1
            returns[i] = log(p[i+1]/p[i])    
        end
    elseif returntype == :simple
        for i = 1:nobs-1
            returns[i] = p[i+1]/p[i]-1    
        end
    end
    return returns
end


