#=
Amortization and Depreciation (6/6):
OK  amortize    - Amortization schedule
OK  depfixdb    - Fixed declining-balance depreciation schedule
OK  depgendb    - General declining-balance depreciation schedule
OK  deprdv      - Remaining depreciable value
OK  depsoyd     - Sum of years' digits depreciation
OK  depstln     - Straight-line depreciation schedule
=#

#=
amortize returns the principal and interest payments of a loan, the remaining balance of the original loan amount, and the periodic payment.

    # [principal, interest, balance, payment] = amortize(rate, numperiods, presentvalue, futurevalue, startofperiod)

Input arguments:
    rate — Interest-rate per period
    numperiods — Number of payment periods
    presentvalue — Present value of the loan
    futurevalue — Future value of the loan. Default: 0
    startofperiod — When payments are due. When payments are due, specified as a scalar integer with value of 0 (end of period) or 1 (beginning of period). Default: 0 (end of period).

Output arguments:
    principal — Principal paid in each period
    interest — Interest paid in each period
    balance — Remaining balance of the loan in each payment period
    payment — Payment per period
=#
function amortize(rate, numperiods, presentvalue, futurevalue, startofperiod::Bool)
    bstart = zeros(numperiods)
    interest = zeros(numperiods)
    principal = zeros(numperiods)
    balance = zeros(numperiods)
    payment = (presentvalue+futurevalue*(1+rate)^(-numperiods))/(((1-(1+rate)^(-numperiods))/rate)*(1+rate)^startofperiod)
    for i = 1:numperiods
        if i == 1
            bstart[i] = presentvalue
        else
            bstart[i] = balance[i-1]
        end
        if startofperiod == true && i == 1
            interest[i] = 0.0
        else
            interest[i] = rate * bstart[i]
        end
        principal[i] = payment - interest[i]
        if i == 1
            balance[i] = presentvalue - principal[i]
        else
            balance[i] = balance[i-1] - principal[i]
        end
    end
    return principal, interest, balance, payment 
end
amortize(rate, numperiods, presentvalue) = amortize(rate, numperiods, presentvalue, 0, false)

#=
depfixdb calculates the fixed declining-balance depreciation schedule.

    depreciation = depfixdb(cost, salvage, life, period, month)

Input arguments:
    cost — Initial value of the asset
    salvage — Salvage value of the asset
    life — Life of the asset in years
    period — Number of years to calculate
    month — Number of months in the first year of asset life. Default: 12

Output arguments:
    depreciation - returned as the fixed declining-balance depreciation for each Period.
=#

function depfixdb(cost, salvage, life, period::Int64, month)
    deprate = 1 - (salvage/cost)^(1/life)
    depschedule = zeros(period)
    netvalue = cost
    for i in 1:period
        if i == 1
            depschedule[i] = netvalue * deprate * month/12
        else
            depschedule[i] = netvalue * deprate
        end
        if i == period
            if month != 12
                depschedule[i] = netvalue * deprate * (1-month/12)
            end    
        end
        netvalue = netvalue - depschedule[i]
    end
    return depschedule
end
depfixdb(cost, salvage, life, period::Int64) = depfixdb(cost, salvage, life, period::Int64, 12)

#=
depgendb computes the declining-balance depreciation for each period.

    depreciation = depgendb(cost, salvage, life, factor) 

Input arguments:
    cost — Initial value of the asset
    salvage — Salvage value of the asset
    life — Life of the asset in years
    factor — Depreciation factor. When factor = 2, then the double-declining-balance method is used.

Output arguments:
    depreciation - returned as the declining-balance depreciation for each period.
=#
function depgendb(cost,salvage,life,factor) 
    csdif = cost-salvage 
    depr = zeros(life)
    for year in 1:life-1
        depr[year] = (cost * factor /life)* (1 - factor/life)^(year-1)
    end
    sumdepr = sum(depr)
    if sumdepr > csdif 
        depr[life] = 0 
    else 
        depr[life] = csdif - sumdepr 
    end 
    return depr
end

#= 
deprdv computes the remaining depreciable value for an asset.
    
    depreciation = deprdv(cost,salvage,accum)

Input arguments:
    cost — Initial value of the asset
    salvage — Salvage value of the asset
    accum - Accumulated depreciation of the asset for prior periods

Output arguments:
    depreciation - the remaining depreciable value for an asset.
=#

function deprdv(cost, salvage, accum)
    return cost - accum - salvage
end


#=
depstln computes the straight-line depreciation for an asset.

    depstln(cost,salvage,life)

Input arguments:
    cost — Initial value of the asset
    salvage — Salvage value of the asset
    life — Life of the asset in years

Output arguments: 
    depreciation - returned as a straight-line depreciation for an asset.
=#

### Cannot reproduce Matlab results when Life is fraccional.

function depstln(cost,salvage,life)
    (cost - salvage)/life
end

#=
depsoyd computes the depreciation for an asset using the sum of years' digits method.

    sum = depsoyd(cost,salvage,life)

Input arguments:
    cost — Initial value of the asset
    salvage — Salvage value of the asset
    life — Life of the asset in years

Output arguments: 
    Depreciation values, returned as a 1-by-Life vector of depreciation values with each element corresponding to a year of the asset's life.
=#

function depsoyd(cost,salvage,life)
    # Computing the years' digits
    clife = convert(Int64,ceil(life))
    dep = zeros(clife)
    for i in 1:clife
        dep[i] = life - i +1
    end
    # The sum of the digits can also be determined by using the formula (n^2+n)/2 where n is equal to - clife - the useful life of the asset in years.
    sumdigits =  (life^2 + life)/2
    dep = dep./sumdigits
    # For clife = 5, depreciation rates are as follows: 5/15 for the 1st year, 4/15 for the 2nd year, 3/15 for the 3rd year, 2/15 for the 4th year, and 1/15 for the 5th year.
    return (cost-salvage) .* dep
end

