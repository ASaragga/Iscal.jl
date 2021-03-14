#=
Coupon Bond Dates (9/9):
OK  accrfrac	- Fraction of coupon period before settlement
OK  cpncount	- Coupon payments remaining until maturity
OK  cpndaten	- Next coupon date for fixed-income security
OK  cpndatenq	- Next quasi-coupon date for fixed-income security
OK  cpndatepq	- Previous quasi-coupon date for fixed-income security
OK  cpndatep	- Previous coupon date for fixed-income security
OK  cpndaysn	- Number of days to next coupon date
OK  cpndaysp	- Number of days since previous coupon date
OK  cpnpersz	- Number of days in coupon period


Calculating coupon dates, either actual or quasi dates, is notoriously complicated. The Financial Toolbox follows the SIA conventions in coupon date calculations.

The first step in finding the coupon dates associated with a bond is to determine the reference, or synchronization date (the sync date). Within the SIA framework, the order of precedence for determining the sync date is 
    (1) the first coupon date, 
    (2) the last coupon date, and finally 
    (3) the maturity date.
In other words, a SIA-compliant function in the Financial Toolbox first examines the FirstCouponDate input. If FirstCouponDate is specified, coupon payment dates and quasi-coupon dates are computed with respect to FirstCouponDate; if FirstCouponDate is unspecified, empty ([]), or NaN, then the LastCouponDate is examined. If LastCouponDate is specified, coupon payment dates and quasi-coupon dates are computed with respect to LastCouponDate. If both FirstCouponDate and LastCouponDate are unspecified, empty ([]), or NaN, the Maturity (a required input argument) serves as the sync date.
=#

#= cpnpersz returns the number of days in the coupon period containing the settlement date. For zero coupon bonds, coupon dates are computed as if the bonds have a semiannual coupon structure.

    cpnpersz(settlement, maturity, period, basis, endmonthrule) - Number of days in coupon period

Note: the Matlab function has more options: function cpnpersz(settlement, maturity, period, basis, endmonthrule, issuedate, firstcoupondate,lastcoupondate)
=#
function cpnpersz(settlement, maturity, period, ::ActualActualMatlab, endmonthrule::Bool = true)
    i = 1
    mtime = 12/period
    if settlement == maturity
        inflim = maturity
    else
        inflim =  eomrule(maturity, - mtime, endmonthrule)
    end
    while inflim > settlement
        i += 1
        inflim = eomrule(maturity, - mtime * i, endmonthrule)
    end
    suplim = eomrule(inflim, mtime, endmonthrule)
    return Dates.value(suplim-inflim)
end


function cpnpersz(settlement, maturity, period, ::ActualActualICMA, endmonthrule::Bool = true)
    i = 1
    mtime = 12/period
    if settlement == maturity
        inflim = maturity
    else
        inflim =  eomrule(maturity, - mtime, endmonthrule)
    end
    while inflim > settlement
        i += 1
        inflim = eomrule(maturity, - mtime * i, endmonthrule)
    end
    suplim = eomrule(inflim, mtime, endmonthrule)
    return Dates.value(suplim-inflim)
end

function cpnpersz(settlement, maturity, period, ::ActualActualISDA, endmonthrule::Bool = true)
    i = 1
    mtime = 12/period
    if settlement == maturity
        inflim = maturity
    else
        inflim =  eomrule(maturity, - mtime, endmonthrule)
    end
    while inflim > settlement
        i += 1
        inflim = eomrule(maturity, - mtime * i, endmonthrule)
    end
    suplim = eomrule(inflim, mtime, endmonthrule)
    return Dates.value(suplim-inflim)
end

function cpnpersz(settlement, maturity, period, ::Thirty360, endmonthrule::Bool = true)
    return 360 / period
end

function cpnpersz(settlement, maturity, period, ::Thirty360PSA, endmonthrule::Bool = true)
    return 360 / period
end

function cpnpersz(settlement, maturity, period, ::Thirty360ICMA, endmonthrule::Bool = true)
    return 360 / period
end

function cpnpersz(settlement, maturity, period, ::Thirty360SIA, endmonthrule::Bool = true)
    return 360 / period
end

function cpnpersz(settlement, maturity, period, ::ThirtyE360, endmonthrule::Bool = true)
    return 360 / period
end

function cpnpersz(settlement, maturity, period, ::ThirtyE360ISDA, endmonthrule::Bool = true)
    return 360 / period
end

function cpnpersz(settlement, maturity, period, ::ThirtyEPlus360, endmonthrule::Bool = true)
    return 360 / period
end

function cpnpersz(settlement, maturity, period, ::Actual360, endmonthrule::Bool = true)
    return 360 / period
end

function cpnpersz(settlement, maturity, period, ::Actual360ICMA, endmonthrule::Bool = true)
    return 360 / period
end

function cpnpersz(settlement, maturity, period, ::Actual365Fixed, endmonthrule::Bool = true)
    return 365 / period
end

function cpnpersz(settlement, maturity, period, ::NL365, endmonthrule::Bool = true)
    return 365 / period
end

function cpnpersz(settlement, maturity, period, ::Actual365ICMA, endmonthrule::Bool = true)
    return 365 / period
end

function cpnpersz(settlement::Date, maturity::Date, period, basis::Bus252, endmonthrule::Bool = true)
    i = 1
    mtime = 12/period
    if settlement == maturity
        inflim = maturity
    else
        inflim =  eomrule(maturity, - mtime, endmonthrule)
    end
    while inflim > settlement
        i += 1
        inflim = eomrule(maturity, - mtime * i, endmonthrule)
    end
    suplim = eomrule(inflim, mtime, endmonthrule)
    return bdayscount(basis.calendar, inflim, suplim)
end


# cpndatepq(settlement,maturity,period,basis) - Previous quasi-coupon date for fixed-income security
function cpndatepq(settlement::Date, maturity::Date, period, endmonthrule::Bool = true)
    i = 1
    mtime = 12/period
    if settlement == maturity
        inflim = maturity
    else
        inflim =   eomrule(maturity, - mtime, endmonthrule)
    end
    while inflim > settlement
        i += 1
        inflim = eomrule(maturity, - mtime * i, endmonthrule)
    end
    return inflim
    #suplim = inflim + Month(mtime)
    #return daysdif(settlement, suplim, basis)
end

# cpndatenq(settlement,maturity,period,basis) - Next quasi-coupon date for fixed-income security
function cpndatenq(settlement::Date, maturity::Date, period, endmonthrule::Bool = true)
    i = 1
    mtime = 12/period
    if settlement == maturity
        inflim = maturity
    else
        inflim =  eomrule(maturity, - mtime, endmonthrule)
    end
    while inflim > settlement
        i += 1
        inflim = eomrule(maturity, - mtime * i, endmonthrule)
    end
    return eomrule(inflim, mtime, endmonthrule)
    #suplim = inflim + Month(mtime)
    #return daysdif(settlement, suplim, basis)
end

# cpndaysp(settlement,maturity,period,basis) cpndaysp - Number of days since previous coupon date. The Matlab function has more options: function cpndaysp(Settle,Maturity,Period,Basis,EndMonthRule,IssueDate,FirstCouponDate,LastCouponDate)
function cpndaysp(settlement::Date, maturity::Date, period, basis, endmonthrule::Bool = true)
    i = 1
    mtime = 12/period
    if settlement == maturity
        inflim = maturity
    else
        inflim =  eomrule(maturity, - mtime, endmonthrule)
    end
    while inflim > settlement
        i += 1
        inflim = eomrule(maturity, - mtime * i, endmonthrule)
    end
    return daysdif(inflim, settlement, basis)
end


# cpndaysn(settlement,maturity,period,basis) - Number of days to next coupon date. The Matlab function has more options: function cpndaysn(Settle,Maturity,Period,Basis,EndMonthRule,IssueDate,FirstCouponDate,LastCouponDate)

function cpndaysn(settlement::Date, maturity::Date, period, basis, endmonthrule::Bool = true)
    daysdif(settlement, cpndatenq(settlement::Date, maturity::Date, period, endmonthrule), basis)
end


#= 
cpndatep returns the previous coupon date on or before settlement for a portfolio of bonds. This function finds the previous coupon date whether or not the coupon structure is synchronized with the maturity date. For zero coupon bonds the previous coupon date is the issue date, if available.

    previouscoupondate = cpndatep(settlement, maturity, period, basis, endmonthrule, issue, firstcoupon, lastcouponDate)

When bonds settled before the first coupon had been paid, cpndatep returns the issue date as the previous coupon date.
=#
function cpndatep(settlement, maturity, period, endmonthrule::Bool = true, issue::Union{Date, Missing} = missing, firstcoupon::Union{Date, Missing} = missing, lastcoupon::Union{Date, Missing} = missing)
    newdate = settlement - Month(12/period)
    cdvec = cfdates(newdate, maturity, period, endmonthrule, firstcoupon, lastcoupon) 

    if settlement < cdvec[1] && ismissing(issue) == false
        return issue
    end
    i = 1
    while cdvec[i] <= settlement
        i += 1
    end
    return cdvec[i-1]
end



#=
cpndaten returns the next coupon date after the settlement date. This function finds the next coupon date whether or not the coupon structure is synchronized with the Maturity date.

nextcoupondate = cpndaten(settlement, maturity, period, endmonthrule, firstcoupon, lastcoupon)

=#
cpndaten(settlement, maturity, period, endmonthrule::Bool = true, firstcoupon::Union{Date,Missing} = missing, lastcoupon::Union{Date,Missing} = missing) = cfdates(settlement, maturity, period, endmonthrule, firstcoupon, lastcoupon)[1]



#=
cpncount returns the whole number of coupon payments between the Settle and Maturity dates for a coupon bond or set of bonds. Coupons falling on or before Settle are not counted, except for the Maturity payment which is always counted.

cpncount(settlement, maturity, period, endmonthrule, firstcoupon, lastcoupon)

=#
function cpncount(settlement, maturity, period, endmonthrule::Bool = true, firstcoupon::Union{Date, Missing} = missing, lastcoupon::Union{Date, Missing} = missing) 
    return size(cfdates(settlement, maturity, period, endmonthrule, firstcoupon, lastcoupon),1)
end


#=
accrfrac returns the fraction of the coupon period before settlement.

accrfrac(settlement, maturity, period, basis, endmonthrule, issue, firstcouponrate,lastcoupon)

 =#
 # helper function
function den(startdate, enddate, period, dc)
    if isa(dc, ActualActualMatlab) == true || isa(dc, ActualActualICMA) == true || isa(dc, ActualActualISDA) == true
        return daysdif(startdate, enddate, dc)
    elseif isa(dc, Thirty360) == true || isa(dc, Thirty360PSA) == true || isa(dc, Thirty360ICMA) == true || isa(dc, Thirty360SIA) == true || isa(dc, ThirtyE360) == true || isa(dc, ThirtyE360ISDA) == true || isa(dc, ThirtyEPlus360) == true || isa(dc, Actual360) == true || isa(dc, Actual360ICMA) == true
        return 360 / period
    elseif isa(dc, Actual365Fixed) == true || isa(dc, NL365) == true || isa(dc, Actual365ICMA) == true
        return 365 / period
    elseif isa(dc, Bus252) == true
        return cpnpersz(startdate, enddate, period, Bus252())
    end
end

function cpnpersz(settlement, maturity, period, ::ActualActualMatlab)
    i = 1
    mtime = 12/period
    if settlement == maturity
        inflim = maturity
    else
        inflim =  maturity - Month(mtime)
    end
    while inflim > settlement
        i += 1
        inflim = maturity - Month(mtime * i)
    end
    suplim = inflim + Month(mtime)
    return Dates.value(suplim-inflim)
end

function cpnpersz(settlement, maturity, period, ::ActualActualICMA)
    i = 1
    mtime = 12/period
    if settlement == maturity
        inflim = maturity
    else
        inflim =  maturity - Month(mtime)
    end
    while inflim > settlement
        i += 1
        inflim = maturity - Month(mtime * i)
    end
    suplim = inflim + Month(mtime)
    return Dates.value(suplim-inflim)
end

function cpnpersz(settlement, maturity, period, ::ActualActualISDA)
    i = 1
    mtime = 12/period
    if settlement == maturity
        inflim = maturity
    else
        inflim =  maturity - Month(mtime)
    end
    while inflim > settlement
        i += 1
        inflim = maturity - Month(mtime * i)
    end
    suplim = inflim + Month(mtime)
    return Dates.value(suplim-inflim)
end

function cpnpersz(settlement, maturity, period, ::Thirty360)
    return 360 / period
end

function cpnpersz(settlement, maturity, period, ::Thirty360PSA)
    return 360 / period
end

function cpnpersz(settlement, maturity, period, ::Thirty360ICMA)
    return 360 / period
end

function cpnpersz(settlement, maturity, period, ::Thirty360SIA)
    return 360 / period
end

function cpnpersz(settlement, maturity, period, ::ThirtyE360)
    return 360 / period
end

function cpnpersz(settlement, maturity, period, ::ThirtyE360ISDA)
    return 360 / period
end

function cpnpersz(settlement, maturity, period, ::ThirtyEPlus360)
    return 360 / period
end

function cpnpersz(settlement, maturity, period, ::Actual360)
    return 360 / period
end

function cpnpersz(settlement, maturity, period, ::Actual360ICMA)
    return 360 / period
end

function cpnpersz(settlement, maturity, period, ::Actual365Fixed)
    return 365 / period
end

function cpnpersz(settlement, maturity, period, ::NL365)
    return 365 / period
end

function cpnpersz(settlement, maturity, period, ::Actual365ICMA)
    return 365 / period
end

function cpnpersz(settlement::Date, maturity::Date, period, basis::Bus252)
    i = 1
    mtime = 12/period
    if settlement == maturity
        inflim = maturity
    else
        inflim =  maturity - Month(mtime)
    end
    while inflim > settlement
        i += 1
        inflim = maturity - Month(mtime * i)
    end
    suplim = inflim + Month(mtime)
    return bdayscount(basis.calendar, inflim, suplim)
end


function accrfrac(settlement, maturity, period, basis, endmonthrule::Bool = true, issue::Union{Date, Missing} = missing, firstcoupon::Union{Date, Missing} = missing, lastcoupon::Union{Date, Missing} = missing)    
    startdate = cpndatep(settlement, maturity, period, endmonthrule, issue, firstcoupon, lastcoupon)
    enddate = cpndaten(settlement, maturity, period, endmonthrule, firstcoupon, lastcoupon)
    num = daysdif(startdate, settlement, basis)

    if ismissing(issue) == true    
        return num/den(startdate, enddate, period, basis)
    else
        if startdate == issue
            return num/den(enddate-Month(12/period), enddate, period, basis)
        else    
            return num/den(startdate, enddate, period, basis)
        end
    end
end

accrfrac(settlement, maturity, period, basis, endmonthrule, issue) = accrfrac(settlement, maturity, period, basis, endmonthrule, issue, missing, missing)

