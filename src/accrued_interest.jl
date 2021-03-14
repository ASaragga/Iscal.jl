#=
Accrued Interest (2/2):
OK  acrubond	- Accrued interest of security with periodic interest payments
OK  acrudisc	- Accrued interest of discount security paying at maturity
=#

#=
acrubond returns the accrued interest for a security with periodic interest payments with standard, short, and long first coupon periods. Note: cfamounts or accrfrac is recommended when calculating accrued interest beyond the first period.

    acrubond(issue, settlement, firstcoupon, face, couponrate, period, basis)

=#
function acrubond(issue, settlement, firstcoupon, face, couponrate, period, basis)
    cflowamounts = cfamounts(couponrate, settlement, max(firstcoupon,settlement) + Month(12), period, basis,true, issue, firstcoupon, missing, face)
    return abs(cflowamounts[1][1])
end



acrudisc(startdate, enddate, face, discount_rate) = acrudisc(startdate, enddate, face, discount_rate, 2, ActualActualMatlab())  # Default in Matlab

function acrudisc(startdate, enddate, face, discount_rate, period, basis::ActualActualMatlab) 
    return daysdif(startdate, enddate, basis) ./ cpnpersz(startdate, enddate, period, basis) .* discount_rate .* face ./ period 
end

function acrudisc(startdate, enddate, face, discount_rate, period, basis::Thirty360SIA) 
    return daysdif(startdate, enddate, basis) ./ cpnpersz(startdate, enddate, period, basis) .* discount_rate .* face ./ period 
end

function acrudisc(startdate, enddate, face, discount_rate, period, basis::Actual360) 
    return daysdif(startdate, enddate, basis) ./ cpnpersz(startdate, enddate, period, basis) .* discount_rate .* face ./ period 
end

function acrudisc(startdate, enddate, face, discount_rate, period, basis::Actual365Fixed) 
    return daysdif(startdate, enddate, basis) ./ cpnpersz(startdate, enddate, period, basis) .* discount_rate .* face ./ period 
end

function acrudisc(startdate, enddate, face, discount_rate, period, basis::Thirty360PSA) 
    return daysdif(startdate, enddate, basis) ./ cpnpersz(startdate, enddate, period, basis) .* discount_rate .* face ./ period 
end

function acrudisc(startdate, enddate, face, discount_rate, period, basis::Thirty360) 
    return daysdif(startdate, enddate, basis) ./ cpnpersz(startdate, enddate, period, basis) .* discount_rate .* face ./ period 
end

function acrudisc(startdate, enddate, face, discount_rate, period, basis::ThirtyE360) 
    return daysdif(startdate, enddate, basis) ./ cpnpersz(startdate, enddate, period, basis) .* discount_rate .* face ./ period 
end

function acrudisc(startdate, enddate, face, discount_rate, period, basis::NL365) 
    return daysdif(startdate, enddate, basis) ./ cpnpersz(startdate, enddate, period, basis) .* discount_rate .* face ./ period 
end

function acrudisc(startdate, enddate, face, discount_rate, period, basis::Bus252) 
    return daysdif(startdate, enddate, basis) ./ cpnpersz(startdate, enddate, period, basis) .* discount_rate .* face ./ period 
end

# ISMA Rule
function acrudisc(settlement, maturity, face, discount_rate, period, basis::ActualActualICMA) 
    return yearfrac(settlement, maturity, basis) * discount_rate * face /period 
end

function acrudisc(settlement, maturity, face, discount_rate, period, basis::Actual360ICMA) 
    return yearfrac(settlement, maturity, basis) * discount_rate * face /period 
end
         
function acrudisc(settlement, maturity, face, discount_rate, period, basis::Actual365ICMA) 
    return yearfrac(settlement, maturity, basis) * discount_rate * face /period 
end

function acrudisc(settlement, maturity, face, discount_rate, period, basis::Thirty360ICMA) 
    return yearfrac(settlement, maturity, basis) * discount_rate * face /period 
end

function acrudisc(settlement, maturity, face, discount_rate, period, basis::ActualActualISDA) 
    return yearfrac(settlement, maturity, basis) * discount_rate * face /period 
end
