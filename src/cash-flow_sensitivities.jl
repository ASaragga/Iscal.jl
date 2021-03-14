#=
Cash Flow Sensitivities (2/4):
OK  cfconv	    - Cash flow convexity
OK  cfdur	    - Cash flow duration and modified duration
    cfprice	    - Compute price for cash flow given yield to maturity
    cfyield	    - Compute yield to maturity for cash flow given price
=#

#=
cfconv calculates the convexity of a cash flow given the cash flow and the periodic yield. 

    convx = cfconv(cashflow,yield) 

=#
function cfconv(cashflow,yield) 
    m = length(cashflow)
    num = 0.0
    den = 0.0
    for i in 1:m
        coef = cashflow[i]/(1.0+yield)^i
        num += (i^2+i) * coef
        den += coef
    end
    return num/(den*(1+yield)^2)
end 



#=
cfdur calculates the duration and modified duration of a cash flow in periods.

    duration, modduration = cfdur(cashflows,yield)


=#
function cfdur(cashflow,yield) 
    m = length(cashflow)
    num = 0.0
    den = 0.0
    for i in 1:m
        coef = cashflow[i]/(1.0+yield)^i
        num += i * coef
        den += coef
    end
    duration = num/den
    mduration = duration/(1+yield)
    return duration, mduration
end 

#= cfprice computes price for cash flow given yield to maturity

    price = cfprice(cflowamounts, cflowdates, yield, settlement, compoundingfreq, basis)

=#