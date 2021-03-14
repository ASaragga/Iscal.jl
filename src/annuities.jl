#=
Anuities (6/6):
OK  annurate    - Periodic interest rate of annuity
OK  annuterm    - Number of periods needed to obtain a future value
OK  payadv      - Periodic payment given number of advance payments
OK  payodd      - Payment of loan or annuity with odd first period
OK  payper      - Periodic payment of loan or annuity
OK  payuni      - Uniform payment equal to varying cash flow
=#

#=
annurate returns the periodic interest rate paid on a loan or annuity.

    rate = annurate(numperiods, payment, presentvalue, futurevalue, startofperiod, rm)

    rate = annurate(numperiods, payment, presentvalue)
    rate = annurate(numperiods, payment, presentvalue, futurevalue)
    rate = annurate(numperiods, payment, presentvalue, startofperiod)
    rate = annurate(numperiods, payment, presentvalue, futurevalue, startofperiod)

Input arguments:
    numperiods — Number of payment periods
    payment — Payment per period
    presentvalue — Present value of the loan
    futurevalue — Future value of the loan. Default: 0
    startofperiod - When payments are due, specified as a scalar integer with value of 0 (end of period) or 1 (beginning of period). Default: 0 (end of period).
    rm - Root finding method {Secant, Steffensen}. Default: Secant(0.05)

Output arguments:
    rate — Periodic interest-rate paid of a loan or annuity
=#


function annurate(numperiods, payment, presentvalue, futurevalue, startofperiod::Bool, rm::Secant)
    f(i) = payment * (1-(1+i)^(-numperiods))/i * (1+i) ^ startofperiod - presentvalue - futurevalue * (1+i)^(-numperiods)
    return find_zero(f, rm.x0, Order1())
end

function annurate(numperiods, payment, presentvalue, futurevalue, startofperiod::Bool, rm::Steffensen)
    f(i) = payment * (1-(1+i)^(-numperiods))/i * (1+i) ^ startofperiod - presentvalue - futurevalue * (1+i)^(-numperiods)
    return find_zero(f, rm.x0, Order2())
end

annurate(numperiods, payment, presentvalue) = annurate(numperiods, payment, presentvalue, 0.0, false, Secant(0.05))
annurate(numperiods, payment, presentvalue, startofperiod::Bool) = annurate(numperiods, payment, presentvalue, 0.0, startofperiod::Bool, Secant(0.05))
annurate(numperiods, payment, presentvalue, futurevalue) = annurate(numperiods, payment, presentvalue, futurevalue, false, Secant(0.05))
annurate(numperiods, payment, presentvalue, futurevalue, startofperiod::Bool) = annurate(numperiods, payment, presentvalue, futurevalue, startofperiod::Bool, Secant(0.05))



#=
annuterm returns the number of payment periods needed to obtain a future value

    numperiods = annurate(rate, payment, presentvalue, futurevalue, startofperiod, rm)

    numperiods = annurate(rate, payment, presentvalue, futurevalue)
    numperiods = annurate(rate, payment, presentvalue, futurevalue, startofperiod)

Input arguments:
    rate — Periodic interest-rate paid of a loan or annuity
    payment — Payment per period
    presentvalue — Present value of the loan
    futurevalue — Future value of the loan. Default: 0
    startofperiod - When payments are due, specified as a scalar integer with value of 0 (end of period) or 1 (beginning of period). Default: 0 (end of period).
    rm - Root finding method {Secant, Steffensen}. Default: Secant(10.0)

Output arguments:
    numperiods — Number of payment periods needed to obtain a future value
=#


function annuterm(rate, payment, presentvalue, futurevalue, startofperiod::Bool, rm::Secant)
    f(n) = presentvalue *(1+rate)^n  + payment * ((1+rate)^n-1)/rate * (1+rate)^startofperiod - futurevalue
    return find_zero(f, rm.x0, Order1())
end

function annuterm(rate, payment, presentvalue, futurevalue, startofperiod::Bool, rm::Steffensen)
    f(n) = presentvalue *(1+rate)^n  + payment * ((1+rate)^n-1)/rate * (1+rate)^startofperiod - futurevalue
    return find_zero(f, rm.x0, Order2())
end

annuterm(rate, payment, presentvalue, futurevalue) = annuterm(rate, payment, presentvalue, futurevalue, false, Steffensen(10.0))
annuterm(rate, payment, presentvalue, futurevalue, startofperiod::Bool) = annuterm(rate, payment, presentvalue, futurevalue, startofperiod::Bool, Steffensen(10.0))


#=

payadv returns the periodic payment given a number of advance payments.

    payment = payadv(rate, numperiods, presentvalue, futurevalue, advance, startofperiod, rm)

Input arguments:
    rate — Periodic interest-rate paid of a loan or annuity
    numperiods — Number of payment periods
    presentvalue — Present value of the loan
    futurevalue — Future value of the loan. Default: 0
    advance - Number of advance payments. If the payments are made at the beginning of the period, 1 is added to Advance.
    startofperiod - When payments are due, specified as a scalar integer with value of 0 (end of period) or 1 (beginning of period). Default: 0 (end of period).
    rm - Root finding method {Secant, Steffensen}. Default: Secant(presentvalue/numperiods)

Output arguments:
        payment — Periodic payment, returned as the periodic payment given a number of advance payments.

=#

function payadv(rate, numperiods, presentvalue, futurevalue, advance::Int64, startofperiod::Bool, rm::Secant)
    if startofperiod == true  # If the payments are made at the beginning of the period, 1 is added to Advance.
        advance += 1
    end 
    f(T) = advance * T + T * (1 - (1+rate)^(-(numperiods-advance)))/rate - presentvalue - futurevalue * (1+rate)^(-numperiods) 
    return find_zero(f, rm.x0, Order1())
end


function payadv(rate, numperiods, presentvalue, futurevalue, advance::Int64, startofperiod::Bool, rm::Steffensen)
    if startofperiod == true  # If the payments are made at the beginning of the period, 1 is added to Advance.
        advance += 1
    end 
    f(T) = advance * T + T * (1 - (1+rate)^(-(numperiods-advance)))/rate - presentvalue - futurevalue * (1+rate)^(-numperiods) 
    return find_zero(f, rm.x0, Order1())
end

payadv(rate, numperiods, presentvalue, futurevalue, advance::Int64) = payadv(rate, numperiods, presentvalue, futurevalue, advance::Int64, false, Secant(presentvalue/numperiods))



#=
payodd - Payment of loan or annuity with odd first period 

    payment = payodd(rate, numperiods, presentvalue, futurevalue, days, rm)

Input arguments:
    rate — Periodic interest-rate paid of a loan or annuity
    numperiods — Number of payment periods
    presentvalue — Present value of the loan
    futurevalue — Future value of the loan. Default: 0
    days - Actual number of days until the first payment is made, specified as an integer.
    rm - Root finding method {Secant, Steffensen}. Default: Secant(10.0)

Output arguments:
    payment — Periodic payment, returned as the periodic payment given a number of advance payments.

=#
function payodd(rate, nperiods, pv, fv, ndays)
    if rate == 0.0
        return (pv + fv)/nperiods
    end
    if ndays < 30 && rate != 0.0
        num = pv * (1 + rate * ndays/30) + fv * (rate+1)^(-nperiods)
        den = (1+rate) * ((1-(1+rate)^(-nperiods))/rate)
        return num/den
    end
    if ndays >= 30 && rate != 0.0
        return payper(rate, nperiods, pv, fv, false) + payper(rate, nperiods, rate * rem(ndays/30,1) * pv, 0, false)
    end
end

#=
payper returns the periodic payment of a loan or annuity.
    
    payment = payper(rate, numperiods, presentvalue, futurevalue, startofperiod) 
    
=#
function payper(rate, numperiods, presentvalue, futurevalue, startofperiod::Bool) 
    return (presentvalue + futurevalue * (1+rate)^(-numperiods))/(((1-(1+rate)^(-numperiods))/rate)*(1+rate)^startofperiod)
end
payper(rate, numperiods, presentvalue, futurevalue) = payper(rate, numperiods, presentvalue, futurevalue, false)
payper(rate, numperiods, presentvalue) = payper(rate, numperiods, presentvalue, 0, false)


#= payuni retuns an uniform payment equal to varying cash flow

    uniform = payuni(cashflows,rate)

Input arguments:
    cashflows - Cash flows, specified as a vector of varying cash flows. Include the initial investment as the initial cash flow value (a negative number).
    rate - interest rate

Output argument:
    uniform - Uniform series, returned as the value of a varying cash flow.
=#

function payuni(cashflows,rate)
    if length(cashflows) == 1
        return cashflows[1]
    else
        return payper(rate, length(cashflows)-1, npv(rate,cashflows))  
    end
end