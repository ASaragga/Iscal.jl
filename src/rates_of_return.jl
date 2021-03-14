#=
Rates of Return (6/6):
OK  effrr       - Effective rate of return
OK  irr	        - Internal rate of return
OK  mirr	    - Modified internal rate of return
OK  nomrr	    - Nominal rate of return
OK  taxedrr	    - After-tax rate of return
OK  xirr	    - Internal rate of return for nonperiodic cash flow
=#




function irr(cashflows, rm::Robust)
    f(x) = npv(x,cashflows)
    return find_zero(f, rm.x0, Order0())
end
irr(cashflows) = irr(cashflows, Robust(0))

function xirr(cashflows, dates, rm::Robust)         
    f(x) = xnpv(x, cashflows, dates)
    return find_zero(f,rm.x0, Order0())
end



function irr(cashflows, rm::Secant)
    f(x) = npv(x,cashflows)
    return find_zero(f, rm.x0, Order1())
end

function xirr(cashflows, dates, rm::Secant)         # Fastest method: Secant method
    f(x) = xnpv(x, cashflows, dates)
    return find_zero(f,rm.x0, Order1())
end
xirr(cashflows, dates) = xirr(cashflows, dates, Secant(0.0))



function irr(cashflows, rm::Steffensen)
    f(x) = npv(x,cashflows)
    return find_zero(f, rm.x0, Order2())
end

function xirr(cashflows, dates, rm::Steffensen)         # Second fastest method: Steffensen method
    f(x) = xnpv(x,cashflows,dates)
    return find_zero(f, rm.x0, Order2())
end


function nomrr(effrr,nperiods)          # effrate: Effective annual percentage rate. nperiods: number of compounding periods per year, an integer.
    return nperiods * ((1+effrr)^(1/nperiods) -1)
end


function effrr(nomrr,nperiods)          # nomrr: Annual percentage rate. nperiods: number of compounding periods per year.
    return (nomrr/nperiods+1)^nperiods-1
end

function mirr(cashflows,financing_rate,reinvestment_rate)
    len = length(cashflows)
    pcashflows = zeros(len)
    ncashflows = zeros(len)
    for i in 1:len
        if cashflows[i]>0
            pcashflows[i] = cashflows[i]
        else
            ncashflows[i] = cashflows[i]
        end
    end
    pv_ncf = npv(financing_rate,ncashflows)
    fv_pcf = nfv(reinvestment_rate,pcashflows) 
    n = len - 1
    return (-fv_pcf/pv_ncf)^(1/n) - 1
end

function taxedrr(pretax_return,taxrate)
    return pretax_return * (1-taxrate)
end

