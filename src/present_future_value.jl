#=
Present and Future Value (5/5):
OK  pvfix	    - Present value with fixed periodic payments
OK  pvvar	    - Present value of varying cash flow
OK  fvdisc	    - Future value of discounted security
OK  fvfix	    - Future value with fixed periodic payments
OK  fvvar	    - Future value of varying cash flow
    =======================================================
OK  npv         - Net present value with periodic payments
OK  xnpv        - Net present value with variable payment periods
OK  nfv         - Net future value with periodic payments
OK  xnfv        - Net future value with variable payment periods
=#

function pvfix(rate, nperiods, payment)   # Returns the present value of a series of equal payments.
    periods = 1:1:nperiods                # Periods starting from 1
    return sum(payment./(1+ .+ rate).^periods) 
end

function fvfix(rate,nperiods,payment, starting_balance = 0)   # Returns the future value of a series of equal payments. Starting balance is the amount available in period 0.
    periods = 0:1:(nperiods-1)              # Periods starting from 0
    return sum(payment.*(1+ .+ rate).^periods) + starting_balance*(1+rate)^nperiods
end

                       

function npv(rate,cashflows)                # Present value of sequence of cashflows given an interest rate
    periods = 0:1:(length(cashflows)-1)     # Periods starting from 0
    return sum(cashflows./(1 .+ rate).^periods)
end
pvvar(cashflows, rate) = npv(rate,cashflows)

function nfv(rate,cashflows)                # Returns the future value of a varying cash flow.
    periods = length(cashflows)-1:-1:0      # Periods starting from length(cashflows)
    return sum(cashflows.*(1 .+ rate).^periods)
end
fvvar(cashflows, rate) = nfv(rate,cashflows)
    
function xnpv(rate,cashflows,dates)
    interval = map(d -> DayCounts.yearfrac(dates[1],d,DayCounts.Actual365Fixed()),dates)
    return sum(cashflows./(1 .+ rate).^interval)
end
pvvar(cashflows, rate, dates) = xnpv(rate, cashflows, dates)

function xnfv(rate,cashflows,dates)
    #interval = map(d -> DayCounts.yearfrac(d,dates[length(dates)],DayCounts.Actual365Fixed()),dates)
    interval = map(d -> DayCounts.yearfrac(d,last(dates),DayCounts.Actual365Fixed()),dates)
    return sum(cashflows.*(1 .+ rate).^interval)
end

#function fvvar(cashflows,rate,dates)
#    #interval = map(d -> DayCounts.yearfrac(d,dates[length(dates)],DayCounts.Actual365Fixed()),dates)
#    interval = map(d -> DayCounts.yearfrac(d,last(dates),DayCounts.ActualActualMatlab()),dates)
#    return sum(cashflows.*(1 .+ rate).^interval)
#end
function fvvar(rate, cfvalues, cfdates, dc)
    ncf = size(cfdates,1)
    tfrac = zeros(ncf)
    for i = 1:ncf
        tfrac[i] = yearfrac(cfdates[1], cfdates[i], dc);
    end
    return sum(cfvalues .* (1+rate).^ (maximum(tfrac).-tfrac))
end

function fvdisc(settlement, maturity, price, discount, basis)
    return price/(1- yearfrac(settlement,maturity,basis)*discount)
end

