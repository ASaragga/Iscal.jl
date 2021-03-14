#=
Analyze Yield Curves (4/9):
OK  disc2zero	- Zero curve given discount curve
OK  fwd2zero	- Zero curve given forward curve
    prbyzero	- Price bonds in portfolio by set of zero curves
    pyld2zero	- Zero curve given par yield curve
    zbtprice	- Zero curve bootstrapping from coupon bond data given price
    zbtyield	- Zero curve bootstrapping from coupon bond data given yield
OK  zero2disc	- Discount curve given zero curve
OK  zero2fwd	- Forward curve given zero curve
    zero2pyld	- Par yield curve given zero curve
=#

#= 
disc2zero computes the zero curve given discount curve

    (zerorates,curvedates) = disc2zero(discountfactors, curvedates, settlement, compounding, basis)

Input arguments:
    discountfactors -Discount factors. In aggregate, the factors in discountfactors constitute a discount curve for the investment horizon represented by curvedates.
    curvedates — Maturity dates that correspond to the discount factors in discountfactors.
    settlement — Common settlement date for discount factors in discountfactors
    compounding — Rate at which output zero rates are compounded when annualized. Default: 2. Allowed values are:
        0 — Simple interest (no compounding)
        1 — Annual compounding
        2 — Semiannual compounding (default)
        3 — Compounding three times per year
        4 — Quarterly compounding
        6 — Bimonthly compounding
        12 — Monthly compounding
        365 — Daily compounding
        -1 — Continuous compounding
    basis — Day-count basis used for annualizing output zero rates
Output arguments:
    zerorates - Zero curve for the investment horizon represented by curvedates. The zero rates are the yields to maturity on theoretical zero-coupon bonds.
    curvedates - Maturity dates that correspond to the zerorates, returned as a column vector. This vector is the same as the input vector curvedates, but the output is sorted by ascending maturity. 
=#
function disc2zero(discountfactors, curvedates, settlement, compounding, basis)
    tenors = yearfrac.(settlement, curvedates, Ref(basis))
    if compounding == -1
        zerorates = -log.(discountfactors) ./ tenors
    elseif compounding == 0
        zerorates = (1.0 ./ discountfactors .- 1.0) ./ tenors
    else
        zerorates = (discountfactors .^ (-1.0 ./ (tenors .* compounding)) .- 1.0) .* compounding 
    end
    return zerorates, curvedates
end


#=
fwd2zero returns a zero curve given an implied forward rate curve and its maturity dates. 

    (zerorates, curvedates) = fwd2zero(forwardrates, curvedates, settlementent, inputcompounding, inputbasis, outputcompounding, outputbasis)

Input arguments:
    forwardrates — Annualized implied forward rates. In aggregate, the rates in forwardrates constitute an implied forward curve for the investment horizon represented by curvedates. The first element pertains to forward rates from the settlement date to the first curve date.
    curvedates — Maturity dates
    settlement — Common settlement date for forwardrates
    inputcompounding - Compounding frequency of input forward rates, specified with allowed values:
        0 — Simple interest (no compounding)
        1 — Annual compounding
        2 — Semiannual compounding (default)
        3 — Compounding three times per year
        4 — Quarterly compounding
        6 — Bimonthly compounding
        12 — Monthly compounding
        365 — Daily compounding
        -1 — Continuous compounding
    inputbasis — Day-count basis of input forward rates. numeric values: 0,1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13
    outputcompounding — Compounding frequency of output zero rates, specified with the allowed values:
        0 — Simple interest (no compounding)
        1 — Annual compounding
        2 — Semiannual compounding (default)
        3 — Compounding three times per year
        4 — Quarterly compounding
        6 — Bimonthly compounding
        12 — Monthly compounding
        365 — Daily compounding
        -1 — Continuous compounding
    outputbasis — Day-count basis of output zero rates
Output arguments:
    zerorates — Zero curve for investment horizon represented by curvedates
    curvedates — Maturity dates that correspond to zerorates
=#
function fwd2zero(forwardrates, curvedates, settlement, inputcompounding, inputbasis, outputcompounding, outputbasis)
    InYearMats = yearfrac.(settlement, curvedates, Ref(inputbasis))
    FwdYearMats = [InYearMats[1]; diff(InYearMats)]
    # Compute forward discounts
    if inputcompounding == -1
        FwdDiscounts = exp.(-forwardrates .* FwdYearMats)
    elseif inputcompounding == 0
        FwdDiscounts = 1.0 ./ (1.0 .+ forwardrates .* FwdYearMats)      
    else
        FwdDiscounts = (1.0 .+ forwardrates/inputcompounding).^(-FwdYearMats * inputcompounding)
    end
    # Compute discounts from settlement
    DiscountFromSettle = cumprod(FwdDiscounts)

    # Convert discounts to zero rates
    OutYearMats = yearfrac.(settlement, curvedates, Ref(outputbasis))
    if outputcompounding == -1
        zerorates = -log.(DiscountFromSettle)./OutYearMats
    elseif outputcompounding == 0
        zerorates = (1.0 ./ DiscountFromSettle .- 1.0) ./ OutYearMats
    else
        zerorates = outputcompounding * (DiscountFromSettle .^ (-1.0 ./ (OutYearMats * outputcompounding)) .- 1.0)
    end
    return zerorates, curvedates
end



#=
zero2disc returns a discount curve given a zero curve and its maturity dates. 

    (discfactors,curvedates) = zero2disc(zerorates, curvedates, settlement, compounding, basis)

=#
function zero2disc(zerorates, curvedates, settlement, compounding, basis)
    InYearMats = yearfrac.(settlement, curvedates, Ref(basis))
    if (compounding == -1) 
        discfactors = exp.(-zerorates .* InYearMats);
    elseif (compounding == 0)
        discfactors = 1.0 ./ (1.0 .+ zerorates .* InYearMats);
    else
        discfactors = (1.0 .+ zerorates ./ compounding) .^ (-InYearMats .* compounding);
    end
    return discfactors, curvedates
end



#=
zero2fwd returns an implied forward rate curve given a zero curve and its maturity dates. 

        (forwardrates,curvedates) = zero2fwd(zerorates, curvedates, settlement, inputcompounding, inputbasis, outputcompounding, outputbasis)

Input arguments:
    zerorates - Annualized zero rates, specified as a NUMBONDS-by-1 vector. In aggregate, the rates constitute an implied zero curve for the investment horizon represented by CurveDates. The first element pertains to forward rates from the settlement date to the first curve date.
    curvedates — Maturity dates
    settlement — Common settlement date for forwardrates
    inputcompounding - Compounding frequency of input forward rates, specified with allowed values:
        0 — Simple interest (no compounding)
        1 — Annual compounding
        2 — Semiannual compounding (default)
        3 — Compounding three times per year
        4 — Quarterly compounding
        6 — Bimonthly compounding
        12 — Monthly compounding
        365 — Daily compounding
        -1 — Continuous compounding
    inputbasis — Day-count basis of input forward rates. numeric values: 0,1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13
    outputcompounding — Compounding frequency of output zero rates, specified with the allowed values:
        0 — Simple interest (no compounding)
        1 — Annual compounding
        2 — Semiannual compounding (default)
        3 — Compounding three times per year
        4 — Quarterly compounding
        6 — Bimonthly compounding
        12 — Monthly compounding
        365 — Daily compounding
        -1 — Continuous compounding
    outputbasis — Day-count basis of output zero rates
Output arguments:
    forwardrates - Forward curve for the investment horizon represented by curvedates, returned as a NUMBONDS-by-1 vector of decimal fractions. In aggregate, the rates in ForwardRates constitute a forward curve over the dates in curvedates. 
    curvedates — Maturity dates that correspond to forwardrates
=#

function zero2fwd(zerorates, curvedates, settlement, inputcompounding, inputbasis, outputcompounding, outputbasis)
    # Get values for T in terms of fractional years from the settlement date
    InYearMats = yearfrac.(settlement, curvedates, Ref(inputbasis));

    # Get discount factors between time zero and curvedates
    if inputcompounding == -1
        DiscountFromSettle = exp.(-zerorates .* InYearMats)
    elseif inputcompounding == 0
        DiscountFromSettle = 1.0 ./ (1.0 .+ zerorates .* InYearMats)
    else
        DiscountFromSettle = (1.0 .+ zerorates ./ inputcompounding).^(-InYearMats .* inputcompounding)
    end
    
    # Compute Forward Discounts
    FwdDiscounts = append!([DiscountFromSettle[1]], DiscountFromSettle[2:end]./DiscountFromSettle[1:end-1])

    # Convert Forward Discounts to Forward Rates
    OutYearMats = yearfrac.(settlement, curvedates, Ref(outputbasis))
    FwdYearMats = [OutYearMats[1]; diff(OutYearMats)]

    if outputcompounding == -1
        forwardrates = -log.(FwdDiscounts)./FwdYearMats;
    elseif outputcompounding == 0
        forwardrates = (1.0 ./FwdDiscounts .- 1.0)./FwdYearMats;
    else
        forwardrates = outputcompounding * (FwdDiscounts.^(-1.0 ./(FwdYearMats * outputcompounding)) .- 1.0)
    end
    return forwardrates, curvedates
end






#= 
pyld2zero	returns a zero curve given a par yield curve and its maturity dates. 

        (zerorates, curvedates) = pyld2zero(parrates, curvedates, settlement, inputcompounding, inputbasis, outputcompounding, outputbasis)

Input arguments:
    parrates -  Annualized par yields. In aggregate, the rates constitute an implied zero curve for the investment horizon represented by curvedates.
    curvedates — Maturity dates
    settlement — Common settlement date for forwardrates
    inputcompounding - Compounding frequency of input forward rates, specified with allowed values:
        0 — Simple interest (no compounding)
        1 — Annual compounding
        2 — Semiannual compounding (default)
        3 — Compounding three times per year
        4 — Quarterly compounding
        6 — Bimonthly compounding
        12 — Monthly compounding
        365 — Daily compounding
        -1 — Continuous compounding
    inputbasis — Day-count basis of input forward rates. numeric values: 0,1, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13
    outputcompounding — Compounding frequency of output zero rates, specified with the allowed values:
        0 — Simple interest (no compounding)
        1 — Annual compounding
        2 — Semiannual compounding (default)
        3 — Compounding three times per year
        4 — Quarterly compounding
        6 — Bimonthly compounding
        12 — Monthly compounding
        365 — Daily compounding
        -1 — Continuous compounding
    outputbasis — Day-count basis of output zero rates
Output arguments:
    zerorates — Zero curve for investment horizon represented by curvedates
    curvedates — Maturity dates that correspond to zerorates
=#


#= 
prbyzero computes the bond prices in a portfolio using a set of zero curves. 

    bondprices = prbyzero(bonds, settlement, zerorates, zerodates, compounding) 

Input argumentos:
    bonds — Coupon bond information to compute prices, specified as a 6-column table or a NumBonds-by-6 matrix of bond information where the table columns or or matrix columns contains:
        maturity (Required) - Maturity date of the bond, as a serial date number. 
        couponrate (Required) Decimal number indicating the annual percentage rate used to determine the coupons payable on a bond.
        face (Optional) - Face or par value of the bond. Default = 100.
        period (Optional) - Coupons per year of the bond. Allowed values are 0, 1, 2 (default), 3, 4, 6, and 12.
        basis (Optional) - Day-count basis of the bond. 
        endmonthrule (Optional) - End-of-month rule. This rule applies only when Maturity is an end-of-month - date for a month having 30 or fewer days. 0 = ignore rule, meaning that a bond's coupon payment date is always the same numerical day of the month. 1 = set rule on (default), meaning that a bond's coupon payment date is always the last actual day of the month
    settlement - Settlement date
    zerorates — Observed zero rates. Each column represents a rate curve. Each row represents an observation date
    zerodates -  Observed dates for zerorates
    compounding — Compounding frequency of input zerorates when annualized, specified using the allowed values:
         1 — Annual compounding
         2 — Semiannual compounding (default)
         3 — Compounding three times per year
         4 — Quarterly compounding
         6 — Bimonthly compounding
        12 — Monthly compounding
Output arguments:
    bondprices — Clean bond prices
=# 