#=
Fixed-Income Securities (7/7):
OK  beytbill	- Bond equivalent yield for Treasury bill
OK  bndyield	- Yield to maturity for fixed-income security
OK  discrate	- Bank discount rate of security
OK  tbl2bond	- Treasury bond parameters given Treasury bill parameters
OK  ylddisc	    - Yield of discounted security
OK  yldmat	    - Yield with interest at maturity
OK  yldtbill	- Yield of Treasury bill
=#


function beytbill(settlement,maturity,discount)
    ndays = Dates.value(maturity-settlement)
    return (365 * discount)/(360-(discount*ndays))
end


function discrate(settlement,maturity,face,price,basis)
    return (1-price/face)*(1/yearfrac(settlement,maturity,basis))
end
discrate(settlement,maturity,face,price) = discrate(settlement,maturity,face,price,DayCounts.Actual360())


function ylddisc(settlement, maturity, face, price, basis)
    return (face/price-1)/yearfrac(settlement,maturity,basis)
end

function yldtbill(settlement, maturity, face, price)
    return (face/price-1)/yearfrac(settlement,maturity,DayCounts.Actual360())
end

function yldmat(settlement, maturity, issue, face, price, couponrate, basis)
    couponvalue = yearfrac(issue, maturity, basis) * couponrate * face
    accruedinterest = yearfrac(issue,settlement, basis) * couponrate * face
    return ((face+couponvalue)/(price+accruedinterest)-1)/yearfrac(settlement,maturity,basis) 
end

#=
tbl2bond restates US Treasury bill market parameters in US Treasury bond form as zero-coupon bonds. This function makes Treasury bills directly comparable to Treasury bonds and notes.

    (tbondmatrix, settlement) = tbl2bond(tbillmatrix)

Input arguments:
    tbillmatrix - Treasury bill parameters, specified as a 5-column table or a N-by-5 matrix of bond information where the table columns or matrix columns contains:
        Maturity - Maturity date of Treasury bills, specified as a serial date number when using a matrix. 
        DaysMaturity - Days to maturity, specified as an integer. Days to maturity are quoted on a skip-day basis; the actual number of days from settlement to maturity is DaysMaturity + 1.
        Bid - Bid bank-discount rate (the percentage discount from face value at which the bill could be bought, annualized on a simple-interest basis), specified as a decimal fraction.
        Asked - Asked bank-discount rate, specified as a decimal fraction.
        AskYield - Asked yield (the bond-equivalent yield from holding the bill to maturity, annualized on a simple-interest basis and assuming a 365-day year), specified as a decimal fraction.
Output arguments:
    tbondmatrix — Treasury bond parameters. The parameters or columns returned for tbondmatrix are:
        CouponRate - (Column 1) Coupon rate, which is always 0 since the Treasury bills are, by definition, a zero coupon instrument.
        Maturity - (Column 2) Maturity date for each bond in the portfolio as a serial date number. The format of the dates matches the format used for Maturity in TBillMatrix (serial date number, date character vector, or datetime array).
        Bid - (Column 3) Bid price based on $100 face value.
        Asked - (Column 4) Asked price based on $100 face value.
        AskYield - (Column 5) Asked yield to maturity: the effective return from holding the bond to maturity, annualized on a compound-interest basis.
    settlement - Settlement dates implied by the maturity dates and the number of days to maturity quote
    =#
    function tbl2bond(TBill)
        nbonds = size(TBill,1)
        tbond = [zeros(nbonds), TBill[:,1], 100 .- TBill[:,3] * 100 .* TBill[:,2]./360, 100 .- TBill[:,4] .* 100 .* TBill[:,2]./360, TBill[:,5]]
        settlement = TBill[:,1] - Day.(TBill[:,2].+1)
        return tbond, settlement
    end
    
    #=
    tr2bonds returns term-structure parameters (Bonds, Prices, and Yields) sorted by ascending maturity date, given Treasury bond parameters. The formats of the output matrix and vectors meet requirements for input to the zbtprice and zbtyield zero-curve bootstrapping functions

        (bonds, prices, yields) = tr2bonds(treasurymatrix, settlement)

    Input arguments:
        treasurymatrix - Treasury bond parameters, specified as a 5-column table or a NumBonds-by-5 matrix of bond information where the table columns or matrix columns contains:
            CouponRate - Coupon rate of the Treasury bond, specified as a decimal indicating the coupon rates for each bond in the portfolio.
            Maturity - Maturity date of the Treasury bond.
            Bid - Bid prices, specified using an integer-decimal form for each bond in the portfolio.
            Asked - Asked prices, specified using an integer-decimal form for each bond in the portfolio.
            AskYield - Quoted ask yield, specified using a decimal form for each bond in the portfolio.
        settlement - (Optional) Settlement date of the Treasury bond

    Output arguments:
        bonds — Coupon bond information. The parameters or columns returned for Bonds are:
            maturity - (Column 1) Maturity date for each bond in the portfolio as a serial date number. The format of the dates matches the format used for Maturity in TreasuryMatrix (serial date number, date character vector, or datetime array).
            couponrate - (Column 2) Coupon rate for each bond in the portfolio in decimal form.
            face - (Column 3, Optional) Face or par value for each bond in the portfolio. The default is 100.
            period - (Column 4, Optional) Number of coupon payments per year for each bond in the portfolio with allowed values: 1, 2, 3, 4, 6, and 12. The default is 2, unless you are dealing with zero coupons, then Period is 0 instead of 2.
            basis - (Column 5, Optional) Day-count basis for each bond in the portfolio
            endmonthrule - (Column 6, Optional) End-of-month rule flag for each bond in the portfolio. 
    =#

    #=
bndyield returns the bond equivalent yields to maturity.

    yield = bndyield(price, couponrate, settlement, maturity, period, basis, endmonthrule, issue, firstcoupon, lastcoupon, face, lastcouponinterest, compoundingfrequency, discountbasis)

Input arguments:
    price - Clean price of the bond (current price without accrued interest),
    compoundingfrequency — Compounding frequency for yield calculation. Values are:
        1 — Annual compounding
        2 — Semiannual compounding
        3 — Compounding three times per year
        4 — Quarterly compounding
        6 — Bimonthly compounding
        12 — Monthly compounding
        Note: By default, SIA bases (0-7) and BUS/252 use a semiannual compounding convention and ICMA bases (8-12) use an annual compounding convention.
    discountbasis — Basis used to compute the discount factors for computing the yield. Values are:
        0 = actual/actual
        1 = 30/360 (SIA)
        2 = actual/360
        3 = actual/365
        4 = 30/360 (PSA)
        5 = 30/360 (ISDA)
        6 = 30/360 (European)
        7 = actual/365 (Japanese)
        8 = actual/actual (ICMA)
        9 = actual/360 (ICMA)
        10 = actual/365 (ICMA)
        11 = 30/360E (ICMA)
        12 = actual/365 (ISDA)
        13 = BUS/252
        Note: If a SIA day-count basis is defined in the Basis input argument and there is no value assigned for DiscountBasis, the default behavior is for SIA bases to use the actual/actual day count to compute discount factors. If an ICMA day-count basis or BUS/252 is defined in the Basis input argument and there is no value assigned for DiscountBasis, the specified bases from the Basis input argument are used.
    lastcouponinterest — Compounding convention for computing yield of a bond in last coupon period. This is based on only the last coupon and the face value to be repaid. Acceptable values are:
        :simple
        :compound (default)
Output arguments:
    yield — Yield to maturity with semiannual compounding
=#
function bndyield(price::Real, couponrate::Real, settlement::Date, maturity::Date, period::Int64, basis, endmonthrule::Bool = true, issue::Union{Date, Missing} = missing, firstcoupon::Union{Date, Missing} = missing, lastcoupon::Union{Date, Missing} = missing, face::Union{Real, Missing} = 100.0, compoundingfrequency = missing, lastcouponinterest::Union{Symbol, Missing} = :compound, discountbasis = missing)
    if ismissing(lastcouponinterest) == true
        lastcouponinterest = :compound
    end
    if ismissing(face) == true
        face = 100
    end

    CFlowAmounts = cfamounts(couponrate, settlement, maturity, period, basis, endmonthrule, issue, firstcoupon, lastcoupon, face, discountbasis)[1]
    TFactors = cfamounts(couponrate, settlement, maturity, period, basis, endmonthrule, issue, firstcoupon, lastcoupon, face, discountbasis)[3]
    if ismissing(compoundingfrequency) == true  # Frequency
        if issia(basis) == true
            compoundingfrequency = 2
        elseif isisma(basis) == true
            compoundingfrequency = 1
        elseif isa(basis, Bus252) == true
            compoundingfrequency = 2
        end
    end
    CFlowAmounts[1] = CFlowAmounts[1] .- price
    x0 = 1 / (1 + couponrate/compoundingfrequency)
    if length(CFlowAmounts) == 2 # settlement in last coupon period
        if lastcouponinterest == :compound 
            disc = (-CFlowAmounts[1]/CFlowAmounts[2])^(1/TFactors[2])
            yield = (1/disc -1) * compoundingfrequency
        elseif lastcouponinterest == :simple
            yield = (-CFlowAmounts[2]/CFlowAmounts[1]-1)/TFactors[2] * compoundingfrequency
        end
    else # settlement before last coupon period
        f(x) = sum(CFlowAmounts .* x .^TFactors)
        disc = find_zero(f, x0, Order1(), xatol = 1e-14, atol = 1e-14) # xatol = absolute tolerance for x values, atol = absolute tolerance for f(x) values
        yield = (1/disc -1) * compoundingfrequency
    end
    return yield
end
