#=
Price Bonds (5/8):
OK  bndprice	    - Price fixed-income security from yield to maturity
    bndspread	    - Static spread over spot curve
    bndtotalreturn  - Total return of fixed-coupon bond
OK  floatmargin	    - Margin measures for floating-rate bond
    floatdiscmargin	- Discount margin for floating-rate bond
OK  prdisc          - Price of discounted security
OK  prmat           - Price with interest at maturity
OK  prtbill         - Price of Treasury bill
=#

#=
prdisc returns the price of a security whose yield is quoted as a bank discount rate (for example, U. S. Treasury bills).

    price = prdisc(settlement, maturity, face, discount, basis) 

Input arguments:
    settlement  — Settlement date
    maturity    — Maturity date
    face        - Redemption (par, face) value
    discount    - Bank discount rate of the security, specified as a decimal fraction value
    basis       - Day-count basis of instrument
Output arguments:
    price       — Price of discounted security
=#
function prdisc(settlement, maturity, face, discountrate, basis)
    face - (discountrate .* face .* yearfrac(settlement,maturity,basis))
end

#= prmat computes the price with interest at maturity. Returns the price (p) and accrued interest (ai) of a security that pays interest at maturity. This function also applies to zero-coupon bonds or pure discount securities by letting coupon = 0.

    [price, accruedinterest] = prmat(settlement, maturity, issue, face, couponrate, yield, basis)

=#
function prmat(settlement, maturity, issuedate, face, couponrate, yield, basis)
    accrint = yearfrac(issuedate,settlement,basis) .* couponrate .* face
    price = (face + (yearfrac(issuedate,maturity,basis) .* couponrate .* face))./(1+(yearfrac(settlement,maturity,basis).*yield))-accrint
    return price, accrint
end

#= prtbill returns the price for a Treasury bill.

        price = prtbill(settlement, maturity, face, discountrate)
        p = Face.*(1-Discount.*daysact(Settle,Maturity)./360)
=#
function prtbill(settlement, maturity, face, discountrate)
    return face .* (1-discountrate .* daysact(settlement,maturity)./360)
end


#=
bndprice given bonds with yields to maturity, returns the clean prices and accrued interest due.

(Price, AccruedInt) = bndprice(Yield, CouponRate, Settle, Maturity, Period, Basis, EndMonthRule, IssueDate, FirstCouponDate, LastCouponDate, StartDate, Face, LastCouponInterest, CompoundingFrequency, DiscountBasis)
=#
function bndprice(Yield, CouponRate, Settle, Maturity, Period, Basis, EndMonthRule, IssueDate, FirstCouponDate, LastCouponDate, Face, CompoundingFrequency, LastCouponInterest, DiscountBasis)

    CFlowAmounts = cfamounts(CouponRate, Settle, Maturity, Period, Basis, EndMonthRule, IssueDate, FirstCouponDate, LastCouponDate, Face, DiscountBasis)[1]
    TFactors = cfamounts(CouponRate, Settle, Maturity, Period, Basis, EndMonthRule, IssueDate, FirstCouponDate, LastCouponDate, Face, DiscountBasis)[3]

    accruedinterest = -CFlowAmounts[1]
    if ismissing(CompoundingFrequency) == true  # Frequency
        if issia(Basis) == true
            CompoundingFrequency = 2
        elseif isisma(Basis) == true
            CompoundingFrequency = 1
        elseif isa(Basis, Bus252) == true
            CompoundingFrequency = 2
        end
    end
    if length(CFlowAmounts) == 2 # settlement in last coupon period
        if LastCouponInterest == :compound 
            coef = 1 /(Yield/CompoundingFrequency +1)
            price = CFlowAmounts[2] * coef^TFactors[2] + CFlowAmounts[1]
        elseif LastCouponInterest == :simple
            price = - (-CFlowAmounts[2]/(Yield * TFactors[2] / CompoundingFrequency +1) - CFlowAmounts[1])
        end
    else
        disc = 1 / (1 + Yield/CompoundingFrequency)
        price = sum(CFlowAmounts .* disc .^ TFactors) 
    end
    return price, accruedinterest
end
bndprice(Yield, CouponRate, Settle, Maturity, Period, Basis) = bndprice(Yield, CouponRate, Settle, Maturity, Period, Basis, true, missing, missing, missing, 100, missing, :compound, missing)


#=
floatmargin computes margin measures for a floating rate bond. Use this function to calculate the following types of margin measures for a floating rate bond:
    - Spread for Life
    - Adjusted Simple Margin
    - Adjusted Total Margin
To calculate the discount margin or zero discount margin, please refer to FLOATDISCMARGIN.

    (Margin, AdjPrice) = floatmargin(Price, Spread, Settle, Maturity, SpreadType, LatestFloatingRate, StubRate, ReferenceRate, Reset, Basis, Principal, EndMonthRule, Holidays, BusinessDayConvention)

Input arguments:
    SpreadType - The type of spread to calculate. Possible values are:
        :spreadforlife (default)
        :adjustedsimple
        :adjustedtotal
    LatestFloatingRate - rate for the next floating payment set at the last reset date. This rate must be specified for SpreadType of 'adjustedsimple' and 'adjustedtotal'
    StubRate - Stub rate between the settlement date and the first coupon rate  
    SpotRate - The reference rate for the term of the floating coupons (for example, the 3 month LIBOR from settlement date for a bond with Reset of 4). StubRate and SpotRate must be specified for SpreadType of 'adjustedsimple' and 'adjustedtotal'
    Reset - Frequency of payments per year. Default is 1.
    Basis - The basis used for time factor calculations. Default is 0 (actual/actual)
    Principal - Notional principal amounts. Default is 100.
    EndMonthRule - End-of-month rule. Default is 1 (in effect).
    Holidays - Holidays to be used in computing business days.
    BusinessDayConvention - business day convention to be used in computing payment dates. Possible values are:
        :actual (default)
        :follow
        :modifiedfollow
        :previous
        :modifiedprevious
    Output arguments:
        Margin - Spreads for floating rate bond.
        AdjPrice - Adjusted price used to calculate spreads for SpreadType of 'adjustedsimple' and 'adjustedtotal'.
=#

# Spread for life
function spreadForLife(Price, Spread, MaturityAF, Principal)
    return (100 * (Principal - Price) / MaturityAF + Spread) * (Principal / Price)
end

# Adjusted price
function adjustedPrice(Price, AccrInt, AccrFactor, SpotRate, StubRate, LatestFloatingRate, Principal)
    return Price - ((LatestFloatingRate * Principal - (Price + AccrInt) * StubRate) * AccrFactor) / (1 + AccrFactor * SpotRate)
end

# Adjusted total margin
function adjustedTotal(Price, Spread, MaturityAF, SpotRate, Principal)
    Margin = (100 * (Principal - Price) / MaturityAF + Spread + 100 * (Principal - Price) * SpotRate) * (Principal / Price)
end

function floatmargin(Price, Spread, Settle, Maturity, SpreadType, LatestFloatingRate, StubRate, SpotRate, Reset, Basis, Principal, EndMonthRule, Holidays, BusinessDayConvention)
    
    MaturityAF = yearfrac(Settle, Maturity, Basis) #  Calculate number of years from settlement to maturity
    if SpreadType == :spreadforlife
        Margin = spreadForLife(Price, Spread, MaturityAF, Principal)
        return Margin
    else
        CFlowAmounts = cfamounts(LatestFloatingRate, Settle, Maturity, Reset, Basis, EndMonthRule, missing, missing, missing, Principal, missing, missing, Holidays, BusinessDayConvention)[1]
        # Obtain accrued interest
        AccrInt = -CFlowAmounts[1]
        CFlowDates = cfamounts(LatestFloatingRate, Settle, Maturity, Reset, Basis, EndMonthRule, missing, missing, missing, Principal, missing, missing, Holidays, BusinessDayConvention)[2]

        # Obtain accrual factor for the first period
        AccrFactor = yearfrac(Settle, CFlowDates[2], Basis)
        # Calculate adjusted price
        AdjPrice = adjustedPrice(Price, AccrInt, AccrFactor, SpotRate, StubRate, LatestFloatingRate, Principal)

        # Adjusted simple margin
        if SpreadType == :adjustedsimple
            Margin = spreadForLife(AdjPrice, Spread, MaturityAF, Principal)
        # Adjusted total margin
        elseif SpreadType == :adjustedtotal
            Margin = adjustedTotal(AdjPrice, Spread, MaturityAF, SpotRate, Principal)
        end
        return Margin, AdjPrice
    end   
end

#=
floatdiscmargin computes the discount margin for a floating rate bond

    Margin = floatdiscmargin(Price, Spread, Settle, Maturity, StubRate, ReferenceRate, LatestFloatingRate, Reset,
    Basis, Principal, EndMonthRule, AdjustCashFlowsBasis, Holidays, BusinessDayConvention) 

=#
function floatdiscmargin(Price, Spread, Settle, Maturity, StubRate, ReferenceRate, LatestFloatingRate, Reset,
    Basis, Principal, EndMonthRule, AdjustCashFlowsBasis, Holidays, BusinessDayConvention) 
end
