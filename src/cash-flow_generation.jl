#=
Cash Flow Generation (5/6):
OK  cfamounts	- Cash flow and time mapping for bond portfolio
    cfport	    - Portfolio form of cash flow amounts
OK  cftimes	    - Time factors corresponding to bond cash flow dates
OK  cfdates	    - Cash flow dates for fixed-income security
OK  cfdatesq	- Quasi-coupon dates for fixed-income security
OK  tmfactor	- Time factors of arbitrary dates
=#


#=
    [CFlowAmounts,CFlowDates,TFactors,CFlowFlags,CFPrincipal] = cfamounts(CouponRate,Settle,Maturity, Period, Basis; EndMonthRule::Bool = true, IssueDate::Date, FirstCouponDate::Date, LastCouponDate::Date, StartDate::Date, Face = 100, AdjustCashFlowsBasis::Bool = false, BusinessDayConvention = :actual, CompoundingFrequency, DiscountBasis = 0, Holidays, PrincipalType = :sinking)

Optional Input Arguments:

'EndMonthRule', default: 1 (in effect) - This rule applies only when Maturity is an end-of-month date for a month having 30 or fewer days.
    0 = Ignore rule, meaning that a bond coupon payment date is always the same numerical day of the month.
    1 = Set rule on, meaning that a bond coupon payment date is always the last actual day of the month.

'IssueDate' — Bond issue date (the date the bond begins to accrue interest).

'FirstCouponDate' — Irregular or normal first coupon date.
    When FirstCouponDate and LastCouponDate are both specified, the FirstCouponDate takes precedence in determining the coupon payment structure. If FirstCouponDate is not specified, then LastCouponDate determines the coupon structure of the bond.

'LastCouponDate' — Irregular or normal last coupon date.

'StartDate' — Forward starting date of coupon payments after the Settle date.

'Face', default: 100 — Face value of bond.

'AdjustCashFlowsBasis', default:false - Adjusts cash flows according to accrual amount based on actual period day count

'BusinessDayConvention' default: 'actual' - Business day conventions. The selection for business day convention determines how nonbusiness days are treated. Nonbusiness days are defined as weekends plus any other date that businesses are not open (for example, statutory holidays). Values are:
    'actual' — Nonbusiness days are effectively ignored. Cash flows that fall on non-business days are assumed to be distributed on the actual date.
    'follow' — Cash flows that fall on a nonbusiness day are assumed to be distributed on the following business day.
    'modifiedfollow' — Cash flows that fall on a non-business day are assumed to be distributed on the following business day. However if the following business day is in a different month, the previous business day is adopted instead.
    'previous' — Cash flows that fall on a nonbusiness day are assumed to be distributed on the previous business day.
    'modifiedprevious' — Cash flows that fall on a nonbusiness day are assumed to be distributed on the previous business day. However if the previous business day is in a different month, the following business day is adopted instead.

'CompoundingFrequency' — Compounding frequency for yield calculation. Default: SIA bases (0-7) and BUS/252 use a semiannual compounding convention and ICMA bases (8-12) use an annual compounding convention. Values are:
    1 — Annual compounding
    2 — Semiannual compounding
    3 — Compounding three times per year
    4 — Quarterly compounding
    6 — Bimonthly compounding
    12 — Monthly compounding

'DiscountBasis', default: 0 — Basis used to compute the discount factors for computing the yield. Values are:
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

If a SIA day-count basis (0-7) is defined in the Basis input argument and there is no value assigned for DiscountBasis, the default behavior is for SIA bases to use the actual/actual day count to compute discount factors. If an ICMA day-count basis (8-12) or BUS/252 is defined in the Basis input argument and there is no value assigned for DiscountBasis, the specified bases from the Basis input argument are used.

'Holidays' — Dates for holidays. Holidays are used in computing business days.

'PrincipalType', default: sinking — Type of principal when a Face schedule is specified
    If 'sinking', principal cash flows are returned throughout the life of the bond.
    If 'bullet', principal cash flow is only returned at maturity.

Output Arguments:

CFlowAmounts — Cash flow amounts. The first entry in each row vector is the accrued interest due at settlement. This amount could be zero, positive or negative. If no accrued interest is due, the first column is zero. If the bond is trading ex-coupon then the accrued interest is negative.

CFlowDates — Cash flow dates for a portfolio of bonds.

TFactors — Matrix of time factors for a portfolio of bonds. Time factors are for price/yield conversion and time factors are in units of whole semiannual coupon periods plus any fractional period using an actual day count. The term time factors refer to the exponent TF in the discounting equation
\[PV = sum_{i=1}^n\left(\frac{CF}{(1+z/f)^\text{TF}}\right)\]
where:
    PV = Present value of a cash flow.
    CF = Cash flow amount.
    z = Risk-adjusted annualized rate or yield corresponding to a given cash flow. The yield is quoted on a semiannual basis.
    f = Frequency of quotes for the yield. Default is 2 for Basis values 0 to 7 and 13 and 1 for Basis values 8 to 12. The default can be overridden by specifying the CompoundingFrequency name-value pair.
    TF = Time factor for a given cash flow. The time factor is computed using the compounding frequency and the discount basis. If these values are not specified, then the defaults are as follows: CompoundingFrequency default is 2 for Basis values 0 to 7 and 13 and 1 for Basis values 8 to 12. DiscountBasis is 0 for Basis values 0 to 7 and 13 and the value of the input Basis for Basis values 8 to 12.

CFlowFlags — Cash flow flags for a portfolio of bonds. Flags identify the type of each cash flow (for example, nominal coupon cash flow, front, or end partial, or "stub" coupon, maturity cash flow).
    0 - Accrued interest due on a bond at settlement.
    1 - Initial cash flow amount smaller than normal due to a “stub” coupon period. A stub period is created when the time from issue date to first coupon date is shorter than normal.
    2 - Larger than normal initial cash flow amount because the first coupon period is longer than normal.
    3 - Nominal coupon cash flow amount.
    4 - Normal maturity cash flow amount (face value plus the nominal coupon amount).
    5 - End “stub” coupon amount (last coupon period is abnormally short and actual maturity cash flow is smaller than normal).
    6 - Larger than normal maturity cash flow because the last coupon period longer than normal.
    7 - Maturity cash flow on a coupon bond when the bond has less than one coupon period to maturity.
    8 - Smaller than normal maturity cash flow when the bond has less than one coupon period to maturity.
    9 - Larger than normal maturity cash flow when the bond has less than one coupon period to maturity.
    10 - Maturity cash flow on a zero coupon bond.
    11 - Sinking principal and initial cash flow amount smaller than normal due to a "stub" coupon period. A stub period is created when the time from issue date to first coupon date is shorter than normal.
    12 - Sinking principal and larger than normal initial cash flow amount because the first coupon period is longer than normal.
    13 - Sinking principal and nominal coupon cash flow amount.

CFPrincipal — Principal cash flows
    If PrincipalType is 'sinking', CFPrincipal output indicates when the principal is returned.
    If PrincipalType is 'bullet', CFPrincipal is all zeros and, at Maturity, the appropriate Face value.
=#

# helper functions: stubs, cfamounts_zero
function stubs(startdate, enddate, period, endmonthrule, basis)
    num = daysdif(startdate, enddate, basis)
    den = daysdif(eomrule(enddate, -12/period, endmonthrule), enddate, basis)
    return num/den
end

function cfamounts_zero(settlement, maturity, face, compoudingfrequency, discountbasis, ncf = 2)

    CFlowAmounts = Array{Union{Missing, Float64}}(missing,1,ncf)
    CFlowDates = Array{Union{Missing, Date}}(missing, 1, ncf)
    TFactors = Array{Union{Missing, Float64}}(missing, 1, ncf)
    CFlowPrincipal = Array{Union{Missing, Float64}}(missing, 1, ncf)
    
    CFlowAmounts[1] = 0
    CFlowAmounts[2] = face

    CFlowDates[1] = settlement 
    CFlowDates[2] = maturity
    
    k = 0
    nxtc = CFlowDates[2]
    while nxtc > settlement + Month(12/compoudingfrequency)  # TODO: Check case of endmonthrule = true
        nxtc = nxtc - Month(12/compoudingfrequency)
        k += 1
    end
    num = daysdif(settlement, nxtc, discountbasis)
    den = daysdif(nxtc - Month(12/compoundingfrequency), nxtc, discountbasis)
    TFactors[1] = 0
    TFactors[2] = k + num/den

    CFlowPrincipal[1] = 0
    CFlowPrincipal[2] = face

    return CFlowAmounts, CFlowDates,  TFactors, CFlowPrincipal
end



#=====================================
firstcoupon = missing
lastcoupon = missing
=====================================#
function cfamounts1(couponrate, settlement, maturity, period, basis, endmonthrule, issue, face, compoundingfrequency, discountbasis, bdayconvention, holidays)
    
    if ismissing(compoundingfrequency) == true  # Frequency
        if issia(basis) == true
            compoundingfrequency = 2
        elseif isisma(basis) == true
            compoundingfrequency = 1
        elseif isa(basis, Bus252) == true
            compoundingfrequency = 2
        end
    end

    if ismissing(discountbasis) == true
        if issia(basis) == true
            discountbasis = ActualActualMatlab()
        elseif isisma(basis) == true || isa(basis, Bus252) == true
            discountbasis = basis
        end
    end

    if isa(discountbasis, Actual360ICMA) == true  # adjustment in TFactors made by Matlab when basis is Actual360ICMA
        coefxx = 365/360
    else
        coefxx = 1.0
    end

    if period == 0  # zero coupon bonds
        ncf = 2
        CFlowAmounts = Array{Union{Missing, Float64}}(missing,1,ncf)
        CFlowDates = Array{Union{Missing, Date}}(missing, 1, ncf)
        TFactors = Array{Union{Missing, Float64}}(missing, 1, ncf)
        CFlowPrincipal = Array{Union{Missing, Float64}}(missing, 1, ncf)

        CFlowAmounts[1] = 0
        CFlowAmounts[2] = face
    
        CFlowDates[1] = settlement 
        CFlowDates[2] = maturity
        CFlowDates = busdate.(CFlowDates, Ref(bdayconvention), Ref(holidays)) # Adjust Business Days Convention

        k = 0
        nxtc = CFlowDates[2]
        while nxtc > settlement + Month(12/compoundingfrequency)
            nxtc = nxtc - Month(12/compoundingfrequency)
            k += 1
        end
        num = daysdif(settlement, nxtc, discountbasis)
        den = daysdif(nxtc - Month(12/compoundingfrequency), nxtc, discountbasis) 

        TFactors[1] = 0
        TFactors[2] = k + num/den 
    
        CFlowPrincipal[1] = 0
        CFlowPrincipal[2] = face
    
        return CFlowAmounts, CFlowDates,  TFactors, CFlowPrincipal
    end
  
    ncoupons = cpncount(settlement, maturity, period, endmonthrule)
    ncf = ncoupons + 1

    # CFlowDates
    CFlowDates = Array{Union{Missing, Date}}(missing, 1, ncf)
    CFlowDates[1] = settlement
    for i = 2:ncf
        CFlowDates[i] = cfdates(settlement, maturity, period, endmonthrule)[i-1]
    end
    CFlowDates = busdate.(CFlowDates, Ref(bdayconvention), Ref(holidays)) # Adjust Business Days Convention


    # CFlowPrincipal
    CFlowPrincipal = Array{Union{Missing, Float64}}(missing,1,ncf)
    for i = 1:ncf-1
        CFlowPrincipal[i] = 0.0
    end
    CFlowPrincipal[ncf] = face

    
    # CFlowAmounts
    CFlowAmounts = Array{Union{Missing, Float64}}(missing, 1, ncf)
    CFlowAmounts[1] = - accrfrac(settlement, maturity, period, basis, endmonthrule, issue) * couponrate * face / period
    for i = 2:ncf-1
        CFlowAmounts[i] = face * couponrate / period
    end
    CFlowAmounts[ncf] = face * (1 + couponrate/period)    

    
    # TFactors
    TFactors = Array{Union{Missing, Float64}}(missing, 1, ncf)
    for i = 1:ncf
        ### Make time factors in units of whole semi-annual coupon periods plus any fractional period using an Actual day count: multiple semi-annual periods (k) + fracctinal period (num/den)
        nxtc = CFlowDates[i]
        k = 0
        while nxtc > settlement + Month(12/compoundingfrequency)
            nxtc = nxtc - Month(12/compoundingfrequency)
            k += 1
        end
        num = daysdif(settlement, nxtc, discountbasis)
        den = daysdif(nxtc - Month(12/compoundingfrequency), nxtc, discountbasis) 
        TFactors[i] = k + num/den * coefxx
    end
  
    return CFlowAmounts, CFlowDates,  TFactors, CFlowPrincipal
end



#=====================================
firstcoupon = not missing
lastcoupon = missing
=====================================#
function cfamounts2(couponrate, settlement, maturity, period, basis, endmonthrule, issue, firstcoupon, lastcoupon, face, compoundingfrequency, discountbasis, bdayconvention, holidays)
    
    if ismissing(compoundingfrequency) == true  # Frequency
        if issia(basis) == true
            compoundingfrequency = 2
        elseif isisma(basis) == true
            compoundingfrequency = 1
        elseif isa(basis, Bus252) == true
            compoundingfrequency = 2
        end
    end

    if ismissing(discountbasis) == true
        if issia(basis) == true
            discountbasis = ActualActualMatlab()
        elseif isisma(basis) == true || isa(basis, Bus252) == true
            discountbasis = basis
        end
    end
    if isa(discountbasis, Actual360ICMA) == true # adjustment in TFactors made by Matlab when basis is Actual360ICMA
        coefxx = 365/360
    else
        coefxx = 1.0
    end

    if period == 0  # zero coupon bonds
        cfamounts_zero(settlement, maturity, face, compoudingfrequency, discountbasis)
    end
    
    ncoupons = cpncount(settlement, maturity, period, endmonthrule, firstcoupon, lastcoupon)
    ncf = ncoupons + 1

    # CFlowDates
    CFlowDates = Array{Union{Missing, Date}}(missing, 1, ncf)
    CFlowDates[1] = settlement
    for i = 2:ncf
        CFlowDates[i] = cfdates(settlement, maturity, period, endmonthrule, firstcoupon, lastcoupon)[i-1]
    end
    CFlowDates = busdate.(CFlowDates, Ref(bdayconvention), Ref(holidays)) # Adjust Business Days Convention

    # CFlowPrincipal
    CFlowPrincipal = Array{Union{Missing, Float64}}(missing,1,ncf)
    for i = 1:ncf-1
        CFlowPrincipal[i] = 0.0
    end
    CFlowPrincipal[ncf] = face
    

    # CFlowAmounts
    CFlowAmounts = Array{Union{Missing, Float64}}(missing, 1, ncf)

    # Compute CFlowAmounts[1], the amount of accrued interest
    CFlowAmounts[1] = - accrfrac(settlement, maturity, period, basis, endmonthrule, issue, firstcoupon, lastcoupon) * couponrate * face / period
    
    # Compute CFlowAmounts[2]. Test if period [2] is regular 
    if CFlowDates[2] == firstcoupon && ismissing(issue) == false # Potentially a non regular period
        CFlowAmounts[2] = stubs(issue, firstcoupon, period, endmonthrule, basis) * face * couponrate / period
    else
        CFlowAmounts[2] = face * couponrate / period
    end
    # Compute CFlowAmounts[3] : CFlowAmounts[ncf-1]
    for i = 3:ncf
        CFlowAmounts[i] = stubs(CFlowDates[i-1], CFlowDates[i], period, endmonthrule, basis) * face * couponrate / period
    end
    CFlowAmounts[ncf] += face 
    

    TFactors = Array{Union{Missing, Float64}}(missing, 1, ncf)
    for i = 1:ncf
        ### Make time factors in units of whole semi-annual coupon periods plus any fractional period using an Actual day count: multiple semi-annual periods (k) + fracctinal period (proportion of periods of 6 months)
        nxtc = CFlowDates[i]
        k = 0
        while nxtc > settlement + Month(12/compoundingfrequency)
            nxtc = nxtc - Month(12/compoundingfrequency)
            k += 1
        end
        TFactors[i] = k + stubs(settlement, nxtc, 2, endmonthrule, discountbasis) * coefxx
    end
  
    return CFlowAmounts, CFlowDates,  TFactors, CFlowPrincipal
end


#=====================================
firstcoupon = missing
lastcoupon = not missing
=====================================#
function cfamounts3(couponrate, settlement, maturity, period, basis, endmonthrule, issue, firstcoupon, lastcoupon, face, compoundingfrequency, discountbasis, bdayconvention, holidays)

    if ismissing(compoundingfrequency) == true  # Frequency
        if issia(basis) == true
            compoundingfrequency = 2
        elseif isisma(basis) == true
            compoundingfrequency = 1
        elseif isa(basis, Bus252) == true
            compoundingfrequency = 2
        end
    end

    if ismissing(discountbasis) == true
        if issia(basis) == true
            discountbasis = ActualActualMatlab()
        elseif isisma(basis) == true || isa(basis, Bus252) == true
            discountbasis = basis
        end
    end

    if isa(discountbasis, Actual360ICMA) == true # adjustment in TFactors made by Matlab when basis is Actual360ICMA
        coefxx = 365/360
    else
        coefxx = 1.0
    end

    if period == 0  # zero coupon bonds
        cfamounts_zero(settlement, maturity, face, compoudingfrequency, discountbasis)
    end
    
    ncoupons = cpncount(settlement, maturity, period, endmonthrule, firstcoupon, lastcoupon)
    ncf = ncoupons + 1

    # CFlowDates
    CFlowDates = Array{Union{Missing, Date}}(missing, 1, ncf)
    CFlowDates[1] = settlement
    for i = 2:ncf
        CFlowDates[i] = cfdates(settlement, maturity, period, endmonthrule, firstcoupon, lastcoupon)[i-1]
    end
    CFlowDates = busdate.(CFlowDates, Ref(bdayconvention), Ref(holidays)) # Adjust Business Days Convention


    # CFlowPrincipal
    CFlowPrincipal = Array{Union{Missing, Float64}}(missing,1,ncf)
    for i = 1:ncf-1
        CFlowPrincipal[i] = 0.0
    end
    CFlowPrincipal[ncf] = face
    

    # CFlowAmounts
    CFlowAmounts = Array{Union{Missing, Float64}}(missing, 1, ncf)

    # Compute CFlowAmounts[1], the amount of accrued interest
    CFlowAmounts[1] = - accrfrac(settlement, maturity, period, basis, endmonthrule, issue, firstcoupon, lastcoupon) * couponrate * face / period
    
    # Compute CFlowAmounts[2]. Test if period [2] is regular 
    if CFlowDates[2] - Month(12/period) < issue && ismissing(issue) == false # Potentially a non regular period
        CFlowAmounts[2] = stubs(issue, CFlowDates[2], period, endmonthrule, basis) * face * couponrate / period
    else
        CFlowAmounts[2] = face * couponrate / period
    end
    # Compute CFlowAmounts[3] : CFlowAmounts[ncf-1]
    for i = 3:ncf
        CFlowAmounts[i] = stubs(CFlowDates[i-1], CFlowDates[i], period, endmonthrule, basis) * face * couponrate / period
    end
    CFlowAmounts[ncf] += face 
    

    # TFactors
    TFactors = Array{Union{Missing, Float64}}(missing, 1, ncf)
    for i = 1:ncf
        ### Make time factors in units of whole semi-annual coupon periods plus any fractional period using an Actual day count: multiple semi-annual periods (k) + fracctinal period (num/den)
        nxtc = CFlowDates[i]
        k = 0
        while nxtc > settlement + Month(12/compoundingfrequency)
            nxtc = nxtc - Month(12/compoundingfrequency)
            k += 1
        end
        TFactors[i] = k + stubs(settlement, nxtc, 2, endmonthrule, discountbasis) * coefxx
    end
    return CFlowAmounts, CFlowDates,  TFactors, CFlowPrincipal
end

function cfamounts(couponrate, settlement::Date, maturity::Date, period::Int64, basis, endmonthrule::Bool = true, issue::Union{Date, Missing} = missing, firstcoupon::Union{Date, Missing} = missing, lastcoupon::Union{Date, Missing} = missing, face::Real = 100.0, compoundingfrequency = missing, discountbasis = missing, bdayconvention = Unadjusted(), holidays::BusinessDays.HolidayCalendar = WeekendsOnly())

    if ismissing(firstcoupon) == true && ismissing(lastcoupon) == true
        return cfamounts1(couponrate, settlement, maturity, period, basis, endmonthrule, issue, face, compoundingfrequency, discountbasis, bdayconvention, holidays)
    elseif ismissing(firstcoupon) == false && ismissing(lastcoupon) == true
        return cfamounts2(couponrate, settlement, maturity, period, basis, endmonthrule, issue, firstcoupon, lastcoupon, face, compoundingfrequency, discountbasis, bdayconvention, holidays)
    elseif ismissing(firstcoupon) == true && ismissing(lastcoupon) == false
        return cfamounts3(couponrate, settlement, maturity, period, basis, endmonthrule, issue, firstcoupon, lastcoupon, face, compoundingfrequency, discountbasis, bdayconvention, holidays)
    elseif ismissing(firstcoupon) == false && ismissing(lastcoupon) == false
        return cfamounts2(couponrate, settlement, maturity, period, basis, endmonthrule, issue, firstcoupon, lastcoupon, face, compoundingfrequency, discountbasis, bdayconvention, holidays)
    end
end
cfamounts(couponrate, settlement, maturity, period, basis, endmonthrule, issue)  = cfamounts(couponrate, settlement, maturity, period, basis, endmonthrule, issue, missing, missing, 100.0, missing, Unadjusted(), WeekendsOnly())
cfamounts(couponrate, settlement, maturity, period, basis)  = cfamounts(couponrate, settlement, maturity, period, basis, true, missing, missing, missing, 100.0, missing, Unadjusted(), WeekendsOnly())


#=
cfdates returns a matrix of cash flow dates for a bond or set of bonds. cfdates determines all cash flow dates for a bond whether or not the coupon payment structure is normal or the first and/or last coupon period is long or short.
    
    cfdates(settlement::Date, maturity::Date, period::Int, eom::Bool, lastcoupon, firstcoupon)

    Note: Signature is different from MATLAB

=#
# helper function
function incr(eom,anchor,gap)
    if eom == true
        return eomrule(anchor,gap)
    else
        return anchor + Month(gap)
    end
end

function cfdates(settlement::Date, maturity::Date, period::Int, eom::Bool = true, firstcoupon::Union{Date, Missing} = missing, lastcoupon::Union{Date, Missing} = missing)
    if ismissing(firstcoupon) == true && ismissing(lastcoupon) == true
        cfdvector = Vector{Date}()
        date = maturity
        nmper = 12/period
        m = 1
        while date > settlement
            pushfirst!(cfdvector, date)   #push!
            if eom == true
                date = eomrule(maturity, - m * nmper)
            else
                date = maturity - Month(m * nmper)
            end 
            m +=1
        end
        return cfdvector
    elseif ismissing(firstcoupon) == true && ismissing(lastcoupon) == false
        cfdvector = [maturity]
        date = lastcoupon
        nmper = 12/period
        m = 1
        while date > settlement
            pushfirst!(cfdvector, date)
            if eom == true
                date = eomrule(lastcoupon, -m * nmper)
            else
                date = lastcoupon - Month(m * nmper)
            end
            m += 1
        end
        return cfdvector
    elseif ismissing(firstcoupon) == false && ismissing(lastcoupon) == true
        if settlement <= firstcoupon
            cfdvector = Vector{Date}()
            date = firstcoupon
            nmper = 12/period
            m = 1
            while date < maturity
                push!(cfdvector, date)
                if eom == true
                    date = eomrule(firstcoupon, m * nmper)
                else
                    date = firstcoupon + Month(m * nmper)
                end
                m += 1
            end
            push!(cfdvector, maturity)
            return cfdvector
        else
            cfdvector = Vector{Date}()
            date = firstcoupon
            nmper = 12/period
            m = 1
            while date < settlement
                if eom == true
                    date = eomrule(firstcoupon, m * nmper)
                else
                    date = firstcoupon + Month(m * nmper)
                end
                m += 1
            end
            while date < maturity
                push!(cfdvector, date)
                if eom == true
                    date = eomrule(firstcoupon, m * nmper)
                else
                    date = firstcoupon + Month(m * nmper)
                end
                m += 1
            end
            push!(cfdvector, maturity)
            return cfdvector    
        end
    elseif ismissing(firstcoupon) == false && ismissing(lastcoupon) == false
        if settlement <= firstcoupon
            cfdvector = Vector{Date}()
            date = firstcoupon
            nmper = 12/period
            m = 1
            while date < lastcoupon
                push!(cfdvector, date)
                if eom == true
                    date = eomrule(firstcoupon, m * nmper)
                else
                    date = firstcoupon + Month(m * nmper)
                end
                m += 1
            end
            push!(cfdvector, lastcoupon, maturity)
            return cfdvector
        else
            cfdvector = Vector{Date}()
            date = firstcoupon
            nmper = 12/period
            m = 1
            while date < settlement
                if eom == true
                    date = eomrule(firstcoupon, m * nmper)
                else
                    date = firstcoupon + Month(m * nmper)
                end
                m += 1
            end
            while date < lastcoupon
                push!(cfdvector, date)
                if eom == true
                    date = eomrule(firstcoupon, m * nmper)
                else
                    date = firstcoupon + Month(m * nmper)
                end
                m += 1
            end
            push!(cfdvector, lastcoupon, maturity)
            return cfdvector    
        end
    end     
end



#=
cfdatesq returns a matrix of quasi-coupon dates.

    cfdatesq(settlement::Date, maturity::Date, period::Int, eom::Bool, firstcoupon, lastcoupon)

=#
function cfdatesq(settlement::Date, maturity::Date, period::Int, eom::Bool, firstcoupon, lastcoupon, periods_before_settlement::Int, periods_after_maturity::Int)
    cfdvector = Vector{Date}()
    nmper = 12/period
    if ismissing(firstcoupon) == true && ismissing(lastcoupon) == true
        date = maturity
        m = 1
        while date > settlement
            pushfirst!(cfdvector, date)  
            date = incr(eom, maturity, - m * nmper)
            m +=1
        end
        m -=1
        for i in 1:periods_before_settlement
            pushfirst!(cfdvector, date)
            date = incr(eom, maturity, - (m+i) * nmper)
        end
        for j in 1:periods_after_maturity
            date = incr(eom, maturity, j * nmper)
            push!(cfdvector, date)
        end
        return cfdvector
    elseif ismissing(firstcoupon) == true && ismissing(lastcoupon) == false
        date = lastcoupon
        m = 1
        while date > settlement
            pushfirst!(cfdvector, date)
            date = incr(eom, lastcoupon, -m * nmper)
            m += 1
        end
        m -=1
        for i in 1:periods_before_settlement
            pushfirst!(cfdvector, date)
            date = incr(eom, lastcoupon, -(m+i) * nmper)
        end
        # Case when some quasi-coupon(s) fall between lastcoupon and maturity
        newdate = incr(eom, lastcoupon, nmper)
        m = 2
        while newdate < maturity
            push!(cfdvector, newdate)
            newdate = incr(eom, lastcoupon, m * nmper)
            m += 1
        end
        m -= 1
        for j in 1:periods_after_maturity
            push!(cfdvector, newdate)
            newdate = incr(eom, lastcoupon, (j+m) * nmper)
        end
        return cfdvector
    elseif ismissing(firstcoupon) == false # && ismissing(lastcoupon) == true
        date = firstcoupon
        m = 1
        while date < settlement
            date = incr(eom, firstcoupon, m * nmper)
            m += 1
        end
        while date < maturity
            push!(cfdvector, date)
            date = incr(eom, firstcoupon, m * nmper)
            m += 1
        end
        m -= 2            
        for j in 1:periods_after_maturity
            date = incr(eom, firstcoupon, (j+m) * nmper)
            push!(cfdvector, date)
        end    

        if settlement <= firstcoupon
            for i in 1:periods_before_settlement
                date = incr(eom, firstcoupon, -i * nmper)
                pushfirst!(cfdvector, date)
            end
        else 
            datebs = firstcoupon
            z = 0
            while datebs < settlement  # Finding z that defines the date before settlement
                z += 1
                datebs = incr(eom, firstcoupon, z * nmper)
            end
            for i in 1:periods_before_settlement
                date = incr(eom, firstcoupon, (z-i) * nmper)
                pushfirst!(cfdvector, date)
            end 
        end
        return cfdvector
    end     
end

cfdatesq(settlement::Date, maturity::Date, period::Int) = cfdatesq(settlement::Date, maturity::Date, period::Int, true, missing, missing, 0, 0)

cfdatesq(settlement::Date, maturity::Date, period::Int, eom::Bool) = cfdatesq(settlement::Date, maturity::Date, period::Int, eom::Bool, missing, missing, 0, 0)

cfdatesq(settlement::Date, maturity::Date, period::Int, eom::Bool, firstcoupon, lastcoupon) = cfdatesq(settlement::Date, maturity::Date, period::Int, eom::Bool, firstcoupon, lastcoupon, 0, 0)


#=
cftimes computes the time factor of a cash flow, which is the difference between the settlement date and the cash flow date, in units of semiannual coupon periods. In computing time factors, use SIA actual/actual day count conventions for all time factor calculations.

    function cftimes(settlement::Date, maturity::Date, period::Int64, basis, endmonthrule::Bool = true, issue::Union{Date, Missing} = missing, firstcoupon::Union{Date, Missing} = missing, lastcoupon::Union{Date, Missing} = missing, face::Real) 

=#
function cftimes(settlement::Date, maturity::Date, period::Int64, basis, endmonthrule::Bool = true, issue::Union{Date, Missing} = missing, firstcoupon::Union{Date, Missing} = missing, lastcoupon::Union{Date, Missing} = missing, face::Real = 100) 
    cfvectors = cfamounts(0.0, settlement, maturity, period, basis, endmonthrule, issue, firstcoupon, lastcoupon, face)[3]
    return cfvectors[2:end]
end

#=
tmfactor determines the time factors from a vector of Settlement dates to a vector of Maturity dates

    tmfactor(seettlement, maturity)

=#
function tmfactor(settlement, maturity)
    tfvectors = cfamounts(0.0, settlement, maturity, 0, ActualActualMatlab(), false, missing, missing, missing, 100)[3]
    return tfvectors[2:end]
end

