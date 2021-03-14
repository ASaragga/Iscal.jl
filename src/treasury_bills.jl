#=
Treasury Bills (6/6):
OK  tbilldisc2yield - Convert Treasury bill discount to equivalent yield
OK  tbillprice	    - Price Treasury bill
OK  tbillrepo	    - Break-even discount of repurchase agreement
OK  tbillval01	    - Value of one basis point
OK  tbillyield	    - Yield on Treasury bill
OK  tbillyield2disc - Convert Treasury bill yield to equivalent discount
=#

#=
tbilldisc2yield converts the discount rate on Treasury bills into their respective money-market or bond-equivalent yields.

    (BEYield, MMYield) = tbilldisc2yield(Discount,Settle,Maturity)

Input arguments:
    discount    - Discount rates of T-Bills
    settlement  - Settlement dates
    maturity    - Maturity of T-Bills
Output arguments:
    BEYield     - Bond equivalent yields of the Treasury bills. 
    MMYield     - Money-market yields of the Treasury bills. 

Note: The Money Market Yield basis is actual/360.
      The Bond-Equivalent Yield basis is actual/365.
      The Discount Rate basis is actual/360.
=#
function tbilldisc2yield(Discount, Settle, Maturity)
    # find "short" and "long" bills, ones with maturity shorter/equal to and longer than 182 days, respectively on 365 day basis
    DSM = daysact(Settle, Maturity)
    if DSM <= 182
        idx = :short 
    elseif DSM > 182
        idx = :long
    end

    if idx == :short  # Short T-Bills 
        MMYield = 360 * Discount / (360 - Discount * DSM)
        BEYield = MMYield * 365/360 
    end

    if idx == :long  # Long T-Bills - all qty are annualized
        Price = 100*(1 - Discount .* DSM/360)    
        MMYield = ((100 / Price)-1) * (360 / DSM);
        X = DSM/365
        BEYield = (-2*X + 2*sqrt( X^2 - (2*X - 1)* (1 - 100 / Price))) / (2*X - 1)
    end
    return BEYield, MMYield
end

#=
tbillyield2disc converts the yield on some Treasury bills into their respective discount rates.

    Discount = tbillyield2disc(Yield, Settle, Maturity, Type)

Input arguments:
    Type — Yield type (determines how to interpret values entered in Yield), with values:  
        :mm = money market (default) 
        :be = bond-equivalent
Output arguments:
    Discount - Discount rates of the Treasury bills,
=#
function tbillyield2disc(Yield, Settle, Maturity, Type = :mm)
    # find "short" and "long" bills, ones with maturity shorter/equal to and longer than 182 days, respectively on 365 day basis
    DSM = daysact(Settle, Maturity)
    if DSM <= 182
        idx = :short 
    elseif DSM > 182
        idx = :long
    end

    if Type == :mm
        A = 360  # US Money-market Yield basis
    elseif Type == :be
        A = 365  # US Bond Equivalent Yield basis
    end

    if idx == :short  # Short T-bills (maturity < 182 days)
        Discount = 360/DSM * (1 - (1 / (1 + Yield * DSM/A)))
    end
    if idx == :long  # Long T-Bills - all qty are annualized
        Discount = 360/DSM * (1 - (1 /(((1 + Yield/2)) * (1 +(2 * DSM/A - 1) * Yield/2))))
    end
    return Discount
end


#=
tbillprice computes the price of a Treasury bill given a yield or discount rate.

    Price = tbillprice(Rate, Settle, Maturity, Type) 

Input arguments:
    Rate — Bond-equivalent yield, money-market yield, or discount rate (defined by the input Type)
    Type - Rate type (determines how to interpret values entered in Rate), specified as:
        :mm  = money market yield 
        :be  = bond-equivalent yield (default)
        :dsc = discount rate
    Note: The bond-equivalent yield basis is actual/365. The money-market yield basis is actual/360. The discount rate basis is actual/360.
=#
function tbillprice(Rate, Settle, Maturity, Type) 
    if Type == :mm || Type == :be # First convert the MMYield and BEYield to Discount Rate
        Rate = tbillyield2disc(Rate,Settle,Maturity,Type)
    end
    # And compute the Price of T-bill now that we have all annualized discount rate.
    Face = 100
    DSM = daysact(Settle, Maturity)
    Price = Face - (Rate * Face * DSM/360)
    return Price
end

#=
tbillyield computes the yield of US Treasury bills given Price, Settle, and Maturity.

    (MMYield, BEYield, Discount) = tbillyield(Price, Settle, Maturity)

Input arguments:
    Price — Price of Treasury bills for every $100 face value
Output arguments:
    MMYield — Money-market yields of Treasury bills 
    BEYield — Bond equivalent yields of Treasury bills
    Discount — Discount rates of Treasury bills
=#
function tbillyield(Price, Settle, Maturity)
    DSM = daysact(Settle, Maturity)
    if DSM <= 182
        idx = :short 
    elseif DSM > 182
        idx = :long
    end
    MMYield = (100/Price - 1) * 360/DSM
    Discount = (1 - Price/100) * 360/DSM
    
    # short T-Bills
    if idx == :short
        BEYield = MMYield * 365/360    
    end
    
    # long T-Bills
    if idx == :long
        X = DSM/365
        BEYield = (-2 * X + 2*sqrt( X.^2 -(2*X - 1) * (1 - 100/Price))) / (2*X - 1)
    end
    return MMYield, BEYield, Discount
end

#= 
tbillrepo computes the true break-even discount of a repurchase agreement.

    TBEDiscount = tbillrepo(RepoRate, InitialDiscount, PurchaseDate, SaleDate, Maturity)

Input arguments:
    RepoRate        — Annualized, 360-day based repurchase rate. ACT/360 based.
    InitialDiscount - Discount on the Treasury bill on the day of purchase
    PurchaseDate    - Date the Treasury bill is purchased
    SaleDate        — Date Treasury bill repurchase term is due
    Maturity        — Maturity date of Treasury bill
Output arguments:
    TBEDiscount     — True break-even discount of repurchase agreement
=#
function tbillrepo(RepoRate, InitialDiscount, PurchaseDate, SaleDate, Maturity)
    T1 = daysdif(PurchaseDate, Maturity, ActualActualMatlab())
    T2 = daysdif(SaleDate, Maturity, ActualActualMatlab())
    tbediscount = (1 - (1 - InitialDiscount * T1/360) * (1 + RepoRate * (T1-T2)/360) ) * 360/T2
    return tbediscount
end

#=
tbillval01 calculates the value of one basis point of $100 Treasury bill face value on the discount rate, money-market yield, or bond-equivalent yield.

    (Val01Disc, Val01MMY, Val01BEY) = tbillval01(Settle, Maturity)

Output arguments:
    Val01Disc   — Value of one basis point of discount rate for every $100 face
    Val01MMY    — Value of one basis point of money-market yield for every $100 face
    Val01BEY    — Value of one basis point of bond-equivalent yield for every $100 face
=#
function tbillval01(Settle, Maturity)
    # Days to maturity in actual basis
    DSM = daysact(Settle, Maturity)
    
    # Rate = 0.0001 (1 bp.) out of $100
    Rate = 0.0001
    Face = 100
    
    CMMY = 360 * Rate / (360 + Rate * DSM)
    CBEY = 360 * Rate / (365 + Rate *DSM)
    
    # We will compute the value of Tbills 01 now that we have all annualized discount rate
    value01Disc = Rate * Face * DSM/360
    value01MMY = CMMY * Face * DSM/360
    value01BEY = CBEY * Face * DSM/360
    
    return value01Disc, value01MMY, value01BEY 
end
