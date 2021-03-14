#=
Certificate of Deposit (3/3):
OK  cdai	    - Accrued interest on certificate of deposit
OK  cdprice	    - Price of certificate of deposit
OK  cdyield	    - Yield on certificate of deposit (CD)
=#

#=
cdai computes the accrued interest on a certificate of deposit. cdai assumes that the certificates of deposit pay interest at maturity. Because of the simple interest treatment of these securities, this function is best used for short-term maturities (less than 1 year). The default simple interest calculation uses the Basis for the actual/360 convention 

    accrint = cdai(couponrate, settlement, maturity, issuedate, basis)

Input arguments:
    couponrate  - Annual interest rate
    settlement  - Settlement date for certificate of deposit
    maturity    - Maturity date for certificate of deposit 
    issuedate   - Issue date for certificate of deposit
    basis       - Day-count basis for certificate of deposit. Default: 2 (actual/360) 

Output arguments:
    accrint - Accrued interest per $100 of face value
=#

function cdai(couponrate, settlement, maturity, issuedate, basis)
    if settlement > maturity
        return error("must have settlement <= maturity")
    end
    return 100 * couponrate .* yearfrac(issuedate, settlement, basis)
end 
cdai(couponrate, settlement, maturity, issuedate) = cdai(couponrate, settlement, maturity, issuedate, Actual360())  # Default



#= cdprice returns the price of CD given its simple-yield. The default simple interest calculation uses the Basis for the actual/360 convention (2).
    
    price = cdprice(yield, couponrate, settlement, maturity, issuedate, basis)

Input arguments:
    yield       - simple yields to maturity, over the basis denominator.
    couponrate  - Annual interest rate
    settlement  - Settlement date for certificate of deposit
    maturity    - Maturity date for certificate of deposit 
    issuedate   - Issue date for certificate of deposit
    basis       - Day-count basis for certificate of deposit. Default: 2 (actual/360) 

Output arguments: 
    price       - clean prices of the CD per $100 face value.
=#

function cdprice(yield, couponrate, settlement, maturity, issuedate, basis)
    issuemat = yearfrac(issuedate, maturity, basis)
    settlmat = yearfrac(settlement, maturity, basis)
    accrint = couponrate .* yearfrac(issuedate, settlement, basis)
    return 100.0 * ((1.0 + couponrate .* issuemat ) ./ (1.0 + yield .* settlmat) - accrint)
end
cdprice(yield, couponrate, settlement, maturity, issuedate) = cdprice(yield, couponrate, settlement, maturity, issuedate, Actual360())

#=
cdyield(Price,CouponRate,Settle,Maturity,IssueDate) computes the yield to maturity of a certificate of deposit given its clean price. cdyield assumes that the certificates of deposit pay interest at maturity. Because of the simple interest treatment of these securities, this function is best used for short-term maturities (less than 1 year). The default simple interest calculation uses the Basis for the actual/360 convention (2).

    yield = cdyield(price, couponrate, settlement, maturity, issuedate, basis)

Input arguments:
    price       - Clean price of certificate of deposit per $100 face
    couponrate  - Annual interest rate
    settlement  - Settlement date for certificate of deposit
    maturity    - Maturity date for certificate of deposit 
    issuedate   - Issue date for certificate of deposit
    basis       - Day-count basis for certificate of deposit. Default: 2 (actual/360) 

Output arguments: 
    yield       - Simple yield to maturity of certificate of deposit
=#
function cdyield(price, couponrate, settlement, maturity, issuedate, basis)
    issuemat = yearfrac(issuedate, maturity, basis)
    settlmat = yearfrac(settlement, maturity, basis)
    accrint = couponrate .* yearfrac(issuedate, settlement, basis)
    dirty = price/100 + accrint
    return ((1.0 + couponrate .* issuemat) ./ dirty  - 1.0) .* (1.0 / settlmat)
end
cdyield(price, couponrate, settlement, maturity, issuedate) = cdyield(price, couponrate, settlement, maturity, issuedate, Actual360())