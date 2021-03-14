module Iscalc

using Dates 
using BusinessDays 
using Roots                                 # rates_of_return.jl, fixed-income_securities.jl 
using StatsBase                             # data_transformation.jl
using Distributions                         # price_derivatives.jl
using TimeSeries                            # timetables.jl
using AlphaVantage                          # financial_data.jl
using CSV, HTTP                             # financial_data.jl
using JSON                                  # financial_data.jl
using DelimitedFiles                        # From stdlib
using Statistics                            # From stdlib
using DataFrames                            # financial_data.jl
using StatsPlots
using ShiftedArrays
using JuMP                                  # portfolio_optimisation.jl
using Ipopt                                 # portfolio_optimisation.jl  (LP, QP, NLP)          / EPL
using Juniper                               # portfolio_optimisation.jl  (MISOCP, MINLP)        / MIT
using Cbc                                   # portfolio_optimisation.jl  (MILP)                 / EPL
using SCS                                   # portfolio_optimisation.jl  (LP, SOCP, SDP)        / MIT
using LaTeXStrings                          # portfolio.optimisation.jl



import ForwardDiff: Dual                    # price_derivatives.jl
import BusinessDays: TARGET, USNYSE, CanadaTSX, AustraliaASX, BrazilB3, USSettlement, USGovernmentBond, UKSettlement, CanadaSettlement, BRSettlement, WeekendsOnly       # daycounts_extra_tests.jl
import SpecialFunctions: erfc               # investment_performance.jl
import Interpolations: LinearInterpolation  # investment_performance.jl
import DSP: filt                            # stats.jl
import ShiftedArrays: lag, lead
import LinearAlgebra: diag, Diagonal, dot   # portfolio.optimisation.jl

export Actual365Fixed, Actual360, ActualActualISDA, Thirty360, ThirtyE360, ThirtyE360ISDA, ActualActualICMA, Actual365Fixed                                # daycounts_kernel.jl
export Thirty360Excel, ActualActualExcel                        # daycounts_excel.jl
export NL365, Bus252, Thirty360PSA, Thirty360SIA, ThirtyEPlus360, ActualActualL, Actual365A, ActualActualY, Thirty360ICMA, Actual365ICMA, Actual360ICMA, ActualActualAFB    # daycounts_extra.jl
export ActualActualMatlab                                       # daycounts_matlab.jl
export Robust, Secant, Steffensen, Brent                        # rates_of_return.jl, price_derivatives.jl

export acrubond, acrudisc # accrued_interest.jl (2/2)
export annurate, annuterm, payadv, payper, payuni, payodd # anuities.jl (6/6)
export amortize, depfixd, deprdv, depstln, depsoyd, depgendb # amortization_depreciation.jl (6/6)
export days252bus, days360, days360e, days360isda, days360psa, days365, daysact, daysdif, busdate, datemnth, datewrkdy, daysadd, fbusdate, isbusday, lbusdate, quarter, thirdwednesday, wrkdydif # business_days.jl (19/19)
export cfdates, cfdatesq, cfamounts, cftimes, tmfactor # cash-flow_generation (5/6)
export cfdur, cfconv # cash-flow_sensitivities.jl (2/4)
export cdai, cdprice, cdyield # certificates_deposit.jl (3/3)
export cpndatenq, cpndatepq, cpndaysn, cpndaysp, cpnpersz, cpndatep, cpndaten, cpncount, accrfrac # coupon_dates.jl (9/9)
export cur2frac, cur2str, frac2cur, thirtytwo2dec, todecimal, toquoted # currency_price_conversion.jl (6/7)
export holdings2weights, weights2holdings, cor2cov, cov2cor, tick2ret, ret2tick, boxcox, arith2geom, geom2arith # data_transformation (9/11)
export today, second, minute, hour, day, month, year, eomdate, yeardays, lweekdate, months, nweekdate, weeknum # date_time.jl (13/13) 
export num29feb, eomrule # date_utils.jl (4/4)
export yearfrac                                 # daycounts_kernel.jl, daycounts_excel.jl, daycounts_extra.jl
export issia, isisma                            # daycounts_extra.jl
export stock, forex, report, fred, yahoo, alphavantage            # financial_data (5/5)
export beytbill, discrate, ylddisc, yldtbill, yldmat, tbl2bond, bndyield # fixed-income_securities.jl (7/7)
export transbasis, datenum, datestr # interopt.jl (1/4)
export maxdrawdown, emaxdrawdown, lpm, elp, inforatio, portalpha, sharpe, omega, sortino, upside # investment_performance.jl (10/10)
export pvfix, pvvar, fvdisc, fvfix, fvvar, npv, xnpv, nfv, xnfv # present_future_value.jl (5/5)
export prdisc, prmat, prtbill, bndprice, floatmargin # price_bonds.jl (5/8)
export opprofit, blsprice, blkprice, blsdelta, blsgamma, blsvega, blsrho, blstheta, blslambda, blsimpv, blkimpv, blspsi, blsvanna, binprice, binprice_eu, bingreeks, mcprice, mcgreeks # price_derivatives.jl (12/12)
export effrr, irr, mirr, nomrr, taxedrr, xirr # rates_of_return.jl (6/6)
export movavg, moving, movstd, movsum, movmean, movmax, movmin, movmedian, movmad, movcor, movslope  # stats.jl (11/11)
export timetable, load # timetables.jl
export adosc, adline, chaikocs, bollinger, onbalvol, chaikvolat, stochosc, tsaccel, tsmom, willad, ema, sma, macd, rsindex, medprice, willpctr, llow, hhigh, pvtrend,  posvolidx, negvolidx, wclose, prcroc, volroc   # technical_indicators.jl (24/24)
export tbilldisc2yield, tbillyield2disc, tbillprice, tbillyield, tbillrepo, tbillval01  # treasury_bills.jl (6/6)
export disc2zero, fwd2zero, zero2disc, zero2fwd # yield-curves.jl (4/9)

export Portfolio, setAssetMoments, setAssetList, setInitPort, setCosts, checkFeasibility
export setDefaultConstraints, setBounds, setBoundsType, setBudget, setTurnover, setOneWayTurnover, setNetReturnTarget, setTrackingPort, setTrackingError, setMinMaxNumAssets, addEquality, addInequality, addGroups, addGroupRatio
export optimalp, estimateMinimumVariancePortfolio, estimateMaximumNetReturnPortfolio, estimateFrontierLimits, estimateFrontierByReturn, estimateFrontierByRisk, estimateFrontier, estimateMaxSharpeRatio
export plotFrontier, plotFrontierAllocations, plotSharpeRatio
export computePortSharpeRatio, computePortMoments, computePortReturn, computePortRisk, portfolioDiff


abstract type DayCount end
    
# daycounts_kernel.jl
struct Actual365Fixed <: DayCount end
struct Actual360 <: DayCount end
struct ActualActualISDA <: DayCount end
struct ActualActualICMA{T<:AbstractVector{Date}} <: DayCount
    frequency::Int
    schedule::T
end
struct Thirty360 <: DayCount end
struct ThirtyE360 <: DayCount end
struct ThirtyE360ISDA <: DayCount  
    maturity::Date
end
# daycounts_excel.jl
struct Thirty360Excel <: DayCount end
struct ActualActualExcel <: DayCount end
# dayconts_extra.jl
struct NL365 <: DayCount end
struct Bus252 <: DayCount
    calendar::HolidayCalendar
end
struct Thirty360PSA <: DayCount end
struct Thirty360SIA <: DayCount
    eom::Bool
end
struct ThirtyEPlus360 <: DayCount end
struct ActualActualL <: DayCount
    frequency::Int
 end
 struct Actual365A <: DayCount end
 struct ActualActualY <: DayCount end
 struct Actual365ICMA <: DayCount end
 struct Actual360ICMA <: DayCount end
 struct Thirty360ICMA <: DayCount end
 struct ThirtyE360ISMA <: DayCount end
 struct ActualActualAFB <: DayCount end 
# daycounts_matlab.jl
 struct ActualActualMatlab <: DayCount end

 
 
 abstract type DateRoll end

 struct Unadjusted <: DateRoll; end
 struct Following <: DateRoll; end
 struct ModifiedFollowing <: DateRoll; end
 struct Preceding <: DateRoll; end
 struct ModifiedPreceding <: DateRoll; end
    


abstract type RootFindingMethod end

struct Robust <: RootFindingMethod 
    x0::Float64
end
struct Secant <: RootFindingMethod          # Secant Method
    x0::Float64
end
struct Steffensen <: RootFindingMethod      # Steffensen Method
    x0::Float64
end
struct Brent <: RootFindingMethod           # Brent Method
    inflimit::Float64
    suplimit::Float64
end




abstract type BinomialTree end

struct CRR <: BinomialTree              # Cox-Ross-Rubinstein model
    steps::Int64
end
struct CRRd <: BinomialTree             # Cox-Ross-Rubinstein model with drift
    steps::Int64
end
struct JRep <: BinomialTree             # Jarrow-Rudd (equal-probability) model
    steps::Int64
end
struct JRrn <: BinomialTree             # Jarrow-Rudd (risk-neutral) model
    steps::Int64
end
struct Tian <: BinomialTree             # Tian model
    steps::Int64
end
struct LR <: BinomialTree               # Leisen-Reimer model
    steps::Int64
end

abstract type MonteCarlo end

struct LS <: MonteCarlo                 # Longstaff–Schwartz least-squares method
    steps:: Int64
    paths:: Int64
end



abstract type AbstractPortfolio end

mutable struct Portfolio <: AbstractPortfolio
# Set Up    
    AssetList:: Vector{Symbol}      # Names or symbols of assets in universe
    InitPort:: Vector{Float64}      # Initial portfolio
    Name:: String                   # Name for instance of Portfolio object
    NumAssets:: Real                # Number of assets in the universe

# Portfolio Forecasting
    AssetMean      # Mean of asset returns
    AssetCovar     # Covariance of asset returns (square matrix)
    RiskFreeRate   # Risk-free rate
    BuyCost        # Proportional purchase costs
    SellCost       # Proportional sales costs

# Portfolio Constraints
    # Bound Constraints
    LowerBound     # Lower-bound constraint, e.g zeros(NumAssets) if portfolio weights must be nonnegative
    UpperBound     # Upper-bound constraint
    BoundType      # : simple or :conditional. A Conditional BoundType is also known as semicontinuous constraint. For example, the weight you invest in each asset is either 0 or between [0.01, 0.5]

    # Budget Constraints
    LowerBudget     # Lower-bound budget constraint. 1 if portfolio weights must sum to 1.
    UpperBudget     # Upper-bound budget constraint. 1 if portfolio weights must sum to 1.
  
    # Net Return Target Constraint
    NetReturn       # Net return inequality constraint.
    
    # Cardinality Constraints
    MinNumAssets    # Cardinality constraints
    MaxNumAssets    # Cardinality constraints

    # Turnover Constraints
    BuyTurnover    # Turnover constraint on purchases
    SellTurnover   # Turnover constraint on sales
    Turnover       # Turnover constraint

    # Tracking Error Constraints
    TrackingPort   # Tracking portfolio for tracking error constraint
    TrackingError  # Upper bound for tracking error constraint

    # Linear Constraints
    AEquality      # Linear equality constraint matrix
    bEquality      # Linear equality constraint vector
    AInequality    # Linear inequality constraint matrix
    bInequality    # Linear inequality constraint vector

    # Quadratic Constraints
    AQuadEquality       # Quadratic equality constraint matrix
    bQuadEquality       # Quadratic equality constraint vector
    AQuadInequality     # Quadratic inequality constraint matrix
    bQuadInequality     # Quadratic inequality constraint vector

    # Group Constraints
    GroupMatrix     # Group membership matrix
    GroupA          # Group A weights to be bounded by weights in group B
    GroupB          # Group B weights 
    UpperGroup      # Upper-bound group constraint
    LowerGroup      # Lower-bound group constraint
    UpperRatio      # Maximum ratio of allocations between Groups A and B 
    LowerRatio      # Minimum ratio of allocations between Groups A and B 
end


mutable struct PortfolioCVar <: AbstractPortfolio
# Set Up    
    AssetList:: Vector{Symbol}      # Names or symbols of assets in universe
    InitPort:: Vector{Float64}      # Initial portfolio
    Name:: String                   # Name for instance of Portfolio object
    NumAssets:: Real                # Number of assets in the universe

# Portfolio Modeling
    RiskFreeRate        # Risk-free rate
    BuyCost             # Proportional purchase costs
    SellCost            # Proportional sales costs
    ProbabilityLevel    # Value-at-risk probability level which is 1 − (loss probability) 
    NumScenarios        # Number of Scenarios

# Portfolio Constraints
    # Bound Constraints
    LowerBound     # Lower-bound constraint, e.g zeros(NumAssets) if portfolio weights must be nonnegative
    UpperBound     # Upper-bound constraint
    BoundType      # : simple or :conditional. A Conditional BoundType is also known as semicontinuous constraint. For example, the weight you invest in each asset is either 0 or between [0.01, 0.5]

    # Budget Constraints
    LowerBudget    # Lower-bound budget constraint. 1 if portfolio weights must sum to 1.
    UpperBudget    # Upper-bound budget constraint. 1 if portfolio weights must sum to 1.
    
    # Cardinality Constraints
    MinNumAssets   # Cardinality constraints
    MaxNumAssets   # Cardinality constraints

    # Turnover Constraints
    BuyTurnover    # Turnover constraint on purchases
    SellTurnover   # Turnover constraint on sales
    Turnover       # Turnover constraint

    # Linear Constraints
    AEquality      # Linear equality constraint matrix
    bEquality      # Linear equality constraint vector
    AInequality    # Linear inequality constraint matrix
    bInequality    # Linear inequality constraint vector

    # Quadratic Constraints
    AQuadEquality       # Quadratic equality constraint matrix
    bQuadEquality       # Quadratic equality constraint vector
    AQuadInequality     # Quadratic inequality constraint matrix
    bQuadInequality     # Quadratic inequality constraint vector
        
    # Group Constraints
    GroupMatrix     # Group membership matrix
    GroupA          # Group A weights to be bounded by weights in group B
    GroupB          # Group B weights 
    UpperGroup      # Upper-bound group constraint
    LowerGroup      # Lower-bound group constraint
    UpperRatio      # Maximum ratio of allocations between Groups A and B 
    LowerRatio      # Minimum ratio of allocations between Groups A and B 
end


include("daycounts_kernel.jl")
include("daycounts_excel.jl")
include("daycounts_extra.jl")
include("daycounts_matlab.jl")

include("accrued_interest.jl")
include("amortization_depreciation.jl")
include("annuities.jl") 
include("business_days.jl")
include("cash-flow_generation.jl")
include("cash-flow_sensitivities.jl")
include("certificates_deposit.jl")
include("coupon_dates.jl")
include("currency_price_conversion.jl")
include("data_transformation.jl")
include("date_time.jl")
include("date_utils.jl")
include("financial_data.jl")
include("fixed-income_securities.jl")
include("interopt.jl")
include("investment_performance.jl")
include("portfolio_optimisation.jl")
include("present_future_value.jl")
include("price_bonds.jl")
include("price_derivatives.jl")
include("rates_of_return.jl")
include("stats.jl")
include("technical_indicators.jl")
include("timearrays.jl")
include("treasury_bills.jl")
include("yield-curves.jl")

include("aliases.jl")


end # module
