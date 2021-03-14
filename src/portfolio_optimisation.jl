#= 
portfolio_optimization (43/47)

## 4.1 - Portfolio Optimization Theory (2/3)
OK  Portfolio       - Create object for mean-variance portfolio optimization and analysis
OK  PortfolioCVaR   - Creates object for conditional value-at-risk portfolio optimization and analysis
    PortfolioMAD    - Create  object for mean-absolute deviation portfolio optimization and analysis


### 4.2 - Create Portfolio (3/3) 5/6
OK  setAssetList    - Set up list of identifiers for assets
OK  setInitPort	    - Set up initial or current portfolio
OK  setCosts        - Set up proportional transaction costs


### 4.3 - Estimate Mean and Covariance for Returns (3/3) 8/9
OK  estimateAssetMoments    - Estimate mean and covariance of asset returns from data
OK  setAssetMoments	        - Set moments (mean and covariance) of asset returns for Portfolio object
OK  getAssetMoments	        - Obtain mean and covariance of asset returns from Portfolio object


### 4.4 - Set Portfolio Constraints (12/12) 20/21
OK  setDefaultConstraints	- Set up portfolio constraints with nonnegative weights that sum to 1
OK  setBounds	            - Set up bounds for portfolio weights
OK  setBudget	            - Set up budget constraints
OK  setTurnover             - Set up maximum portfolio turnover constraint
OK  setOneWayTurnover	    - Set up one-way portfolio turnover constraints
OK  setTrackingPort         - Set up benchmark portfolio for tracking error constraint
OK  setTrackingError        - Set up maximum portfolio tracking error constraint
OK  setMinMaxNumAssets      - Set cardinality constraints on the number of assets invested
OK  setEquality             - Set up linear equality constraints for portfolio weights
OK  setInequality	        - Set up linear inequality constraints for portfolio weights
OK  setGroups	            - Set up group constraints for portfolio weights
OK  setGroupRatio	        - Set up group ratio constraints for portfolio weights


### 4.5 - Add Portfolio Constraints (4/4) 24/25
OK  addEquality     - Add linear equality constraints for portfolio weights to existing constraints
OK  addInequality   - Add linear inequality constraints for portfolio weights to existing constraint
OK  addGroups       - Add group constraints for portfolio weights to existing group constraints
OK  addGroupRatio   - Add group ratio constraints for portfolio weights to existing group ratio constraints
    
    
### 4.6 - Display Portfolio Constraints (8/8) 32/33
OK  getBounds	        - Obtain bounds for portfolio weights 
OK  getBudget	        - Obtain budget constraint bounds 
OK  getCosts            - Obtain buy and sell transaction costs 
OK  getEquality         - Obtain equality constraint arrays 
OK  getGroupRatio	    - Obtain group ratio constraint arrays 
OK  getGroups	        - Obtain group constraint arrays 
OK  getInequality	    - Obtain inequality constraint arrays 
OK  getOneWayTurnover	- Obtain one-way turnover constraints 


### 4.7 - Validate Portfolio (1/2) 33/35
OK  checkFeasibility    - Check feasibility of input portfolios against input portfolio object.
    estimateBounds      - Estimate global lower and upper bounds for set of portfolios


### 4.8 - Estimate Efficient Frontiers and Portfolios

#### Efficient Frontiers (5/5) 38/40
OK  estimateFrontierLimits      - Estimate optimal portfolios at endpoints of efficient frontier
OK  estimateFrontierByReturn    - Estimate optimal portfolios with targeted portfolio returns
OK  estimateFrontierByRisk      - Estimate optimal portfolios with targeted portfolio risks
OK  estimateFrontier            - Estimate specified number of optimal portfolios on the efficient frontier
OK  plotFrontier                - Plot efficient frontier
===============================
OK plotFrontierAllocations       - Plot stacked bars with efficient frontier allocations
OK plotSharpeRatio              - Plot maximum Sharpe ratio portfolio

#### Efficient Portfolios (4/5) 42/45
    estimateMaxSharpeRatio  - Estimate efficient portfolio to maximize Sharpe ratio
OK  computePortSharpeRatio	- Estimate Sharpe ratio of given portfolio weights
OK  computePortMoments	    - Estimate moments of portfolio returns
OK  computePortReturn      - Estimate mean of portfolio returns
OK  computePortRisk        - Estimate portfolio risk according to risk proxy associated with corresponding object


### 4.9 - Set Solver (1/2) 43/47
    setSolver	    - Choose main solver and specify associated solver options for portfolio optimization
OK  setSolverMINLP  - Choose mixed integer nonlinear programming (MINLP) solver for portfolio optimization


### 4.10 - Utilities* (not in MATLAB)
OK  portfolioDiff   - Shows purchases and sales relative to initial portfolio.
OK  optimalp        - Estimates optimal portfolios with (i) targeted portfolio risks or returns, or (ii) a specified number of optimal portfolios on the efficient frontier, including endpoints

=#


#=
Create Portfolio object for mean-variance portfolio optimization and analysis

    Portfolio()

=#
function Portfolio()
    Portfolio(
    # Set Up
        Symbol[],   # Names or symbols of assets in universe
        Float64[],  # Initial portfolio
        "",         # Name for instance of Portfolio object
        NaN,        # Number of assets in the universe
    # Portfolio Forecasting
        Float64[],  # Mean of asset returns
        Array{Float64}(undef, 0, 0),    # Covariance of asset returns (square matrix)
        NaN,        # Risk-free rate
        Float64[],  # Proportional purchase costs
        Float64[],  # Proportional sales costs
    # Portfolio Constraints
        # Bound Constraints
        Float64[],  # Lower-bound constraint, e.g zeros(NumAssets) if portfolio weights must be nonnegative
        Float64[],  # Upper-bound constraint
        # Semicontinuous Constraints
        Symbol[],
        # Budget Constraints
        NaN,  # Lower-bound budget constraint. 1 if portfolio weights must sum to 1.
        NaN,  # Upper-bound budget constraint. 1 if portfolio weights must sum to 1.
        # Net Return Target Constraint
        NaN,  # Taget net return
        # Cardinality Constraints
        NaN,  # Cardinality constraints
        NaN,  # Cardinality constraints
        # Turnover Constraints
        NaN,  # Turnover constraint on purchases
        NaN,  # Turnover constraint on sales
        NaN,  # Turnover constraint
        # Tracking Error Constraints
        Float64[],  # Tracking portfolio for tracking error constraint
        NaN,        # Upper bound for tracking error constraint
        # Linear Constraints
        Array{Float64}(undef, 0, 0),    # Linear equality constraint matrix
        Float64[],                      # Linear equality constraint vector
        Array{Float64}(undef, 0, 0),    # Linear inequality constraint matrix
        Float64[],                      # Linear inequality constraint vector
        # Quadratic Constraints
        Array{Float64}(undef, 0, 0, 0),     # Quadratic equality constraint matrix
        Float64[],                          # Quadratic equality constraint vector
        Array{Float64}(undef, 0, 0, 0),     # Quadratic inequality constraint matrix
        Float64[],                          # Quadratic inequality constraint vector
        # Group Constraints
        Array{Bool}(undef, 0, 0),   # Group membership matrix 
        Array{Bool}(undef, 0, 0),   # Group A weights to be bounded by weights in group B
        Array{Bool}(undef, 0, 0),   # Group B weights 
        Float64[],      # Upper-bound group constraint
        Float64[],      # Lower-bound group constraint
        Float64[],      # Maximum ratio of allocations between Groups A and B 
        Float64[]       # Minimum ratio of allocations between Groups A and B
    )
end 


#=
Creates PortfolioCVaR object for conditional value-at-risk portfolio optimization and analysis

    PortfolioCVar()

=#


#=
setAssetList sets up list of identifiers for assets

    p = setAssetList(p, A::Vector{Symbol})


If an asset list is entered as an input, this function overwrites an existing asset list in the Portfolio if one exists. If no asset list is entered as an input, three actions can occur:
    If NumAssets is nonempty and AssetList is empty, AssetList becomes a numbered list of assets with default names according to the hidden property in defaultforAssetList ('Asset').
    If NumAssets is nonempty and AssetList is nonempty, nothing happens.
    If NumAssets is empty and AssetList is empty, the default NumAssets =1 is set and a default asset list is created ('Asset1').
=#
function setAssetList(p, defaultforAssetList::Symbol = :Asset)
    if isnan(p.NumAssets) == false && p.AssetList == Symbol[]
        p.AssetList = Array{Symbol}(undef, p.NumAssets)
        for i =1:p.NumAssets
            p.AssetList[i] = Symbol(String(defaultforAssetList) * "$i")
        end
    end
    return p
end
function setAssetList(p, A::Vector{Symbol})
    p.AssetList = A
    if isnan(p.NumAssets) == true
        p.NumAssets = length(A)
    end
    return p
end

#=
setInitPort	sets up initial or current portfolio

    p = setInitPort(p, A::Vector{Float64})

=#
function setInitPort(p, A::Vector{Float64})
    p.InitPort = A
    if isnan(p.NumAssets) == true
        p.NumAssets = length(A)
    end
    return p
end


#=
setCosts sets up proportional transaction costs

    p = setCosts(p, BC::Vector{Float64}, SC::Vector{Float64})

Input arguments:
    BC - Buy costs
    SC - Sell costs
Note: Given proportional transaction costs and an initial portfolio in the variables BuyCost, SellCost, and InitPort, the transaction costs for any portfolio Port reduce expected portfolio return by:
    BuyCost' * max{0, Port - InitPort} + SellCost' * max{0, InitPort - Port}
=#
function setCosts(p, BuyCost::Vector{Float64}, SellCost::Vector{Float64})
    if length(BuyCost) != length(SellCost)
        error("The vectors Buy Cost and Sell Cost must be of the same length")
    end
    p.BuyCost = BuyCost
    p.SellCost = SellCost
    lBC = length(BuyCost)
    if isnan(p.NumAssets) == true
        p.NumAssets = lBC
    elseif lBC != p.NumAssets
        error("BuyCost ($lBC) is not conformable with NumAssets ($(p.NumAssets)")
    end
    return p
end
function setCosts(p, BuyCost::Vector{Float64}, SellCost::Vector{Float64}, InitPort::Vector{Float64})
    if length(BuyCost) != length(SellCost)
        error("The vectors Buy Cost and Sell Cost must be of the same length")
    end
    p.BuyCost = BuyCost
    p.SellCost = SellCost
    p.InitPort = InitPort
    lBC = length(BuyCost)
    if isnan(p.NumAssets) == true
        p.NumAssets = lBC
    elseif lBC != p.NumAssets
        error("BuyCost ($lBC) is not conformable with NumAssets ($(p.NumAssets)")
    end
    return p
end


#= 
setNetReturnTarget sets the net return target (gross return less buy and sell transaction costs) for portfolio p 

    setNetReturnTarget(p, netreturn)

Inequality constraint: Portfolio net return >= netreturn
=#
function setNetReturnTarget(p, netreturn)
    p.NetReturn = netreturn
    if p.BuyCost == Float64[]
        p.BuyCost = zeros(p.NumAssets)
    end
    if p.SellCost == Float64[]
        p.SellCost = zeros(p.NumAssets)
    end
    if p.InitPort == Float64[]
        p.InitPort = zeros(p.NumAssets)
    end
end


#=
Estimate mean and covariance of asset returns from data in the following formats: Matrix, DataFrame or TimeArray

    p = estimateAssetMoments(p, Returns)
    p = estimateAssetMoments(p, Prices, ReturnType)

Input arguments:
    p           - Portfolio 
    Returns     - Return data
    Prices      - Price data
    ReturnType  - Return type, either 
                    :log 
                    :simple

Also sets the number of assets (NumAssets)
=#
# helper functions
function estimateMomentsMatrix(p, Returns::Matrix)
    p.AssetMean = transpose(mean(Returns, dims = 1))
    p.AssetCovar = cov(Returns)
    if isnan(p.NumAssets) == true
        p.NumAssets = size(Returns,2)
    end
end
function estimateMomentsMatrix(p, Prices::Matrix, ReturnType::Symbol)
    nobs = size(Prices,1)
    ncols = size(Prices,2)
    B = zeros(nobs-1,ncols)
    if ReturnType == :log
        for i = 1:nobs-1
            B[i,:] = log.(Prices[i+1,:]./Prices[i,:])
        end
        estimateAssetMoments(p, B)
    elseif ReturnType == :simple
        for i = 1:nobs-1
            B[i,:] = Prices[i+1,:]./Prices[i,:] .- 1.0
        end
        estimateAssetMoments(p, B)
    else 
        throw(UndefVarError(Symbol("$ReturnType")))
    end
    if isnan(p.NumAssets) == true
        p.NumAssets = size(Prices,2)
    end
end
function estimateAssetMoments(p, Returns)
    if isa(Returns, Matrix) == true
        estimateMomentsMatrix(p, Returns)
    elseif isa(Returns, DataFrames.DataFrame) == true
        estimateMomentsMatrix(p, Matrix(Returns))
    elseif isa(Returns, TimeSeries.TimeArray) == true
        estimateMomentsMatrix(p, values(Returns))
    end
    return p
end
function estimateAssetMoments(p, Prices, ReturnType::Symbol)
    if isa(Prices, Matrix) == true
        estimateMomentsMatrix(p, Prices, ReturnType)
    elseif isa(Prices, DataFrames.DataFrame) == true
        estimateMomentsMatrix(p, Matrix(Prices), ReturnType)
    elseif isa(Prices, TimeSeries.TimeArray) == true
        estimateMomentsMatrix(p, values(Prices), ReturnType)
    end
    return p
end


#=
setAssetMoments sets the mean and covariance of asset returns for portfolio p

    p = setAssetMoments(p,AssetMean,AssetCovar)
    p = setAssetMoments(p, AssetMean, AssetStd, AssetCorrel)

Also sets the number of assets (NumAssets)
=#
function setAssetMoments(p, AssetMean::Vector{Float64}, AssetCovar)
    p.AssetMean = AssetMean
    p.AssetCovar = AssetCovar
    if size(p.AssetMean,1) != size(p.AssetCovar,1)
        error("Number of assets in Mean vector$(size(p.AssetMean)) and Covariance matrix$(size(p.AssetCovar)) mismatch")
    end
    if isnan(p.NumAssets) == true
        p.NumAssets = length(p.AssetMean)
    end
    return p
end
function setAssetMoments(p, AssetMean::Vector{Float64}, AssetStd, AssetCorrel)
    A = Diagonal(AssetStd)
    AssetCovar = A * AssetCorrel * A
    setAssetMoments(p, AssetMean, AssetCovar)
end

#=
setBounds sets up bounds for portfolio weights for a portfolio object

    p = setBounds(p, LBound, UBound = Float64[], BoundT = :simple)

Input arguments:
    LBound - Lower bound constraints (vector or scalar that will be expanded to vector)
    UBound - Upper bound constraints (vector or scalar that will be expanded to vector)
    BoundT - Bound type (vector or scalar that will be expanded to vector). Possible values are :simple, :conditional
Given bound constraints LowerBound and UpperBound and 'Simple' BoundType, every weight in a portfolio Port must satisfy the following:
    LowerBound <= Port <= UpperBound
Given bound constraints LowerBound and UpperBound, and 'Conditional' BoundType, every weight in a portfolio Port must satisfy the following:
    Port = 0 or LowerBound <= Port <= UpperBound
If either LowerBound or UpperBound are input as empties with [], the corresponding attributes in the Portfolio, PortfolioCVaR, or PortfolioMAD object are cleared and set to [].
If LowerBound or UpperBound are specified as scalars and NumAssets exists or can be computed, then they undergo scalar expansion. The default value for NumAssets is 1.
If both LowerBound and UpperBound exist and they are not ordered correctly, the setBounds function switches bounds if necessary.
If 'Conditional'BoundType is specified, the LowerBound cannot be a negative value.
=#
function setBounds(p, LBound::Vector{<:Real}, UBound::Vector{<:Real})
    if LBound == Float64[] && UBound == Float64[]
        return p
    end
    if LBound != Float64[]
        p.LowerBound = Float64.(LBound)
        if isnan(p.NumAssets) == true 
            p.NumAssets = length(LBound)
        end
    end
    if UBound != Float64[]
        p.UpperBound = Float64.(UBound)
        if isnan(p.NumAssets) == true 
            p.NumAssets = length(UBound)
        end
    end
    return p
end

function setBounds(p, LBound::Vector{<:Real})
    setBounds(p, LBound, Float64[])
end

function setBounds(p, LBound::Real, UBound::Real, NAssets = p.NumAssets)
    p.LowerBound = Array{Real}(undef, NAssets)
    p.LowerBound .= LBound
    p.UpperBound = Array{Real}(undef, NAssets)
    p.UpperBound .= UBound 
    if isnan(p.NumAssets) == true
        p.NumAssets = NAssets
    end
    return p
end


function setBounds(p, LBound::Real, NAssets = p.NumAssets)
    p.LowerBound = Array{Real}(undef, NAssets)
    p.LowerBound .= LBound
    if isnan(p.NumAssets) == true
        p.NumAssets = NAssets
    end
    return p
end


#=
setBoundsType set the type of boud. Possible values are :simple or :conditional. 

    setBoundsType(p, BoundT::Vector{Symbol})

=#
function setBoundsType(p, BoundT::Vector{Symbol})
    p.BoundType = BoundT
    if isnan(p.NumAssets) == true
        p.NumAssets = length(BoundT)
    end
end
function setBoundsType(p, BoundT::Symbol, NAssets = p.NumAssets)
    p.BoundType = Array{Symbol}(undef, NAssets)
    p.BoundType .= BoundT
end


#=
setBudget sets up budget constraints

    p = setBudget(p,LowerBudget)
    p = setBudget(p,LowerBudget,UpperBudget)

Note: Given bounds for a budget constraint in either LowerBudget or UpperBudget, budget constraints require any portfolio in Port to satisfy:
    LowerBudget <= sum(Port) <= UpperBudget
One or both constraints may be specified. The usual budget constraint for a fully invested portfolio is to have LowerBudget = UpperBudget = 1. However, if the portfolio has allocations in cash, the budget constraints can be used to specify the cash constraints. For example, if the portfolio can hold between 0% and 10% in cash, the budget constraint would be set up with
    p = setBudget(p, 0.9, 1)
=#
function setBudget(p, LBudget)
    p.LowerBudget = LBudget
    return p
end
function setBudget(p, LBudget, UBudget)
    p.LowerBudget = LBudget
    p.UpperBudget = UBudget
    return p
end


#=
setDefaultConstraints sets up portfolio constraints with nonnegative weights that sum to 1

    p = setDefaultConstraints(p)    
    p = setDefaultConstraints(p,NumAssets)

=#
function setDefaultConstraints(p) 
    p.LowerBound = zeros(p.NumAssets)
    p.LowerBudget = 1.0
    p.UpperBudget = 1.0
    p.BoundType = Array{Symbol}(undef, p.NumAssets)
    p.BoundType .= :simple
    return p
end
function setDefaultConstraints(p, NumAssets)
    p.NumAssets = NumAssets 
    setDefaultConstraints(p)
    return p
end


#=
setTurnover sets up maximum portfolio turnover constraint

    p = setTurnover(p,Turnover)
    p = setTurnover(p,Turnover,InitPort)

Given an upper bound for portfolio turnover in Turnover and an initial portfolio in InitPort, the turnover constraint requires any portfolio in Port to satisfy the following:
    1' *1/2* | Port - InitPort | <= Turnover
=#
function setTurnover(p, turnover)
    p.Turnover = turnover
    if p.InitPort == Float64[] && isnan(p.NumAssets) == false 
        p.InitPort = zeros(p.NumAssets)
    end
    return p
end
function setTurnover(p, turnover, inip)
    p.Turnover = turnover
    p.InitPort = inip
    return p
end


#=
setOneWayTurnover sets up one-way portfolio turnover constraints

    p = setOneWayTurnover(p,BuyTurnover)
    p = setOneWayTurnover(p,BuyTurnover,SellTurnover,InitPort)

Given an initial portfolio in InitPort and an upper bound for portfolio turnover on purchases in BuyTurnover or sales in SellTurnover, the one-way turnover constraints require any portfolio Port to satisfy the following:
    1' * max{0, Port - InitPort} <= BuyTurnover
    1' * max{0, InitPort - Port} <= SellTurnover
One-way turnover constraints depend upon an initial or current portfolio, which is assumed to be zero if not set when the turnover constraints are set.
=#
function setOneWayTurnover(p, bturnover::Float64)
    p.BuyTurnover = bturnover
    if p.InitPort == Float64[] && isnan(p.NumAssets) == false 
        p.InitPort = zeros(p.NumAssets)
    end
    return p
end
function setOneWayTurnover(p, bturnover::Float64, sturnover::Float64)
    p.BuyTurnover = bturnover
    p.SellTurnover = sturnover
    if p.InitPort == Float64[] && isnan(p.NumAssets) == false 
        p.InitPort = zeros(p.NumAssets)
    end
    return p
end
function setOneWayTurnover(p, bturnover::Float64, sturnover::Float64, inip::Vector{Float64})
    p.BuyTurnover = bturnover
    p.SellTurnover = sturnover
    p.InitPort = inip
    if isnan(p.NumAssets) == true
        p.NumAssets = length(inip)
    end
    return p
end


#=
setTrackingPort sets up benchmark portfolio for tracking error constraint. 
    
    p = setTrackingPort(p,TracPort)

Input arguments:
    TracPort - Tracking portfolio weights, specified using a vector. 
=#
function setTrackingPort(p, trackport0::Vector{Float64})
    p.TrackingPort = trackport0
    if isnan(p.NumAssets) == true
        p.NumAssets = length(trackport0)
    end
    return p
end

#=
setTrackingError sets up maximum portfolio tracking error constraint

    p = setTrackingError(p, TrackingError)
    p = setTrackingError(p, TrackingError, TrackingPort)

Given an upper bound for portfolio tracking error in TrackingError and a tracking portfolio in TrackingPort, the tracking error constraint requires any portfolio in Port to satisfy
    (Port - TrackingPort)' * AssetCovar * (Port - TrackingPort) <= TrackingError^2 

Note: The tracking error constraints can be used with any of the other supported constraints in the Portfolio object without restrictions. However, since the portfolio set necessarily and sufficiently must be a non-empty compact set, the application of a tracking error constraint can result in an empty portfolio set. Use estimateBounds to confirm that the portfolio set is non-empty and compact.
=#
function setTrackingError(p, trackerror)
    p.TrackingError = trackerror
    return p
end
function setTrackingError(p, trackerror, trackport0)
    p.TrackingError = trackerror
    p.TrackingPort = trackport0
    if isnan(p.NumAssets) == true
        p.NumAssets = length(trackport0)
    end
    return p
end


#=
setMinMaxNumAssets sets cardinality constraints on the number of assets invested

    p = setMinMaxNumAssets(p, MinNumAssets, MaxNumAssets)

When working with a Portfolio object, the setMinMaxNumAssets function enables you to set up limits on the number of assets invested. These limits are also known as cardinality constraints. When managing a portfolio, it is common that you want to invest in at least a certain number of assets. In addition, you should also clearly define the weight requirement for each invested asset. You can do this using setBounds with a 'Conditional' BoundType. If you do not specify a 'Conditional' BoundType, the optimizer cannot understand which assets are invested assets and cannot formulate the MinNumAssets constraint.
=#
function setMinMaxNumAssets(p, minassets, maxassets)
    p.MinNumAssets = minassets
    p.MaxNumAssets = maxassets
    return p
end


#=
MINLP returns true if BoundType for some asset is set to :conditional , MinNumAssets > 0, or MaxNumAssets > 0, and false otherwise

    MINLP(p)

=#
function MINLP(p)
    if p.MinNumAssets > 0
        return true
    elseif p.MaxNumAssets > 0
        return true
    elseif p.BoundType != Symbol("")
        for i in 1:length(p.BoundType)
            if p.BoundType[i] == :conditional
                return true
            end
        end
        return false
    else
        return false
    end
end


#=
setEquality sets up linear equality constraints for portfolio weights

    p = setEquality(p, AEqual, bEqual)

Given linear equality constraint matrix AEquality and vector bEquality, every weight in a portfolio Port must satisfy the following:
    AEqual * Port = bEqual
=#
function setEquality(p, AEqual, bEqual)
    p.AEquality = AEqual
    p.bEquality = bEqual
    if isnan(p.NumAssets) == true
        p.NumAssets = size(AEqual,2)
    end
    return p
end


#=
setQuadEquality sets up quadratic equality constraints for portfolio weights

    p = setEquality(p, AQuadEqual, bQuadEqual)

Given quadratic equality constraint matrix AQuadEqual and vector bQuadEqual, every weight in a portfolio Port must satisfy the following:
    Port' * AEqual * Port = bQuadEqual
=#
function setQuadEquality(p, AQuadEqual, bQuadEqual)
    p.AQuadEquality = AQuadEqual
    p.bQuadEquality = bQuadEqual
    if isnan(p.NumAssets) == true
        p.NumAssets = size(AQuadEqual,2)
    end
    return p
end


#=
setInequality sets up linear inequality constraints for portfolio weights

    p = setInequality(p, AInequal, bInequal)

Given a linear inequality constraint matrix AInequality and vector bInequality, every weight in a portfolio Port must satisfy the following:
    AInequality * Port <= bInequality
=#
function setInequality(p, AInequal, bInequal)
    p.AInequality = AInequal
    p.bInequality = bInequal
    if isnan(p.NumAssets) == true
        p.NumAssets = size(AInequal,2)
    end
    return p
end


#=
setQuadInequality sets up quadratic inequality constraints for portfolio weights

    p = setQuadInequality(p, AQuadInequal, bQuadInequal)

Given a quadratic inequality constraint matrix AQuadInequality and vector bQuadInequality, every weight in a portfolio Port must satisfy the following:
    Port' * AInequality * Port <= bInequality
=#
function setQuadInequality(p, AQuadInequal, bQuadInequal)
    p.AQuadInequality = AQuadInequal
    p.bQuadInequality = bQuadInequal
    if isnan(p.NumAssets) == true
        p.NumAssets = size(AQuadInequal,2)
    end
    return p
end


#=
setGroups sets up group constraints for portfolio weights

    p = setGroups(p, GroupMatrix, LGroup)
    p = setGroups(p, GroupMatrix, LGroup, UGroup)

The group matrix GroupMatrix is usually an indicator of membership in groups, which means that its elements are usually either 'true' or 'false'. Given GroupMatrix and either LowerGroup or UpperGroup, a portfolio Port must satisfy the following:
    LowerGroup <= GroupMatrix * Port <= UpperGroup

Example: Suppose you have a portfolio of five assets and you want to ensure that the first three assets constitute at most 30% of your portfolio. Given a Portfolio object p, set the group constraints with the following
    G = [ true true true false false ]
    p = Portfolio()
    p = setGroups(p, G, [], 0.3)
=#
function setGroups(p, GMatrix, LGroup)
    p.GroupMatrix = GMatrix
    p.LowerGroup = LGroup
    if isnan(p.NumAssets) == true
        p.NumAssets = size(p.GroupMatrix,2)
    end
    return p
end
function setGroups(p, GMatrix, LGroup, UGroup)
    p.GroupMatrix = GMatrix
    p.LowerGroup = LGroup
    p.UpperGroup = UGroup
    if isnan(p.NumAssets) == true
        p.NumAssets = size(p.GroupMatrix,2)
    end
    return p
end


#=
setGroupRatio sets up group ratio constraints for portfolio weights

    p = setGroupRatio(p, GroupA, GroupB, LowerRatio)
    p = setGroupRatio(p, GroupA, GroupB, LowerRatio, UpperRatio)

Given base and comparison group matrices GroupA and GroupB and LowerRatio or UpperRatio bounds, group ratio constraints require any portfolio in Port to satisfy the following:
    (GroupB * Port) .* LowerRatio <= GroupA * Port <= (GroupB * Port) .* UpperRatio

Example: Suppose you want to ensure that the ratio of financial to nonfinancial companies in your portfolio never exceeds 50%. Assume you have six assets with three financial companies (assets 1-3) and three nonfinancial companies (assets 4-6). Group ratio constraints can be set with:
    GA = [true true true false false false]    % financial companies
    GB = [false false false true true true]    % nonfinancial companies
    p = Portfolio()
    p = setGroupRatio(p, GA, GB, [], 0.5)
=#
function setGroupRatio(p, GA, GB, LRatio)
    p.GroupA = GA
    p.GroupB = GB
    p.LowerRatio = LRatio
    if isnan(p.NumAssets) == true
        p.NumAssets = size(GA,2)
    end
    return p
end 
function setGroupRatio(p, GA, GB, LRatio, URatio)
    p.GroupA = GA
    p.GroupB = GB
    p.LowerRatio = LRatio
    p.UpperRatio = URatio
    if isnan(p.NumAssets) == true
        p.NumAssets = size(GA,2)
    end
    return p
end 


#=
addEquality adds linear equality constraints for portfolio weights to existing constraints

    p = addEquality(p, AEqual, bEqual)

Given a linear equality constraint matrix AEquality and vector bEquality, every weight in a portfolio Port must satisfy the following:
    AEquality * Port = bEquality
This function "stacks" additional linear equality constraints onto any existing linear equality constraints that exist in the input portfolio object. If no constraints exist, this method is the same as setEquality.
=#
function addEquality(p, AEqual, bEqual)
    if p.AEquality == Array{Float64}(undef, 0, 0) || p.bEquality == Float64[]
        p.AEquality = AEqual
        p.bEquality = bEqual
    else
        p.AEquality = cat(p.AEquality, AEqual, dims = 1)
        p.bEquality = cat(p.bEquality, bEqual, dims = 1)
    end
    if isnan(p.NumAssets) == true
        p.NumAssets = size(AEqual,2)
    end
    return p
end


#=
addQuadEquality adds quadratic equality constraints for portfolio weights to existing constraints

    p = addQuadEquality(p, AQuadEqual, bQuadEqual)

Given a quadratic equality constraint matrix AQuadEquality and vector bQuadEquality, every weight in a portfolio Port must satisfy the following:
    Port' * AEquality * Port = bEquality
This function "stacks" additional quadratic equality constraints onto any existing quadratic equality constraints that exist in the input portfolio object. If no constraints exist, this method is the same as setQuadEquality.
=#
function addQuadEquality(p, AQuadEqual, bQuadEqual)
    if p.AQuadEquality == Array{Float64}(undef, 0, 0, 0) || p.bQuadEquality == Float64[]
        p.AQuadEquality = AQuadEqual
        p.bQuadEquality = bQuadEqual
    else
        p.AQuadEquality = cat(p.AQuadEquality, AQuadEqual, dims = 3)    # Note dims = 3
        p.bQuadEquality = cat(p.bQuadEquality, bQuadEqual, dims = 1)
    end
    if isnan(p.NumAssets) == true
        p.NumAssets = size(AQuadEqual,2)
    end
    return p
end


#=
addInequality adds linear inequality constraints for portfolio weights to existing constraint

    p = addInequality(p, AInequal, bInequal)

Given a linear inequality constraint matrix AInequality and vector bInequality, every weight in a portfolio Port must satisfy the following:
    AInequality * Port <= bInequality
This function "stacks" additional linear inequality constraints onto any existing linear inequality constraints that exist in the input portfolio object. If no constraints exist, this function is the same as setInequality.
=#
function addInequality(p, AInequal, bInequal)
    if p.AInequality == Array{Float64}(undef, 0, 0) || p.bInequality == Float64[]
        p.AInequality = AInequal
        p.bInequality = bInequal
    else
        p.AInequality = cat(p.AInequality, AInequal, dims = 1)
        p.bInequality = cat(p.bInequality, bInequal, dims = 1)
    end
    if isnan(p.NumAssets) == true
        p.NumAssets = size(AInequal,2)
    end
    return p
end


#=
addQuadInequality adds quadratic inequality constraints for portfolio weights to existing constraint

    p = addQuadInequality(p, AQuadInequal, bQuadInequal)

Given a quadratic inequality constraint matrix AQuadInequality and vector bQuadInequality, every weight in a portfolio Port must satisfy the following:
    Port' * AInequality * Port <= bInequality
This function "stacks" additional quadratic inequality constraints onto any existing quadratic inequality constraints that exist in the input portfolio object. If no constraints exist, this function is the same as setQuadInequality.
=#
function addQuadInequality(p, AQuadInequal, bQuadInequal)
    if p.AQuadInequality == Array{Float64}(undef, 0, 0, 0) || p.bQuadInequality == Float64[]
        p.AQuadInequality = AQuadInequal
        p.bQuadInequality = bQuadInequal
    else
        p.AQuadInequality = cat(p.AQuadInequality, AQuadInequal, dims = 3)  # Note dims = 3
        p.bQuadInequality = cat(p.bQuadInequality, bQuadInequal, dims = 1)
    end
    if isnan(p.NumAssets) == true
        p.NumAssets = size(AQuadInequal,2)
    end
    return p
end


#=
addGroups adds group constraints for portfolio weights to existing group constraints

    p = addGroups(p, GMatrix, LGroup)
    p = addGroups(p, GMatrix, LGroup, UGroup)

Given GroupMatrix and either LowerGroup or UpperGroup, a portfolio Port must satisfy the following:
    LowerGroup <= GroupMatrix * Port <= UpperGroup
=#
function addGroups(p, GMatrix, LGroup)
    if p.GroupMatrix == Array{Bool}(undef, 0, 0) || p.LowerGroup == Float64[]
        p.GroupMatrix = GMatrix
        p.LowerGroup = LGroup
    else
        p.GroupMatrix = cat(p.GroupMatrix, GMatrix, dims = 1)
        p.LowerGroup = cat(p.LowerGroup, LGroup, dims = 1)
    end
    if isnan(p.NumAssets) == true
        p.NumAssets = size(GMatrix,2)
    end
    return p
end
function addGroups(p, GMatrix, LGroup, UGroup)
    if p.GroupMatrix == Array{Bool}(undef, 0, 0) || p.LowerGroup == Float64[] || p.UpperGroup == Float64[] 
        p.GroupMatrix = GMatrix
        p.LowerGroup = LGroup
        p.UpperGroup = UGroup
    else
        p.GroupMatrix = cat(p.GroupMatrix, GMatrix, dims = 1)
        p.LowerGroup = cat(p.LowerGroup, LGroup, dims = 1)
        p.UpperGroup = cat(p.UpperGroup, UGroup, dims = 1)
    end
    if isnan(p.NumAssets) == true
        p.NumAssets = size(GMatrix,2)
    end
    return p
end


#=
addGroupRatio adds group ratio constraints for portfolio weights to existing group ratio constraints

    p = addGroupRatio(p, GroupA, GroupB, LowerRatio)
    p = addGroupRatio(p, GroupA, GroupB, LowerRatio, UpperRatio)
=#
function addGroupRatio(p, GA, GB, LRatio)
    if p.GroupA == Array{Bool}(undef, 0, 0) || p.GroupB == Array{Bool}(undef, 0, 0) || p.LowerRaio == Float64[]
        p.GroupA = GA
        p.GroupB = GB
        p.LowerRatio = LRatio
    else
        p.GroupA = cat(p.GroupA, GA, dims = 1)
        p.GroupB = cat(p.GroupB, GB, dims = 1)
        p.LowerRatio = cat(p.LowerRatio, LRatio, dims = 1)
    end
    if isnan(p.NumAssets) == true
        p.NumAssets = size(GA,2)
    end
    return p
end   
function addGroupRatio(p, GA, GB, LRatio, URatio)
    if p.GroupA == Array{Bool}(undef, 0, 0) || p.GroupB == Array{Bool}(undef, 0, 0) || p.LowerRatio == Float64[] || p.UpperRatio == Float64[] 
        p.GroupA = GA
        p.GroupB = GB
        p.LowerRatio = LRatio
        p.UpperRatio = URatio
    else
        p.GroupA = cat(p.GroupA, GA, dims = 1)
        p.GroupB = cat(p.GroupB, GB, dims = 1)
        p.LowerRatio = cat(p.LowerRatio, LRatio, dims = 1)
        p.UpperRatio = cat(p.UpperRatio, URatio, dims = 1)
    end
    if isnan(p.NumAssets) == true
        p.NumAssets = size(GA,2)
    end
    return p
end    

#=
getAssetList obtain the list of assets

    AssetList = getAssetList(p)

=#
function getAssetList(p)
    return p.AssetList
end


#=
getAssetMoments obtains mean and covariance of asset returns from Portfolio object

    (AssetMean, AssetCovar) = getAssetMoments(p)

=#
function getAssetMoments(p)
    return p.AssetMean, p.AssetCovar
end


#=
getBounds obtains bounds for portfolio weights 

    (LowerBound, UpperBound) = getBounds(p)

=#
function getBounds(p)
    return p.LowerBound, p.UpperBound
end


#=
getBudget obtains budget constraint bounds 

    (LowerBudget, UpperBudget) = getBudget(p)

=#
function getBudget(p)
    return p.LowerBudget, p.UpperBudget
end


#=
getCosts obtains buy and sell transaction costs 

    (BuyCost, SellCost) = getCosts(p)

=#
function getCosts(p)
    return p.BuyCost, p.SellCost
end
    
    
#=
getEquality obtains equality constraint arrays 

    (AEquality, bEquality) = getEquality(p)

=#
function getEquality(p)
    return p.AEquality, p.bEquality
end


#=
getQuadEquality obtains quadratic equality constraint arrays 

    (AQuadEquality, bQuadEquality) = getQuadEquality(p)

=#
function getQuadEquality(p)
    return p.AQuadEquality, p.bQuadEquality
end


#=
getGroupRatio obtains group ratio constraint arrays 

    (GroupA, GroupB, LowerRatio, UpperRatio) = getGroupRatio(p)

=#
function getGroupRatio(p)
    return p.GroupA, p.GroupB, p.LowerRatio, p.UpperRatio
end


#=
getGroups obtains group constraint arrays 

    (GroupMatrix, LowerGroup, UpperGroup) = getGroups(p)

=#
function getGroups(p)
    return p.GroupMatrix, p.LowerGroup, p.UpperGroup
end


#= 
getInequality obtains inequality constraint arrays 

    (AInequality,bInequality) = getInequality(p)

=#
function getInequality(p)
    return p.AInequality, p.bInequality
end


#= 
getQuadInequality obtains inequality constraint arrays 

    (AQuadInequality,bQuadInequality) = getQuadInequality(p)

=#
function getQuadInequality(p)
    return p.AQuadInequality, p.bQuadInequality
end


#=
getOneWayTurnover obtains one-way turnover constraints

    (BuyTurnover,SellTurnover) = getOneWayTurnover(p)

=#
function getOneWayTurnover(p)
    return p.BuyTurnover, p.SellTurnover
end


#=
checkFeasibility checks feasibility of input portfolios against input portfolio object.

    status = checkFeasibility(p, x, toler = 1.0e-8)

Input arguments:
    x       — Portfolio to check 
    toler   - feasibility tolerance from AbstractPortfolio
=#
function checkFeasibility(p, x, toler = 1.0e-8) 
    n = length(x)        # number assets
    status = true
    if p.AEquality != Array{Float64}(undef, 0, 0) && p.bEquality != Float64[]
        for k in 1:length(p.bEquality)
            if abs(p.AEquality[k,:]' * x - p.bEquality[k]) > toler
                status = false
                @warn "Portfolio is not feasible: check Linear Equality constraint $(k)"
            end
        end
    end
    if p.AQuadEquality != Array{Float64}(undef, 0, 0, 0) && p.bQuadEquality != Float64[]
        for k in 1:length(p.bQuadEquality)
            if abs(x' * p.AQuadEquality[:, :, k] * x - p.bQuadEquality[k]) > toler
                status = false
                @warn "Portfolio is not feasible: check Quadratic Equality constraint $(k)"
            end
        end
    end
    if p.AInequality != Array{Float64}(undef, 0, 0) && p.bInequality != Float64[]
        for k in 1:length(p.bInequality)
            if p.AInequality[k,:]' * x - p.bInequality[k] > toler
                status = false
                @warn "Portfolio is not feasible: check Linear Inequality constraint $(k)"
            end
        end
    end
    if p.AQuadInequality != Array{Float64}(undef, 0, 0, 0) && p.bQuadInequality != Float64[]
        for k in 1:length(p.bQuadInequality)
            if x' * p.AQuadInequality[:, :, k] * x - p.bQuadInequality[k] > toler
                status = false
                @warn "Portfolio is not feasible: check Quadratic Inequality constraint $(k)"
            end
        end
    end
    if p.LowerBound != Float64[]
        k = 1
        while p.LowerBound[k] - x[k] <= toler
            k += 1
            if k > n break
            end
        end
        if k <= n
            status = false
            @warn "Portfolio is not feasible: check portfolio weigth[$k] against LowerBond constraint"
        end
    end
    if p.UpperBound != Float64[]
        k = 1
        while x[k] - p.UpperBound[k] <= toler
            k += 1
            if k > n break
            end
        end
        if k <= n
            status = false
            @warn "Portfolio is not feasible: check portfolio weigth[$k] against UpperBond constraint"
        end
    end
    if isnan(p.LowerBudget) == false
        if p.LowerBudget - sum(x) > toler
            status = false
            @warn "Portfolio is not feasible: check Lower Budget constraint"
        end
    end
    if isnan(p.UpperBudget) == false
        if sum(x) - p.UpperBudget > toler
            status = false
            @warn "Portfolio is not feasible: check Upper Budget constraints"
        end
    end
    if p.GroupMatrix != Matrix{Bool}(undef, 0, 0) && p.LowerGroup != Float64[]
        Z = p.LowerGroup .- p.GroupMatrix * x
        for i = 1:length(p.LowerGroup)
            if Z[i] > toler
                status = false
                @warn "Portfolio is not feasible: check portfolio Group constraint $i against Lower Bound $(p.LowerGroup[i])"
            end
        end
    end
    if p.GroupMatrix != Matrix{Bool}(undef, 0, 0) && p.UpperGroup != Float64[]
        Z = p.GroupMatrix * x .- p.UpperGroup
        for i = 1:length(p.UpperGroup)
            if Z[i] > toler
                status = false
                @warn "Portfolio is not feasible: check portfolio Group constraint $i against Upper Bound $(p.UpperGroup[i])"
            end
        end
    end
    if p.GroupA != Array{Bool}(undef, 0, 0) && p.GroupB != Array{Bool}(undef, 0, 0) && p.LowerRatio != Float64[]
        Z = (p.GroupB .* p.LowerRatio - p.GroupA) * x  
        for i = 1:length(p.LowerRatio)                                  
            if Z[i] > toler
                status = false
                @warn "Portfolio is not feasible: check portfolio Group Ratio constraint $i against Lower Ratio $(p.LowerRatio[i])"
            end
        end
    end
    if p.GroupA != Array{Bool}(undef, 0, 0) && p.GroupB != Array{Bool}(undef, 0, 0) && p.UpperRatio != Float64[]
        Z = (p.GroupA -  p.GroupB .* p.UpperRatio) * x  
        for i = 1:length(p.UpperRatio)
            if Z[i] > toler
                status = false
                @warn "Portfolio is not feasible: check portfolio Group Ratio constraint $i against Upper Ratio $(p.UpperRatio[i])"
            end
        end
    end
    if isnan(p.Turnover) == false
        if p.InitPort == Float64[]
            x0 = zeros(n)
        else
            x0 = p.InitPort
        end
        if 0.5 * sum(abs.(x - x0)) - p.Turnover > toler
            status = false
            @warn "Portfolio is not feasible: check portfolio Turnover constraint"
        end
    end
    if isnan(p.BuyTurnover) == false
        if p.InitPort == Float64
            x0 = zeros(n)
        else
            x0 = p.InitPort
        end
        if sum(max.(0, (x - x0))) - p.BuyTurnover > toler
            status = false
            @warn "Portfolio is not feasible: check portfolio Buy Turnover constraint"
        end
    end
    if isnan(p.SellTurnover) == false
        if p.InitPort == Float64
            x0 = zeros(n)
        else
            x0 = p.InitPort
        end
        if sum(max.(0, (x0 - x))) - p.SellTurnover > toler
            status = false
            @warn "Portfolio is not feasible: check portfolio Sell Turnover constraint"
        end
    end
    if p.TrackingPort != Float64[] && isnan(p.TrackingError) == false
        if (x - p.TrackingPort)' * p.AssetCovar * (x - p.TrackingPort) - p.TrackingError^2 > toler
            status = false
            @warn "Portfolio is not feasible: check Tracking Portfolio constraint"
        end
    end
    return status
end


function set_model(p, optimizer)
    # Initialising optimising model
    model = Model()

    # Choosing the optimiser and its setting options
    # Setting environment variables for FICO Xpress and Gurobi in ~/.zshenv
    #       export XPRESSDIR=/Applications/FICO\ Xpress/xpressmp
    #       export GUROBI_HOME=/Library/gurobi911/mac64
    if MINLP(p) == true # Choose a mixed integer nonlinear programming (MINLP) solver
        if optimizer == :CPLEX
            set_optimizer(model, CPLEX.Optimizer)
        elseif optimizer == :Gurobi
            set_optimizer(model, Gurobi.Optimizer)
        elseif optimizer == :Juniper
            nl_solver= optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
            mip_solver = optimizer_with_attributes(Cbc.Optimizer, "logLevel" => 0)
            set_optimizer(model, optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>nl_solver, "mip_solver"=>mip_solver))
        elseif optimizer == :Mosek
            set_optimizer(model, Mosek.Optimizer)
        elseif optimizer == :SCIP
            set_optimizer(model, SCIP.Optimizer)
        elseif optimizer == :Xpress
            set_optimizer(model, Xpress.Optimizer)
        else
            error("$optimizer is not a valid optimiser choice. Possible choices are: :CPLEX, :Gurobi, :Juniper, :Mosek, :SCIP and :Xpress")
        end    
    else
        if optimizer == :COSMO
            set_optimizer(model, COSMO.Optimizer)
        elseif optimizer == :CPLEX
            set_optimizer(model, CPLEX.Optimizer)
        elseif optimizer == :Gurobi
            set_optimizer(model, Gurobi.Optimizer)
        elseif optimizer == :Ipopt
            set_optimizer(model, Ipopt.Optimizer)
        elseif optimizer == :Juniper
            nl_solver= optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0)
            mip_solver = optimizer_with_attributes(Cbc.Optimizer, "logLevel" => 0)
            set_optimizer(model, optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>nl_solver, "mip_solver"=>mip_solver))
        elseif optimizer == :Mosek
            set_optimizer(model, Mosek.Optimizer)
        elseif optimizer == :SCIP
            set_optimizer(model, SCIP.Optimizer)
        elseif optimizer == :SCS
            set_optimizer(model, SCS.Optimizer)
        elseif optimizer == :Xpress
            set_optimizer(model, Xpress.Optimizer)
        else
            error("$optimizer is not a valid optimiser choice. Possible choices are: :COSMO, :CPLEX, :Gurobi, :Ipopt, :Juniper, :Mosek, :SCS, :SCIP and :Xpress")
        end
    end
    # Setting solver options. The string names of the attributes are specific to each solver. One should consult the solver's documentation to find the attributes of interest. 
    set_silent(model)
    # set_optimizer_attribute(model, "tol", 1e-6)  # ERROR: LoadError: type Settings has no field tol
    # set_optimizer_attribute(model, "max_iter", 100)

    # Variables w[1:n]: n = p.NumAssets portfolio weigths
    n = p.NumAssets
    @variable(model, w[1:n])

    # Variables ν[1:n] only take binary values 0 and 1 to indicate whether the corresponding asset is invested (1) or not invested (0).
    if MINLP(p) == true
        @variable(model, v[1:n], binary = true)
    end

    # Setting constraints
    # Cardinality constraints and bound constraints with potentially semicontinuous variables
    if MINLP(p) == true
        # Cardinality constraints: minimum and maximum number of assets
        if isnan(p.MinNumAssets) == false
            @constraint(model, sum(v) >= p.MinNumAssets)
        end
        if isnan(p.MaxNumAssets) == false
            @constraint(model, sum(v) <= p.MaxNumAssets)
        end
        # Bound constraints with semicontinuous variables: lower and upper bounds
        if p.LowerBound != Float64[]
            for i in 1:n
                if p.BoundType[i] == :simple && p.LowerBound[i] > 0.0
                    @constraint(model, v[i] == 1)
                end
                @constraint(model, p.LowerBound[i] * v[i] <= w[i])
            end
        else
            for i in 1:n
                @constraint(model, floatmax() * v[i] <= w[i])
            end
        end
        if p.UpperBound != Float64[]
            for i in 1:n
                @constraint(model, w[i] <= p.UpperBound[i] * v[i])
            end
        else
            for i in 1:n
                @constraint(model, w[i] <= floatmax()* v[i])
            end
        end
    else # No cardinality constraints and no bound constraints with potentially semicontinuous variables
        if p.LowerBound != Float64[]
            for i = 1:n
                @constraint(model, w[i] .>= p.LowerBound[i]) 
            end
        end
        if p.UpperBound != Float64[]
            for i = 1:n
                @constraint(model, w[i] .<= p.UpperBound[i]) 
            end
        end
    end

    # Budget constraints
    if  p.LowerBudget == p.UpperBudget
        @constraint(model, budget, sum(w) .== p.LowerBudget)  # naming constraint budget, one will need to remove it for maximising the Sharpe ratio
    end
    if isnan(p.LowerBudget) == false && p.LowerBudget != p.UpperBudget
        @constraint(model, sum(w) .>= p.LowerBudget)
    end
    if isnan(p.UpperBudget) == false && p.LowerBudget != p.UpperBudget
        @constraint(model, sum(w) .<= p.UpperBudget)
    end

    # Turnover constraints
    if isnan(p.Turnover) == false
        #@constraint(model, 0.5 * sum(abs.(w - p.InitPort)) .<= p.Turnover)
        @variable(model, t[1:n])
        @constraint(model, t .>= w - p.InitPort)
        @constraint(model, t .>= p.InitPort - w)
        @constraint(model, 0.5 * sum(t) <= p.Turnover)
    end
    if isnan(p.BuyTurnover) == false
        #@constraint(model, sum(max.(0, (w - p.InitPort))) .<= p.BuyTurnover)
        @variable(model, t1[1:n])
        @constraint(model, t1 .>= w - p.InitPort)
        @constraint(model, t1 .>= 0)
        @constraint(model, sum(t1) <= p.BuyTurnover)
    end
    if isnan(p.SellTurnover) == false
        #@constraint(model, sum(max.(0, (p.InitPort - w))) .<= p.SellTurnover)
        @variable(model, u1[1:n])
        @constraint(model, u1 .>= p.InitPort - w)
        @constraint(model, u1 .>= 0)
        @constraint(model, sum(u1) <= p.SellTurnover)
    end

    # Net return target (net of buy and sell costs) inequality constraint 
    if isnan(p.NetReturn) == false 
        @variable(model, t2[1:n])
        @constraint(model, t2 .>= w - p.InitPort)
        @constraint(model, t2 .>= 0)
        @variable(model, u2[1:n])
        @constraint(model, u2 .>= p.InitPort - w)
        @constraint(model, u2 .>= 0)
        @constraint(model, sum(p.AssetMean[i] * w[i] - (p.BuyCost[i] * t2[i] + p.SellCost[i] * u2[i]) for i in 1:n) >= p.NetReturn)
    end

    # Tracking error constraint
    if p.TrackingPort != Float64[] && isnan(p.TrackingError) == false
        @constraint(model, (w - p.TrackingPort)' * p.AssetCovar * (w - p.TrackingPort) .<= p.TrackingError^2) 
    end

    # Linear equality constraints
    if p.AEquality != Array{Float64}(undef, 0, 0) && p.bEquality != Float64[]
        nconstraints = length(p.bEquality)
        if nconstraints == 1
            @constraint(model, p.AEquality' * w .== p.bEquality)
        else
            for i in 1:nconstraints
                @constraint(model, p.AEquality[i,:]' * w .== p.bEquality[i])
            end
        end
    end

    # Linear inequality constraints
    if p.AInequality != Array{Float64}(undef, 0, 0) && p.bInequality != Float64[]
        nconstraints = length(p.bEquality)
        if nconstraints == 1
            @constraint(model, p.AInequality' * w .== p.bEquality)
        else
            for i in 1:length(p.bInequality)
                @constraint(model, p.AInequality[i,:]' * w .<= p.bInequality[i])
            end
        end
    end

    # Quadratic equality constraints
    if p.AQuadEquality != Array{Float64}(undef, 0, 0, 0) && p.bQuadEquality != Float64[]
        for i in 1:length(p.bQuadEquality)
            @constraint(model, w' * p.AQuadEquality[:, :, i] * w .== p.bQuadEquality[i])
        end
    end

    # Quadratic inequality constraints
    if p.AQuadInequality != Array{Float64}(undef, 0, 0, 0) && p.bQuadInequality != Float64[]
        for i in 1:length(p.bQuadInequality)
            @constraint(model, w' * p.AQuadInequality[:, :, i] * w .<= p.bQuadInequality[i])
        end
    end

    # Constraints by group of assets
    if p.GroupMatrix != Matrix{Bool}(undef, 0, 0) && p.LowerGroup != Float64[]
        ZL = p.LowerGroup .- p.GroupMatrix * w
        for i = 1:length(p.LowerGroup)
            @constraint(model, ZL[i] .<= 0)
        end
    end
    if p.GroupMatrix != Matrix{Bool}(undef, 0, 0) && p.UpperGroup != Float64[]
        ZU = p.GroupMatrix * w .- p.UpperGroup
        for i = 1:length(p.UpperGroup)
            @constraint(model, ZU[i] .<= 0)
        end
    end
    if p.GroupA != Array{Bool}(undef, 0, 0) && p.GroupB != Array{Bool}(undef, 0, 0) && p.LowerRatio != Float64[]
        ZLR = (p.GroupB .* p.LowerRatio - p.GroupA) * w
        for i = 1:length(p.LowerRatio)                                  
            ZLR[i] .<= 0
        end
    end
    if p.GroupA != Array{Bool}(undef, 0, 0) && p.GroupB != Array{Bool}(undef, 0, 0) && p.UpperRatio != Float64[]
        ZUR = (p.GroupA -  p.GroupB .* p.UpperRatio) * w  
        for i = 1:length(p.UpperRatio)
            ZUR[i] .<= 0
        end
    end
    
    if MINLP(p) == true
        return model, w, v
    else
        return model, w 
    end
end


function optimise_model(model, w, objtype, f)
# Setting the objective function. Use the @objective macro to set linear and quadratic objective functions in a JuMP model. Depending on the problem, there is also @NLObjective.
    if objtype == :min
        @objective(model, Min, f(w))
    elseif objtype == :max
        @objective(model, Max, f(w))
    end
    # Getting optimal solution 
    optimize!(model)
    return value.(w), objective_value(model)
end

function optimise_model(model, w, v, objtype, f)
    # Setting the objective function. Use the @objective macro to set linear and quadratic objective functions in a JuMP model. Depending on the problem, there is also @NLObjective.
        if objtype == :min
            @objective(model, Min, f(w,v))
        elseif objtype == :max
            @objective(model, Max, f(w,v))
        end 
        # Getting optimal solution 
        optimize!(model)
        return value.(w), objective_value(model), value.(v)
end


#=
optimalp estimates optimal portfolios with (i) targeted portfolio risks or returns, or (ii) a specified number of optimal portfolios on the efficient frontier, including endpoints

    optimalp(p, f, optype, optimizer = COSMO)

Input arguments:
    f           - Function to optimise 
    optype      - Type of optimisation. The options are:
                    :max
                    :min
    optimizer   - The optimizer to be used. Possible choices are:
                    COSMO       (LP, QP, SOCP, SDP)    / Apache
                    SCS         (LP, SOCP, SDP)        / MIT
                    Ipopt       (LP, QP, NLP)          / EPL
                    Juniper     (MISOCP, MINLP)        / MIT                
                    SCIP        (MILP, MINLP)          / ZIB
Examples:
    Minimum variance portfolio
        f(z) = z' * p.AssetCovar * z
        optimalp(p, f, :min)
    Maximum return portfolio
        f(z) = sum((p.AssetMean[i] * z[i]) for i = 1:length(p.AssetMean))
        optimalp(p, f, :max)
=#    
function optimalp(p, objtype, f, optimizer)
    if MINLP(p) == true
        (model, w, v) = set_model(p, optimizer)
        optimise_model(model, w, v, objtype, f)
    else
        (model, w) = set_model(p, optimizer)
        optimise_model(model, w, objtype, f)
    end
end
function optimalp(p, objtype, f)
    if MINLP(p) == true
        (model, w, v) = set_model(p, :SCIP)
        optimise_model(model, w, v, objtype, f)
    else
        (model, w) = set_model(p, :Ipopt)
        optimise_model(model, w, objtype, f)
    end
end


#=
estimateMinimumVariancePortfolio estimates the minimum variance portfolio given portfolio constraints)

    (pwgt, portfolio_std, portfolio_return) = estimateMinimumVariancePortfolio(p, optimizer = :COSMO)

=#
function estimateMinimumVariancePortfolio(p, optimiser = :Ipopt) 
    if MINLP(p) == true
        if optimiser != :SCIP && optimiser != :Mosek && optimiser != :Xpress && optimiser != :Gurobi
            optimiser = :SCIP
        end
        f1(z,v) = z' * p.AssetCovar * z
        mvp = optimalp(p, :min, f1, optimiser)
    else
        f2(z) = z' * p.AssetCovar * z
        mvp = optimalp(p, :min, f2, optimiser)
    end
    mvp_weigths  = mvp[1]
    # Compute portfolio net return considering possible buy and sell transaction costs
    if p.BuyCost == Float64[] || p.SellCost == Float64 || p.InitPort == Float64
        mvp_return  = dot(p.AssetMean, mvp_weigths)[1]
    else
        mvp_return  = dot(p.AssetMean, mvp_weigths)[1] - (p.BuyCost' * max.(0, mvp_weigths - p.InitPort) + p.SellCost' * max.(0, p.InitPort - mvp_weigths))
    end
    # Compute portfolio risk (std)
    mvp_std = sqrt(mvp[2])
    # Return weigths, std and net return
    return mvp_weigths, mvp_std, mvp_return
end


#=
estimateMaximumNetReturnPortfolio estimates the maximum net return portfolio given portfolio constraints)

    (pwgt, portfolio_std, portfolio_return) = estimateMaximumNetReturnPortfolio(p, optimizer = :COSMO)

=#
function estimateMaximumNetReturnPortfolio(p, optimiser = :Ipopt)
    if MINLP(p) == true
        if optimiser != :SCIP && optimiser != :Mosek && optimiser != :Xpress && optimiser != :Gurobi && optimiser != :Juniper
            optimiser = :Juniper
        end
        (model, w, v) = set_model(p, optimiser)
    else
        (model, w) = set_model(p, optimiser) 
    end
    # Define objective according to the existence or not of buy and sell transaction costs
    if p.BuyCost == Float64[] || p.SellCost == Float64 
        @objective(model, Max, p.AssetMean' * w)
    else
        n = p.NumAssets
        @variable(model, t2[1:n])
        @constraint(model, t2 .>= w - p.InitPort)
        @constraint(model, t2 .>= 0)
        @variable(model, u2[1:n])
        @constraint(model, u2 .>= p.InitPort - w)
        @constraint(model, u2 .>= 0)
        @objective(model, Max, p.AssetMean' * w - (p.BuyCost' * t2 + p.SellCost' * u2))
    end
    optimize!(model)
    mrp_weigth = value.(w)
    mrp_return = objective_value(model)
    mrp_std     = sqrt(mrp_weigth' * p.AssetCovar * mrp_weigth)
    # Return weigths, std and net return
    return mrp_weigth, mrp_std, mrp_return
end


#=
estimateFrontierLimits estimates optimal portfolios at endpoints of efficient frontierm(minimum variance portfolio and maximum return portfolio given portfolio constraints)

    (pwgt, portfolio_std, portfolio_return) = estimateFrontierLimits(p, optimizer = :COSMO)

=#
function estimateFrontierLimits(p, optimiser = :Ipopt)
    sol1 = estimateMinimumVariancePortfolio(p, optimiser)
    sol2 = estimateMaximumNetReturnPortfolio(p, optimiser)
    return [sol1[1] sol2[1]], [sol1[2] sol2[2]], [sol1[3] sol2[3]]
end 


#=
estimateFrontierByReturn estimates optimal portfolios with targeted portfolio returns

    (pwgt, portfolio_std, portfolio_return) = estimateFrontierByReturn(p, target_return, optimiser = :COSMO)
		
When any one, or any combination of the constraints from 'Conditional' BoundType, 'MinNumAssets', and 'MaxNumAssets' are active, the portfolio problem is formulated as mixed integer programming problem and the MINLP solver is used.

TargetReturn specifies target returns for portfolios on the efficient frontier. If any TargetReturn values are outside the range of returns for efficient portfolios, the TargetReturn is replaced with the minimum or maximum efficient portfolio return, depending upon whether the target return is below or above the range of efficient portfolio returns.
=#
function estimateFrontierByReturn(p, target_return, optimiser = :Ipopt)
    if MINLP(p) == true
        if optimiser != :SCIP && optimiser != :Mosek && optimiser != :Xpress && optimiser != :Gurobi && optimiser != :Juniper
            optimiser = :Juniper
        end
    end
    limits = estimateFrontierLimits(p, optimiser)
    if target_return <= limits[3][1]
        @warn("Target return value is < than the lower limit of the feasible range [$(limits[3][1]), $(limits[3][2])]. Returning optimal portfolio weigths at this lower limit")    
        return limits[1][:,1], limits[2][1], limits[3][1]
    elseif target_return >= limits[3][2]
        @warn("Target return value is > than the upper limit of the feasible range [$(limits[3][1]), $(limits[3][2])]. Returning optimal portfolio weigths at this upper limit")
        limits[1][:,2], limits[2][2], limits[3][2]
    else
        q = deepcopy(p)
        setNetReturnTarget(q, target_return)
        if MINLP(q) == true
            f1(z,v) = z' * p.AssetCovar * z
            sol = optimalp(q, :min, f1, optimiser)
        else
            f2(z) = z' * p.AssetCovar * z
            sol = optimalp(q, :min, f2, optimiser)
        end
        return sol[1], sqrt(sol[2]), target_return  # weights, std, return
    end
end

#=
estimateFrontierByRisk estimates optimal portfolios with targeted portfolio risks

    (pwgt, portfolio_std, portfolio_return) = estimateFrontierByRisk(p, target_std, optimiser = :Ipopt)

If any target_std values are outside the range of risks for efficient portfolios, the target risk is replaced with the minimum or maximum efficient portfolio risk, depending on whether the target risk is below or above the range of efficient portfolio risks.
=#
function estimateFrontierByRisk(p, target_std, optimiser = :SCS)
    # Choose :SCS because :Ipopt has a lower precision and :COSMO appears to be unreliable
    if MINLP(p) == true
        if optimiser != :SCIP && optimiser != :Mosek && optimiser != :Xpress && optimiser != :Gurobi && optimiser != :Juniper
            optimiser = :Juniper
        end
    end
    limits = estimateFrontierLimits(p, optimiser)
    if target_std <= limits[2][1]
        @warn("Target std value is <= than the lower limit of the feasible range [$(limits[2][1]), $(limits[2][2])]. Returning optimal portfolio weigths at this lower limit")
        return limits[1][:,1], limits[2][1], limits[3][1]
    elseif target_std >= limits[2][2]
        @warn("Target std value is >= than the upper limit of the feasible range [$(limits[2][1]), $(limits[2][2])]. Returning optimal portfolio weigths at this upper limit")
        return limits[1][:,2], limits[2][2], limits[3][2]
    else
        q = deepcopy(p)
        addQuadInequality(q, p.AssetCovar, target_std^2)
        if MINLP(q) == true
            (model, w, v) = set_model(q, optimiser)
        else
            (model, w) = set_model(q, optimiser)
        end
        if q.BuyCost == Float64[] || q.SellCost == Float64 
            @objective(model, Max, q.AssetMean' * w)
        else
            n = q.NumAssets
            @variable(model, t2[1:n])
            @constraint(model, t2 .>= w - q.InitPort)
            @constraint(model, t2 .>= 0)
            @variable(model, u2[1:n])
            @constraint(model, u2 .>= q.InitPort - w)
            @constraint(model, u2 .>= 0)
            @objective(model, Max, q.AssetMean' * w - (q.BuyCost' * t2 + q.SellCost' * u2))
        end
        optimize!(model)
        pweigth = value.(w)
        preturn = objective_value(model)
        mrp_std = sqrt(pweigth' * p.AssetCovar * pweigth)
        return pweigth, target_std, preturn
    end
end


#=
estimateFrontier estimates a specified number of optimal portfolios on the efficient frontier

    (pwgt, portfolio_std, portfolio_return) = estimateFrontier(p, NumPorts)

Output arguments:
    pwgt - optimal portfolios on efficient frontier with specified number of portfolios spaced equally from minimum to maximum portfolio return (NumAssets x NumPorts)
=#
function estimateFrontier(p, NumPorts = 10, optmiser = :Ipopt)
    if MINLP(p) == true
        if optimiser != :SCIP && optimiser != :Mosek && optimiser != :Xpress && optimiser != :Gurobi && optimiser != :Juniper
            optimiser = :Juniper
        end
    end
    if NumPorts == 1
        return estimateMinimumVariancePortfolio(p, optmiser)
    elseif NumPorts == 2
        return estimateFrontierLimits(p, optmiser)
    else
        flim = estimateFrontierLimits(p, optmiser)
        gap = (flim[3][2] - flim[3][1])/(NumPorts - 1)
        weigths = flim[1][:,1]
        stds    = [flim[2][1]]
        returns = [flim[3][1]]
        for h in 1:(NumPorts-2)
            rpor = flim[3][1] + h * gap
            fr = estimateFrontierByReturn(p, rpor, optmiser)
            weigths = [weigths fr[1]]
            stds = push!(stds, fr[2])
            returns = push!(returns, fr[3])
        end
        weigths = [weigths flim[1][:,2]]
        stds = push!(stds, flim[2][2])
        returns = push!(returns, flim[3][2])
    end
    return weigths, stds, returns
end


#=
plotFrontier plots efficient frontier

    [prsk,pret] = plotFrontier(p, NumPorts = 10)

If the portfolio object has a name in the Name property, the name is displayed as the title of the plot. Otherwise, the plot is labeled “Efficient Frontier.”
If the portfolio object has an initial portfolio in the InitPort property, the initial portfolio is plotted and labeled.
If portfolio risks and returns are inputs, make sure that risks come first in the calling sequence. In addition, if portfolio risks and returns are not sorted in ascending order, this method performs the sort. On output, the sorted moments are returned.

Input arguments:
    NumPorts    — Number of points to obtain on efficient frontier (default value is 10)
Output arguments:
    prsk        — Estimated efficient portfolio risks (standard deviation of returns)
    pret        — Estimated efficient portfolio returns 
=#
function plotFrontier(p, NumPorts = 10, optimiser = :Ipopt)
    if MINLP(p) == true
        if optimiser != :SCIP && optimiser != :Mosek && optimiser != :Xpress && optimiser != :Gurobi && optimiser != :Juniper
            optimiser = :Juniper
        end
    end
    A = estimateFrontier(p, NumPorts, optimiser)
    prsk = A[2]
    pret = A[3]

    arsk = sqrt.(diag(p.AssetCovar))
    aret = p.AssetMean

    if p.Name != ""
        ptitle = p.Name
    else
        ptitle = "Efficient Frontier"
    end
    fig = plot(prsk, pret, title = ptitle, xlabel = "Standard Deviation of Portfolio Returns", ylabel = "Mean of Portfolio Returns", label = "Efficient Frontier", xlim = (0, maximum(arsk) * 1.1), ylim = (0, maximum(pret)*1.1), legend = :bottomright)

    fig = scatter!([prsk[1]], [pret[1]], label = "MVP")
    
    if p.InitPort != Float64[]
        initprsk = sqrt.(p.InitPort' * p.AssetCovar * p.InitPort)
        initpret = dot(p.AssetMean, p.InitPort)
        println()
        fig = scatter!([initprsk],[initpret], label = "Initial Portfolio")
    end

    fig = scatter!(arsk, aret, label = "Assets")
    display(fig)  
end

#=
plotFrontierAllocations plots asset weights along the efficient frontiers using stacked bars

    weights = plotFrontierAllocations(p, NumPorts = 10)

=#
function plotFrontierAllocations(p, NumPorts = 10, optimiser = :Ipopt)
    if MINLP(p) == true
        if optimiser != :SCIP && optimiser != :Mosek && optimiser != :Xpress && optimiser != :Gurobi && optimiser != :Juniper
            optimiser = :Juniper
        end
    end
    weights = estimateFrontier(p, NumPorts, optimiser)[1]
    if p.Name != ""
        ptitle = p.Name
    else
        ptitle = "Efficient Frontier Allocations"
    end
    if getAssetList(p) == Symbol[]
        setAssetList(p)
    end
    ticklabel = "P" .* string.(collect(1:NumPorts))
    ticklabel[1] = "MVP"
    ticklabel[1NumPorts] = "MRP"
    fig = groupedbar(weights', bar_position = :stack, bar_width=1.0, title = ptitle, xlabel = "Portfolios", xticks=(1:NumPorts, ticklabel), ylabel = "Asset Allocation", label = permutedims(String.(p.AssetList)))
    display(fig)
end


#=
estimateMaxSharpeRatio estimates efficient portfolio to maximize Sharpe ratio

    (weights, risk, return, shaperatio) = estimateMaxSharpeRatio(p, optimiser = :Ipopt)

The estimateMaxSharpeRatio function maximizes the Sharpe ratio among portfolios on the efficient frontier. In the case of Portfolio with a risk-free asset, there are multiple efficient portfolios that maximize the Sharpe ratio on the capital asset line. 

The risk-free rate is obtained from the property RiskFreeRate in the Portfolio object. If you leave the RiskFreeRate unset, it is assumed to be 0. If the max return of portfolio is less than the RiskFreeRate, the solution is set as pwgt at max return and the resulting Sharpe ratio will be negative.
=#
function estimateMaxSharpeRatio(p, optimiser = :Ipopt)
    if MINLP(p) == true
        if optimiser != :SCIP && optimiser != :Mosek && optimiser != :Xpress && optimiser != :Gurobi && optimiser != :Juniper
            optimiser = :Juniper
        end
    end
    q = deepcopy(p)
    addEquality(q, q.AssetMean, 1)
    if MINLP(q) == true
        (model, w, v) = set_model(q, optimiser) 
        if q.LowerBudget == q.UpperBudget
            delete(model, model[:budget])
        else
            error("At present, maximising the Sharpe ratio is only available for portfolios with identical lower and upper budget constraints")
        end       
        f1(z,v) = z' * q.AssetCovar * z
        weights = optimise_model(model, w, v, :min, f1)[1]
        weights = weights / sum(weights)
    else
        (model, w) = set_model(q, optimiser)
        if q.LowerBudget == q.UpperBudget
            delete(model, model[:budget])
        else
            error("At present, maximising the Sharpe ratio is only available for portfolios with identical lower and upper budget constraints")
        end      
        f2(z) = z' * q.AssetCovar * z
        weights = optimise_model(model, w, :min, f2)[1]
        weights = weights / sum(weights)
    end
    if isnan(q.RiskFreeRate) == true 
        q.RiskFreeRate = 0.0
    end

    return weights, sqrt(weights' * q.AssetCovar * weights), q.AssetMean' * weights, ((q.AssetMean .- q.RiskFreeRate)' * weights ) / sqrt(weights' * q.AssetCovar * weights)
end


#=
plotSharpeRatio plots the maximum Sharpe ratio portfolio and the efficient frontier

    plotSharpeRatio(p, NumPorts = 10, optimiser = :Ipopt)

=#
function plotSharpeRatio(p, NumPorts = 10, optimiser = :Ipopt)
    if MINLP(p) == true
        if optimiser != :SCIP && optimiser != :Mosek && optimiser != :Xpress && optimiser != :Gurobi && optimiser != :Juniper
            optimiser = :Juniper
        end
    end
    A = estimateFrontier(p, NumPorts, optimiser)
    prsk = A[2]
    pret = A[3]
    
    (weights, riskShapeRatio, returnSharpeRatio, SharpeRatio) = estimateMaxSharpeRatio(p, optimiser)

    arsk = sqrt.(diag(p.AssetCovar))
    aret = p.AssetMean

    if p.Name != ""
        ptitle = p.Name
    else
        ptitle = "Maximum Sharpe Ratio Portfolio ($(round(SharpeRatio, digits = 4)))"
    end

    fig = plot(prsk, pret, title = ptitle, xlabel = "Standard Deviation of Portfolio Returns", ylabel = "Mean of Portfolio Returns", label = "Efficient Frontier", xlim = (0, maximum(arsk)*1.1), ylim = (0, maximum(pret)*1.1), legend = :topright)

    fig = scatter!([prsk[1]], [pret[1]], label = "Minimum Risk Portfolio")
    fig = scatter!([riskShapeRatio], [returnSharpeRatio], label = "Maximum Sharpe Ratio Portfolio")
 
    if p.InitPort != Float64[]
        initprsk = sqrt.(p.InitPort' * p.AssetCovar * p.InitPort)
        initpret = dot(p.AssetMean, p.InitPort)
        println()
        fig = scatter!([initprsk],[initpret], label = "Initial Portfolio")
    end

    fig = scatter!(arsk, aret, label = "Assets")
    
    display(fig)  
end


#=
computePortSharpeRatio estimates the Sharpe ratio of given portfolio weights

    psharpe = computePortSharpeRatio(p,pwgt)

The risk-free rate is obtained from the property RiskFreeRate in the Portfolio object. If you leave the RiskFreeRate unset, it is assumed to be 0.
=#
function computePortSharpeRatio(p, pwgt)
    if isnan(p.RiskFreeRate) == true
        return dot(p.AssetMean, pwgt)[1] / sqrt.(pwgt' * p.AssetCovar * pwgt)[1] # Assume risk-free rate = 0
    else
        return (dot(p.AssetMean, pwgt) .-  p.RiskFreeRate)[1]/ sqrt.(pwgt' * p.AssetCovar * pwgt)[1]
    end
end


#=
computePortMoments	estimates moments (std and mean) of portfolio returns

    (pstd,pret) = computePortMoments(p,pwgt)

=#
function computePortMoments(p,pwgt)
    pstd = sqrt.(pwgt' * p.AssetCovar * pwgt)[1]
    pret = dot(p.AssetMean, pwgt)[1]
    return pstd, pret
end

#=
computePortReturn estimate mean of portfolio returns

    pret = computePortReturn(p,pwgt)

=#
function computePortReturn(p,pwgt)
    pret = dot(p.AssetMean, pwgt)[1]
    return pret
end


#=
computePortRisk estimate portfolio risk according to risk proxy associated with corresponding object

    pstd = computePortRisk(p,pwgt)

=#
function computePortRisk(p,pwgt)
    pstd = sqrt.(pwgt' * p.AssetCovar * pwgt)[1]
    return pstd
end


#=
portfolioDiff shows purchases and sales relative to initial portfolio. If no initial portfolio is specified in p.InitPort, that value is assumed to be 0.

    portfolioDiff(pwgt)

=#
function portfolioDiff(pwgt)
    if p.InitPort != Float64[]
        ppurchases = max.(0, pwgt - p.InitPort)
        psales = max.(0, p.InitPort - pwgt)
    else
        ppurchases = pwgt
        psales = zeros(length(pwgt))
    end
    return ppurchases, psales
end

