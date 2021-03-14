#=
Investment Performance Metrics (10i/10)
OK  elpm            - Compute expected lower partial moments for normal asset returns
OK  emaxdrawdown    - Compute expected maximum drawdown for Brownian motion
OK  inforatio	    - Calculate the information ratio
OK  lpm             - Compute sample lower partial moments of data
OK  maxdrawdown     - Compute maximum drawdown for one or more price series
OK  portalpha	    - Compute risk-adjusted alphas and returns
OK  sharpe          - Compute the Sharpe ratio
    =================================================================================
OK  omega           - Compute the Omega investment ratio  
OK  sortino         - Compute the Sortino investment ratio 
OK  upside          - Computes the Upside Potential investment ratio 
=#


#=
emaxdrawdown computes expected maximum drawdown for a Brownian motion

    EDD = emaxdrawdown(Mu, Sigma, T)

Input arguments:
    Mu      - The drift term of a Brownian motion with drift.
    Sigma   - The diffusion term of a Brownian motion with drift.
    T       - A time period of interest or a vector of times.
Outputs arguments:
    EDD - Expected maximum drawdown for each time period in T.

Notes: If a geometric Brownian motion with stochastic differential equation
    dS(t) = Mu0 * S(t) * dt + Sigma0 * S(t) * dW(t) ,
convert to the form here by Ito's Lemma with X(t) = log(S(t)) such that
    Mu = M0 - 0.5 * Sigma0^2
    Sigma = Sigma0
=#
function emaxdrawdown(Mu, Sigma, T)
    if Sigma < eps()
	    if Mu >= 0.0
		    EDD = 0.0
	    else
		    EDD = - Mu * T
	    end
	    return EDD
    end
    if abs(Mu) <= eps()
        EDD = (sqrt(pi/2) * Sigma) * sqrt(T)
    elseif Mu > eps()
        Alpha = Mu/(2.0 * Sigma^2)
        EDD = emaxddQp(Alpha * Mu * T)/Alpha
    else
        Alpha = Mu/(2.0 * Sigma * Sigma)
        EDD = - emaxddQn(Alpha * Mu * T)/Alpha
    end
    return EDD
end

function emaxddQp(x)
    A = [	0.0005;		0.001;		0.0015;		0.002;		0.0025;		0.005;
            0.0075;		0.01;		0.0125;		0.015;		0.0175;		0.02;
            0.0225;		0.025;		0.0275;		0.03;		0.0325;		0.035;
            0.0375;		0.04;		0.0425;		0.045;		0.0475;		0.05;
            0.055;		0.06;		0.065;		0.07;		0.075;		0.08;
            0.085;		0.09;		0.095;		0.1;		0.15;		0.2;
            0.25;		0.3;		0.35;		0.4;		0.45;		0.5;
            1.0;		1.5;		2.0;		2.5;		3.0;		3.5;
            4.0;		4.5;		5.0;		10.0;		15.0;		20.0;
            25.0;		30.0;		35.0;		40.0;		45.0;		50.0;
            100.0;		150.0;		200.0;		250.0;		300.0;		350.0;
            400.0;		450.0;		500.0;		1000.0;		1500.0;		2000.0;
            2500.0;		3000.0;		3500.0;		4000.0;		4500.0;		5000.0 ]
    B = [	0.01969;	0.027694;	0.033789;	0.038896;	0.043372;	0.060721;
            0.073808;	0.084693;	0.094171;	0.102651;	0.110375;	0.117503;
            0.124142;	0.130374;	0.136259;	0.141842;	0.147162;	0.152249;
            0.157127;	0.161817;	0.166337;	0.170702;	0.174924;	0.179015;
            0.186842;	0.194248;	0.201287;	0.207999;	0.214421;	0.220581;
            0.226505;	0.232212;	0.237722;	0.24305;	0.288719;	0.325071;
            0.355581;	0.382016;	0.405415;	0.426452;	0.445588;	0.463159;
            0.588642;	0.668992;	0.72854;	0.775976;	0.815456;	0.849298;
            0.878933;	0.905305;	0.92907;	1.088998;	1.184918;	1.253794;
            1.307607;	1.351794;	1.389289;	1.42186;	1.450654;	1.476457;
            1.647113;	1.747485;	1.818873;	1.874323;	1.919671;	1.958037;
            1.991288;	2.02063;	2.046885;	2.219765;	2.320983;	2.392826;
            2.448562;	2.494109;	2.532622;	2.565985;	2.595416;	2.621743 ]       
    if x > 5000
        Q = 0.25 * log(x) + 0.49088
    elseif x < 0.0005
        Q = 0.5 * sqrt(pi*x)
    elseif x >= 0.0005 && x <= 5000
        interp_linear = LinearInterpolation(A, B)
        Q = interp_linear(x)
    end
    return Q
end

function emaxddQn(x)    
    A = [	0.0005;		0.001;		0.0015;		0.002;		0.0025;		0.005;
            0.0075;		0.01;		0.0125;		0.015;		0.0175;		0.02;
            0.0225;		0.025;		0.0275;		0.03;		0.0325;		0.035;
            0.0375;		0.04;		0.0425;		0.045;		0.0475;		0.05;
            0.055;		0.06;		0.065;		0.07;		0.075;		0.08;
            0.085;		0.09;		0.095;		0.1;		0.15;		0.2;	
            0.25;		0.3;		0.35;		0.4;		0.45;		0.5;
            0.75;		1.0;		1.25;		1.5;		1.75;		2.0;
            2.25;		2.5;		2.75;		3.0;		3.25;		3.5;
            3.75;		4.0;		4.25;		4.5;		4.75;		5.0 ]
    B = [	0.019965;	0.028394;	0.034874;	0.040369;	0.045256;	0.064633;
            0.079746;	0.092708;	0.104259;	0.114814;	0.124608;	0.133772;
            0.142429;	0.150739;	0.158565;	0.166229;	0.173756;	0.180793;
            0.187739;	0.194489;	0.201094;	0.207572;	0.213877;	0.220056;
            0.231797;	0.243374;	0.254585;	0.265472;	0.27607;	0.286406;
            0.296507;	0.306393;	0.316066;	0.325586;	0.413136;	0.491599;
            0.564333;	0.633007;	0.698849;	0.762455;	0.824484;	0.884593;
            1.17202;	1.44552;	1.70936;	1.97074;	2.22742;	2.48396;
            2.73676;	2.99094;	3.24354;	3.49252;	3.74294;	3.99519;
            4.24274;	4.49238;	4.73859;	4.99043;	5.24083;	5.49882 ]
    if x > 5.0
        Q = x + 0.5
    elseif x < 0.0005
        Q = 0.5 * sqrt(pi*x)
    elseif x >= 0.0005 && x <= 5
        interp_linear = LinearInterpolation(A, B)
        Q = interp_linear(x)
    end
    return Q
end


#=
Calculate maximum drawdown for one or more price series

    [MaxDD, MaxDDIndex] = maxdrawdown(Data, Format)

Input arguments:
    Data - T x N matrix with T samples of N total equity time series with earliest data in row
    T(1,:) and most recent data in row T(end,:).
    Format - Indicates format of data. Options are:
		:return (default)   - Compute maximum drawdown as a maximum percentage drop from a peak
		:arithmetic         - Compute maximum drawdown of a Brownian motion with drift (differences of data from peak to trough)
		:geometric - Compute maximum drawdown of a geometric Brownian motion with drift (differences of log of data from peak to trough).
Output arguments:
    MaxDD - 1 x N vector with maximum drawdown for each of N time series.
    MaxDDIndex - 2 x N vector of latest start and earliest end indexes for each maximum drawdown period for each total equity time series, where the first row contains the start indexes and the second row contains the end indexes of each maximum drawdown period.
=#
function maxdrawdown(data, format = :return)
    if format == :geometric
        data = log.(data)
    end
    # each peak can potentially be the start of a max drawdown
    peaks = accumulate(max,data)     
    if format == :return 
        drdn = (peaks - data)./peaks
    elseif format == :arithmetic || format == :geometric
        drdn = (peaks .- data)
    end
    MaxDD = drdn[1]
    MaxDDIndex = zeros(2,1)
    obs = size(drdn,1)
    for i in 2:obs
        if drdn[i] > MaxDD
            MaxDD = drdn[i]
            MaxDDIndex[2] = i
        end
    end
    k = convert(Int64,MaxDDIndex[2])
    while data[k-1] >= data[k] && k > 1
        MaxDDIndex[1] = k - 1
        k -= 1
    end
    return MaxDD, convert.(Int64, MaxDDIndex)
end

#= lpm computes sample lower partial moments from data

moment = lpm(data, mar = 0, order = 0)

Input arguments:
    data    - Observations of asset returns
	mar     - Scalar minimum acceptable return (default MAR = 0). This is a cutoff level of return such that all returns above MAR contribute nothing to the lower partial moment.
    order   - Non-negative integer moment orders. If no order specified, default order = 0, which is the shortfall probability. Although this function will work for non-integer orders and, in some cases, for negative orders, this falls outside customary usage.
Output arguments:
    Moment - Lower partial moments for a given order.
Notes: To compute upper partial moments, just reverse the signs of both Data and MAR (do not reverse the sign of the output). With lpm, you can compute various investment ratios such as Omega, Sortino, Kappa, and Upside Potential, where:
    Omega = lpm(-Data, -MAR, 1) ./ lpm(Data, MAR, 1)
    Sortino = (mean(Data) - MAR) ./ sqrt(lpm(Data, MAR, 2))
    Upside = lpm(-Data, -MAR, 1) ./ sqrt(lpm(Data, MAR, 2))
=#


function lpm(data, mar = 0, order = 0)
    m = size(data,1)
    p = length(order)
    if p > 1
        moment = zeros(p)
    elseif p == 1
        moment = 0
    end
    ii = 0
    for i = 1:m
        ii = ii + 1;
        if data[i] <= mar
            moment = moment .+ (mar - data[i]) .^ order
        end
    end
    moment = (1/ii) .* moment
    return moment
end

#=
omega computes the Omega investment ratio 

    omega(data, mar)

=#
function omega(data, mar)
    return lpm(-data, -mar, 1) / lpm(data, mar, 1)
end

#=
sortino computes the Sortino investment ratio

    sortino(data, mar)

=#
function sortino(data, mar) 
    return (mean(data) - mar) ./ sqrt(lpm(data, mar, 2))
end

#=
upside computes the Upside Potential investment ratio 

    upside(data, mar)

=#
function upside(data, mar)
    return lpm(-data, -mar, 1) ./ sqrt(lpm(data, mar, 2))
end

#=
elpm computes expected lower partial moments relative to a minimum acceptable return (mar) and moment orders (order)

    moment = elpm(mean, sigma, mar, order)

Input arguments:
    mar     - Minimum acceptable return. 'mar' is a cutoff level of return such that all returns above 'mar' contribute nothing to the lower partial moment.
    order - Expected lower partial moments. 

    Note: Use lower partial moments to examine what is colloquially known as “downside risk.” The main idea of the lower partial moment framework is to model moments of asset returns that fall below a minimum acceptable level of return. 
=#
function elpm(mean, sigma, mar = 0, order = 0)
    p = length(order)
    if p == 1
        MaxOrder = order
    elseif p > 1
        MaxOrder = order[1]
        for i in 1:p
            if order[i] > MaxOrder
                MaxOrder = order[i]
            end
        end
    end
    println(MaxOrder)
    moment = zeros(p)

	AllMoments = elpmrange(mean, sigma, mar, MaxOrder);
	for j = 1:p
		moment[j] = AllMoments[order[j]+1]
	end
    return moment
end

# elpmrange - Recursion to compute all expected lower partial moments for all orders 0, 1, ... , MaxOrder.
function elpmrange(mean, sigma, mar, MaxOrder)
    moments = zeros(MaxOrder + 1,1);
    if sigma == 0
	    if mean <= mar
		    moments = (mar - mean) .^ (0:MaxOrder)
	    end
    else
	    alpha = (mean - mar)/(sigma * sqrt(2))
	    moments[1] = 0.5 * erfc(alpha)
	    if MaxOrder == 0
		    return
	    end
	    moments[2] = sigma * exp(-alpha^2)/sqrt(2*pi) - (mean - mar) * moments[1]
	    if MaxOrder == 1
		    return
	    end
	    MM = mean - mar
	    SS = sigma^2

	    for i = 0:(MaxOrder - 2)
		    moments[i + 3] = (i + 1) * SS * moments[i + 1] - MM * moments[i + 2]
	    end
    end
    return moments
end


#=
inforatio compute the information ratio and tracking error for one asset relative to the benchmark

    (ratio, te) = inforatio(asset, benchmark)
Input arguments:
    asset       - Observations of asset returns
	benchmark   - Returns for a benchmark asset. The periodicity must be the same as the periodicity of Asset
Output arguments
	ratio - Information ratio for Asset series. Any series in Asset with a tracking error of zero will have a NaN value for its information ratio.
    te - Tracking errors, i.e., the standard deviation of Asset relative to Benchmark Returns
=#
function inforatio(asset, benchmark)
    asset = asset - benchmark
    te = std(asset, corrected = false)
    ratio = 1 ./ te
    ratio = ratio .* mean(asset)
    return ratio, te
end


#=
 portalpha compute risk-adjusted alphas and returns for one asset. Given one asset with NUMSAMPLES returns, a NUMSAMPLES vector of Benchmark returns, and either a scalar Cash return or a NUMSAMPLES vector of Cash returns, compute risk-adjusted alphas and returns for one or more methods specified by Choice.

    [Alpha, RAReturn] = portalpha(Asset, Benchmark, Cash, Choice)

 Input arguments:
    Asset - Observations of asset returns
 	Benchmark - Returns for a benchmark asset. The periodicity must be the same as the periodicity of asset
    Cash - Either a scalar return for a riskless asset or a vector of asset returns to be a proxy for a "riskless" asset. In either case, the periodicity must be the same as the periodicity of asset. If no value is supplied, the default value for Cash returns is 0. 
	Choice - Indicator for one or more measures to be computed from among a number of risk-adjusted alphas and return measures. The current list of choices is given in the following table:
		Code		Description
 		------		----------------------------------
 		:xs	    	Excess Return (no risk adjustment)
 		:sml		Security Market Line
 		:capm		Jensen's Alpha
 		:mm 		Modigliani & Modigliani
 		:gh1		Graham-Harvey 1
 		:gh2		Graham-Harvey 2
 		------		----------------------------------
 		:all		Compute all measures
    Choices are specified by the Code from the table (e.g., to select the Modigliani & Modigliani measure, Choice = 'mm'). A single choice is a single Code from the table. Multiple choices can be selected with a array of choice Codes (e.g., to select both Graham-Harvey measures, Choice = {:gh1, :gh2}). To select all choices, specify Choice = :all. If no value is supplied, the default choice is to compute the excess return with Choice = :xs.
Output arguments:
    Alpha - Risk-adjusted alphas for series in asset with each row corresponding to a specified measure in Choice.
    RAReturn - Risk-adjusted returns for each series in Asset with each row corresponding to a specified measure in Choice.
=#
function portalpha(Asset, Benchmark, Cash = 0.0, Choice = :xs)
    # Step 1- initialization
    if Choice == :all
        Choice = [:xs, :sml, :capm, :mm, :gh1, :gh2]
    end
    if isa(Choice, Symbol) == true
        m = 1
        Choice = [Choice]
    elseif isa(Choice, Vector{Symbol}) == true
        m = length(Choice)
    end
    Alpha = zeros(m)
    RAReturn = zeros(m)
    mF = 0
    sF = 0
    rFC = 0
    rFM = 0
    if isa(Cash, Real) == true
        Cash = Cash * ones(length(Asset))
    end
    # Step 2 - compute basic statistics
    mM = mean(Benchmark)
    sM = std(Benchmark)

    mC = mean(Cash)
    sC = std(Cash)

    if sM == 0 || sC == 0
        rMC = 0
    else
        rMC = cor(Benchmark, Cash)
    end

    mF = mean(Asset)
    sF = std(Asset)
    if sF == 0 || sC == 0
        rFC = 0
    else
        rFC = cor(Asset, Cash)
    end
    if sF == 0 || sM == 0
        rFM = 0
    else
        rFM = cor(Asset, Benchmark)
    end

    # Step 3 - compute measures

    ii = 0
    while ii < m
        ii = ii + 1
        if Choice[ii] == :xs
            Alpha[ii] = mF - mM
            RAReturn[ii] = mF
        elseif Choice[ii] == :sml
            if sM == sC
                warn("SML Non Existent!")
                Alpha[ii] = NaN
            else
                Alpha[ii] = (mF - mC) - (sF - sC) * (mM - mC)/(sM - sC)
            end
            RAReturn[ii] = mF - Alpha[ii]
        elseif Choice[ii] == :capm
            if maximum(Benchmark) == minimum(Benchmark)
                warn("CAPM Non Existent!")
                Alpha[ii] = NaN
            else
                Alpha[ii] = (mF - mC) - (rFM*sF/sM) * (mM - mC)
            end
            RAReturn[ii] = mF - Alpha[ii]
        elseif Choice[ii] == :mm
            if maximum(Asset) == minimum(Asset)
                warn("MM Non Existent!")
                Alpha[ii] = NaN
            else
                Alpha[ii] = (sM/sF) * (mF - mC) - (mM - mC)
            end
            RAReturn[ii] = Alpha[ii] + mM
        elseif Choice[ii] == :gh1
            L1 = sC*sC - rMC*sM*sC;
            L2 = sM*sM - rMC*sM*sC;
            if abs(L1 + L2) < eps()
                # Benchmark and Cash have equal risks with either +1 correlation or zero risk.
                warn("GH1 Non Existent")
                Alpha[ii] = NaN
            else
                Disc = L1*L1 + (L1 + L2) * (sF^2 - sC^2)
                if Disc < 0
                    warn("GH1 Non Existent Low Risk")
                    Alpha[ii] = NaN
                else
                    Factor = (L1 + sqrt(abs(Disc)))/(L1 + L2)
                    Alpha[ii] = (mF - mC) - Factor*(mM - mC);
                end
            end
            RAReturn[ii] = mF - Alpha[ii]
        elseif Choice[ii] == :gh2
            L1 = sC^2 - rFC * sF * sC
            L2 = sF^2 - rFC * sF * sC
            Disc = L1^2 + (L1 + L2) * (sM^2 - sC^2) 
            if abs(L1 + L2) < eps()
                # Asset and Cash have equal risks with either +1 correlation or zero risk.
                warn("GH2 Non Existent")
                Alpha[ii] = NaN
            elseif Disc < 0
                warn("GH2 Non Existent Low Risk")
                Alpha[ii] = NaN
            else
                Factor = (L1 + sqrt(abs(Disc)))/(L1 + L2)
                Alpha[ii] = Factor*(mF - mC) - (mM - mC);
            end
            RAReturn[ii] = Alpha[ii] + mM
        else
            error("Invalid Choice Element")
        end
    end
    return Alpha, RAReturn
end


#=
sharpe compute the Sharpe ratio for n asset

    ratio = sharpe(asset, cash)

Input arguments:
    asset - Observations of asset returns
    cash - Either a scalar return for a riskless asset or a vector of asset returns to be a proxy for a riskless asset. In either case, the periodicity must be the same as the periodicity of asset. If no value is supplied, the default value for Cash returns is 0.0
Output arguments:
    ratio - Sharpe ratio. If the asset standard deviation of returns equals to 0 it will have a NaN value for its Sharpe ratio
Note: If 'cash' is a vector, 'asset' and 'cash' need not have the same number of returns but must have the same periodicity of returns. Note that the classic Sharpe ratio assumes that 'cash' is riskless. In reality, a short-term cash rate is not necessarily riskless
=#
function sharpe(asset, cash = 0.0)
    if isa(cash, Real) == true 
        C0 = cash
    elseif isa(cash, Vector{<:Real})
        C0 = mean(cash)
    else
        error("Cash type invalid")
    end
    denom = std(asset, corrected = false)
    if maximum(asset) == minimum(asset)
        ratio = NaN
    else
        ratio = 1 ./ denom
        ratio = ratio .* (mean(asset) - C0)
    end
    return ratio
end
