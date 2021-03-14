#= Currency and Price Conversion (6/7)
OK  cur2frac        - Decimal currency values to fractional values
OK  cur2str         - Bank-formatted text
    dec2thirtytwo	- Decimal to thirty-second quotation
OK  frac2cur        - Fractional currency value to decimal value
OK  thirtytwo2dec	- Thirty-second quotation to decimal
OK  todecimal	    - Fractional to decimal conversion
OK  toquoted        - Decimal to fractional conversion
=#

#=
cur2frac(Decimal,Denominator) converts decimal currency values to fractional values. Fraction is returned as a character vector.

    Fraction = cur2frac(Decimal,Denominator)

Input arguments:
    Decimal     — Decimal currency value
    Denominator — Denominator of the fractions
Output arguments:
    Fraction — Fractional values
=#
function cur2frac(Decimal, Denominator) # f = cur2frac(p, d)
    # Check for negative values
    if Decimal < 0
        coeff = -1
    else
        coeff = 1
    end

    frac = abs(round(rem(Decimal, 1) * Denominator))
    if frac == Denominator
        crt = coeff
        frac = 0
    else
        crt = 0
    end

    if frac < 10 && Denominator >10
        #remain = string(frac * 10^floor(log10(Denominator)))
        remain = "0" * string(convert(Int64,frac))
    else
        remain = string(convert(Int64,frac))
    end

    return string(coeff * floor(Int64,abs(Decimal)+crt)) * "." * remain
end

#=
cur2str returns the given value in bank format. The output format for BankText is a numerical format with dollar sign prefix, two decimal places, and negative numbers in parentheses; for example, ($123.45).

    BankText = cur2str(Value, Digits = 2, Currency = :USD)

Input arguments:
    Digits — Number of significant digits. A negative Digits rounds the value to the left of the decimal point. For example, cur2str(-826444.4456,-2) returns ($826400)
    Currency (not in Matlab) - Defines Value in terms of a currency. Possible symbol values are:
        :USD
        :EUR
        :GBP
        :JPY
        or, alternativelly, a char. For example, cur2str(-826444.4456,3,'₩'), where ₩ is the currency symbol for the Korean won, returns (₩826444.446)
=#
function cur2str(Value, Digits = 2, Currency::Union{Symbol, Char} = :USD) # cur2str(p,d)
    if Currency == :USD
        symb = '$'
    elseif Currency == :EUR 
        symb = '€'
    elseif Currency == :GBP
        symb = '£'
    elseif Currency == :JPY
        symb = '¥'
    else
        symb = Currency 
    end
    if Digits > 0
        if Value > 0
            st = symb * string(round(Value, digits = Digits))
        else
            st = '(' * symb * string(round(abs(Value), digits = Digits)) * ')'
        end
    elseif Digits < 0
        a = round(abs(Value) * 10.0^Digits) 
        b = string(convert(Int64,a))
        if Value < 0 # Negative numbers are enclosed in parenthesis
            pad = ""
            for i = 1:abs(Digits)
                pad = pad * "0"
            end
            st = '(' * symb * b * pad * ')'
        else
            it = symb * b
            st = it * string(zeros(1,abs(Digits)))
        end
    else
        if Value > 0
            st = symb * string(convert(Int64,round(Value, digits = Digits))) * '.'
        else
            st = '(' * symb * string(convert(Int64,round(abs(Value),digits = Digits))) * '.' * ')'
        end
    end
    return st
end
#=
dec2thirtytwo changes a decimal price quotation for a bond or bond future to a fraction with a denominator of 32.

    (OutNumber, Fraction) = dec2thirtytwo(InNumber, Accuracy) 

Input arguments:
    InNumber — Input number
    Accuracy — Rounding, with numeric values of 1, 2, 4 or 10. The values are: 
         1 (default) round down to nearest thirty second, 
         2 nearest half
         4 nearest quarter
        10 nearest decile
Output arguments:
    OutNumber — Output number which is InNumber rounded downward to closest integer 
    Fractions — Fractional part in units of thirty-second. The Fractions output conforms to accuracy as prescribed by the input Accuracy.
=#


#=
frac2cur converts a fractional currency value to a decimal value. Fraction is the fractional currency value input as a character vector, and Denominator is the denominator of the fraction.

    Decimal = frac2cur(Fraction, Denominator) 

=#
function frac2cur(Fraction, Denominator) # frac2cur(p, f)
    # Check for '-' signs
    if occursin('-', Fraction) == true
        coeff = -1
    else
        coeff = 1
    end
    
    # Find the decimal point
    b = findfirst(isequal('.'), Fraction)

    if b === nothing 
        decimal = parse(Int64, Fraction)
    else
        n = length(Fraction)
        decimal = parse(Float64,Fraction[1:b-1]) + coeff * parse(Float64,Fraction[b+1:n])/Denominator
    end
    return decimal
end

#=
thirtytwo2dec changes the price quotation for a bond or bond future from a fraction with a denominator of 32 to a decimal.

    OutNumber = thirtytwo2dec(InNumber, InFraction)

Input arguments:
    InNumber - Input number, representing price without the fractional components
    InFraction - Fractional portions of each element in InNumber
Output arguments:
    OutNumber - Output number that represents sum of InNumber and InFraction, returned as a decimal.
=#
function thirtytwo2dec(InNumber,InFraction)
    OutNumber = InNumber + InFraction/32
end


#=
todecimal returns the decimal equivalent, usddec, of a security whose price is normally quoted as a whole number and a fraction (quote). 

    usddec = todecimal(quote,fracpart)

fracpart indicates the fractional base (denominator) with which the security is normally quoted (default = 32).
=#
function todecimal(usdquote, fracpart = 32)
    # Get the integer part of the figures.
    whole = floor(usdquote)

    # Get the fractional part of the figures and convert them to integer to become the numerator of the fraction.
    numerator = round((usdquote-whole) * 100)

    # Calculate the decimal equivalent of the fractions
    decpart = numerator/fracpart

    # Combine with the integer part to create the figures' decimal equivalents.
    return whole + decpart
end

#=
toquoted returns the fractional equivalent, quote, of the decimal figure, usddec, based on the fractional base (denominator), fracpart. 

    quote = toquoted(usddec, fracpart = 32)

The fractional bases are the ones used for quoting equity prices in the United States (denominator 2, 4, 8, 16, or 32). If fracpart is not entered, the denominator 32 is assumed.

Note: The convention of using . (period) as a substitute for : (colon) in the output is adopted from Excel® software.
=#
function toquoted(usddec, fracpart = 32)
    # Get the integer part of the figures.
    whole = floor(usddec);
    
    # Get the decimal part of the figures.
    decpart = usddec - whole;
    
    # Calculate the integer numerator for the fractional representations. Should I use FLOOR, CEIL, or ROUND here?  Does it matter with regards to bid/ask prices?
    numerator = round(decpart * fracpart)
    
    # Combine with the integer part to create the figures' fractional equivalents.
    usdquote = whole + (numerator/100);
end