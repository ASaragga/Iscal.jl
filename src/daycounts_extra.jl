

"""
    NL365()

**NL/365** or **Actual/365 Japanese** day count convention. 

Number of days in a period is equal to the actual number of days, except for leap days (29th February) which are ignored

The year fraction is
```math
\\frac{\\text{# of days excluding 29 February}}{365}
```

# Reference
"""
function yearfrac(startdate::Date, enddate::Date, dc::NL365)
    if startdate > enddate
        return -yearfrac(enddate, startdate, dc)
    end
    if month(enddate) == 2 && day(enddate) == 29
        # Since the startdate is included and the enddate is excluded, when the enddate is 29th February one should't consider this ocorrence in num29feb(startdate,enddate).
        if enddate == startdate
            return 0.0
        else 
            return (Dates.value(enddate - startdate) - num29feb(startdate, enddate) + 1)/365
        end
    else
        return (Dates.value(enddate - startdate) - num29feb(startdate, enddate))/365
    end
end

"""
    Bus252(calendar::HolidayCalendar)
    Bus252() - default Bus252(USNYSE()) 

**Bus/252** day count convention. The number of days in a period is equal to the actaul number of business days, depending on some calendar of holidays. There are available the following calendars:
* USSettlement()
* CanadaSettlement()
* UKSettlement() 
* BRSettlement()
* USGovernmentBond() 
* USNYSE()
* TARGET2()
* CanadaTSX() 
* AustraliaASX()
* BrazilExchange()
* WeekendsOnly() (excludes weekends only, does not exclude holydays) 
One obtains Matlab bus/252 day count basis using Bus252(USNYSE()) or Bus252()

The year fraction is
```math
\\frac{\\text{# of business days}}{252}
```

# Reference
"""
Bus252() = Bus252(USNYSE())
function yearfrac(startdate::Date, enddate::Date, dc::Bus252)
    if startdate > enddate
        return -yearfrac(enddate, startdate, dc)
    elseif startdate == enddate 
        return 0.0
    end 
    return bdayscount(dc.calendar, startdate, enddate)/252 
end



"""
    Thirty360PSA()

**30/360 (PSA)** or **30/360 (BMA)** day count convention (PSA, the Public Securites Association, was the predecessor of BMA, the Bond Market Association).

The year fraction is computed as:
```math
\\frac{360 \\times (y_2 - y_1) + 30 \\times (m_2 - m_1) + (d_2 - d_1)}{360}
```
where
- ``y_1`` and ``y_2`` are the years of the start and end date, respectively.
- ``m_1`` and ``m_2`` are the months of the start and end date, respectively.
- ``d_1`` is the day of the month at the start date, unless it is:
  * the last day of February, or
  * the 31st day of the month,
  in which case it is 30.
- ``d_2`` is the day of the month at the end date, unless it is:
  * the 31st day of the month, and 
  * ``d1`` is 30th or 31st day of the month, 
  in which case it is 30.

# Reference
"""
const Thirty360BMA = Thirty360PSA
function yearfrac(startdate::Date, enddate::Date, dc::Thirty360PSA)
    if startdate > enddate
        return -yearfrac(enddate, startdate, dc)
    elseif startdate == enddate 
        return 0.0
    end
    dy = year(enddate)-year(startdate)
    m1 = month(startdate)
    m2 = month(enddate)
    d1 = day(startdate)
    d2 = day(enddate)
    
    if startdate == lastdayofmonth(startdate) && m1 == 2 || d1 == 31
        d1 = 30
    end
    if day(enddate)==31 && d1 >= 30
        d2 = 30
    end
    return thirty360(dy,m2-m1,d2-d1)  # Defined in DayCounts.jl: thirty360(dy,dm,dd) = (360*dy + 30*dm + dd)/360
end


"""
    Thirty360SIA(eom:Bool)
    Thirty360SIA() - default Thirty360SIA(true)

**30/360 (SIA)**, **30U/360** or **30US/360 day count convention.

The year fraction is computed as:
```math
\\frac{360 \\times (y_2 - y_1) + 30 \\times (m_2 - m_1) + (d_2 - d_1)}{360}
```
where
- ``y_1`` and ``y_2`` are the years of the start and end date, respectively.
- ``m_1`` and ``m_2`` are the months of the start and end date, respectively.
- ``d_1`` is the day of the month at the start date, unless it is:
  * the last day of February, or
  * the 31st day of the month,
  in which case it is 30.
- ``d_2`` is the day of the month at the end date, unless it is:
  * the 31st day of the month, and ``d1`` is 30th or 31st day of the month, or
  * both dates were originally the last day of February, 
  in which case is it is 30.  
These rules assume that the security follows the end-of-month (eom) rule. If the security does not follow the end-of-month rule, then it returns the same as the 30/360 (also known as BondBasis or 30/360 ISDA) convention should be used instead.

# Reference
"""
Thirty360SIA() = Thirty360SIA(true)
const ThirtyU360 = ThirtyUS360 = Thirty360SIA

function yearfrac(startdate::Date, enddate::Date, dc::Thirty360SIA)
    if startdate > enddate
        return -yearfrac(enddate, startdate, dc)
    end
    dy = year(enddate)-year(startdate)
    m1 = month(startdate)
    m2 = month(enddate)
    d1 = day(startdate)
    d2 = day(enddate)
    if dc.eom == true
        if startdate == lastdayofmonth(startdate) && m1 == 2
            d1 = 30
            if enddate == lastdayofmonth(enddate) && m2 == 2
                d2 = 30
            end
        end
    end
    if d1 >= 30
        d1 = 30
        if d2 >= 30
            d2 = 30
        end
    end
    return thirty360(dy,m2-m1,d2-d1)  # Defined in DayCounts.jl: thirty360(dy,dm,dd) = (360*dy + 30*dm + dd)/360
end


"""
    ThirtyEPlus360()

**30E+/360** day count convention.

The year fraction is computed as:
```math
\\frac{360 \\times (y_2 - y_1) + 30 \\times (m_2 - m_1) + (d_2 - d_1)}{360}
```
where
- ``y_1`` is the year of the start date.
- ``y_2`` is the year of the end date, unless the date is 31st of December, in which case ``y_2`` will be changed to the following year (see ``m_2`` and ``d_2`` rules bellow).
- ``m_1`` is the month of the start date. 
- ``m_2``is the month of the end date, unless the day of that month is the 31st, in which case ``m_2`` will be changed to the following month (see ``d_2`` rule bellow).
- ``d_1`` is the day of the month at the start date, unless it is the 31st day of the month, in which case it is 30.
- ``d_2`` is the day of the month at the end date, unless it is the 31st of a month, in which case it will be changed to the 1st day of the next month.

# Reference
"""
function yearfrac(startdate::Date, enddate::Date, dc::ThirtyEPlus360)
    if startdate > enddate
        return -yearfrac(enddate, startdate, dc)
    elseif startdate == enddate
        return 0.0
    end
    y1 = year(startdate)
    y2 = year(enddate)
    m1 = month(startdate)
    m2 = month(enddate)
    d1 = day(startdate)
    d2 = day(enddate)
    if d1 == 31
        d1 = 30
    end
    if d2 == 31
        if m2 == 12
            d2 = 1
            m2 = 1
            y2 += 1
        else
            d2 = 1
            m2 +=1
        end
    end
    return thirty360(y2-y1,m2-m1,d2-d1)  # Defined in DayCounts.jl: thirty360(dy,dm,dd) = (360*dy + 30*dm + dd)/360
end


"""
    ActualActualL(frequency::Int)
    ActualActualL() - default ActualActualL(1)

**Actual/365L** or **Actual/365 (Leap Year)** day count convention.

The year fraction is computed as:
```math
\\frac{\\text{# of days}}{D}
```
where
* If the frequency is annual, ``D`` is 366 if the period contais February 29th, otherwise it is 365.
* If the frequency is not annual, ``D`` is 366 if the period end date is in a leap year, otherwise it is 365.

# Reference
- 2006 ISDA definitions 4.16i and ICMA rule 251.1(i) part 2.
"""
ActualActualL() = ActualActualL(1)   # Default is frequency equal to 1.
function yearfrac(startdate::Date, enddate::Date, dc::ActualActualL)
    if startdate > enddate
        return -yearfrac(enddate, startdate, dc)
    end
    y2 = year(enddate)
    if dc.frequency == 1
        if num29feb(startdate, enddate) > 0
            return Dates.value(enddate-startdate)/366
        else
            return Dates.value(enddate-startdate)/365
        end
    else
        if isleapyear(y2)
            return Dates.value(enddate-startdate)/366
        else 
            return Dates.value(enddate-startdate)/365
        end
    end
end


"""
    Actual365A

**Actual/365A** or **Actual/365 (Actual)** day count convention.

The year fraction is computed as:
```math
\\frac{\\text{# of days}}{D}
```
where
* ``D`` is 366 if a leap day is contained, or 
* ``D`` is 365 if not.

# Reference
"""
function yearfrac(startdate::Date, enddate::Date, dc::Actual365A)
    if startdate > enddate
        return -yearfrac(enddate, startdate, dc)
    end
    if num29feb(startdate, enddate) > 0
        return Dates.value(enddate-startdate)/366
    else
        return Dates.value(enddate-startdate)/365
    end
end


"""
    ActualActualY

**Actual/Actual (Year)** day count convention.

The year fraction is computed as:
```math
W+\\frac{N}{D}
```
where
- ``W`` is the number of whole years if the period is over one year.  
- ``N`` The numerator is the actual number of days in the **remaining** period, which is computed by adding a number of years (if the period is over one year) to the start date to reduce the remaining period to less than a year. If the start date is February 29th, then each time a year is added the last valid day in February is chosen.
- ``D`` is the actual number of days in the year from the adjusted start date.

# Reference
"""
function yearfrac(startdate::Date, enddate::Date, dc::ActualActualY)
    if startdate > enddate
        return -yearfrac(enddate, startdate, dc)
    elseif startdate == enddate
        return 0.0
    end

    incdate = startdate
    incw = 0
    if startdate == lastdayofmonth(startdate) && day(startdate) == 29
        while incdate < enddate
            global save_date = incdate
            global save_w = incw
            incdate = lastdayofmonth(Date(year(incdate)+1,month(incdate),28))
            incw += 1
        end
    else
        while incdate < enddate
            global save_date = incdate
            global save_w = incw
            incdate = Date(year(incdate)+1, month(incdate), day(incdate))
            incw += 1
        end
    end
    n = Dates.value(enddate-save_date)
    if save_date == lastdayofmonth(save_date) && day(save_date) == 29
        d = Dates.value(Date(year(save_date)+1, month(save_date), 28)-save_date)
    else    
        d = Dates.value(Date(year(save_date)+1, month(save_date), day(save_date))-save_date)    
    end
    return save_w + n/d
end


"""
    Actual365ICMA()

**Actual/365 (ICMA)** day count convention. Returns the same value as the **Actual365Fixed** day count convention when computing the yearfrac() function. However the computation of accrued interest is diferent. 

The year fraction is computed as:
```math
\\frac{\\text{# of days}}{365}
```

# Reference
 - 2006 ISDA definitions, §4.16 (d)
"""
function yearfrac(startdate::Date, enddate::Date, ::Actual365ICMA)
    return Dates.value(enddate-startdate)/365
end



"""
    Actual360ICMA()

**Actual/360 ICMA** day count convention. Returns the same value as the **Actual/360** day count convention when computing the yearfrac() function. However the computation of accrued interest is diferent. 

The year fraction is computed as:
```math
\\frac{\\text{# of days}}{360}
```

# Reference
 - 2006 ISDA definitions, §4.16 (e)

"""
function yearfrac(startdate::Date, enddate::Date, ::Actual360ICMA)
    return Dates.value(enddate-startdate)/360
end




"""
    Thirty360ICMA()

**30/360 ICMA** day count convention. Returns the same value as the **ThirtyE360** day count convention when computing the yearfrac() function. However the computation of accrued interest is diferent. 

The year fraction is computed as:
```math
\\frac{360 \\times (y_2 - y_1) + 30 \\times (m_2 - m_1) + (d_2 - d_1)}{360}
```
where
- ``y_1`` and ``y_2`` are the years of the start and end date, respectively.
- ``m_1`` and ``m_2`` are the months of the start and end date, respectively.
- ``d_1`` is the day of the month at the start date, unless it is 31st day of the month,
  in which case it is 30.
- ``d_2`` is the day of the month at the end date,  unless it is 31st day of the month,
  in which case it is 30.

# Reference
 - 2006 ISDA definitions, §4.16 (g)
"""
function yearfrac(startdate::Date, enddate::Date, dc::Thirty360ICMA)
    if startdate > enddate
        return -yearfrac(enddate, startdate, dc)
    end
    dy = year(enddate)-year(startdate)
    dm = month(enddate)-month(startdate)
    d1 = min(day(startdate),30)
    d2 = min(day(enddate),30)
    return thirty360(dy,dm,d2-d1)
end


"""
    ActualActualAFB()

**Actual/Actual AFB** day count convention.

# Reference
"""
function yearfrac(startdate::Date,  enddate::Date, dc::ActualActualAFB)
    if (startdate > enddate)
        return -yearfrac(enddate, startdate, dc)
    elseif startdate == enddate
        return 0.0
    end

    newD2 = enddate
    temp = enddate
    sum = 0.0
    while (temp > startdate)
	    temp = newD2 - Year(1)
	    if (day(temp) == 28 && month(temp) == 2 && isleapyear(temp))
	        temp + Day(1)
	    end
	    if (temp>=(startdate))
	        sum += 1.0
	        newD2 = temp
	    end
    end

    den = 365.0

    if isleapyear(newD2)
	    if (newD2 > Date(year(newD2), February, 29) && startdate < Date(year(newD2), February, 29))
	        den += 1.0
	    end
    elseif (isleapyear(startdate))
	    if (newD2 > (ymd(year(startdate), February, 29)) && startdate < (ymd(year(startdate), February, 29)))
	        den += 1.0
	    end
    end
    return sum + Dates.value(newD2 - startdate) / den
end

# Function testing if a day count convention is a SIA base 
function issia(dc)
    if isa(dc, ActualActualMatlab) == true || isa(dc, Thirty360SIA) == true || isa(dc, Actual360) == true || isa(dc, Actual365Fixed) == true || isa(dc, Thirty360PSA) == true || isa(dc, Thirty360) == true || isa(dc, ThirtyE360) == true || isa(dc, NL365) == true 
        return true
    else
        return false
    end
end

# Function testing if a day count convention is a ISCMA base 
function isisma(dc)
    if isa(dc, ActualActualICMA) == true || isa(dc, Actual360ICMA) == true || isa(dc, Actual365ICMA) == true || isa(dc, Thirty360ICMA) == true || isa(dc, ActualActualISDA) == true 
        return true
    else
        return false
    end
end

