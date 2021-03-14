"""
    ActualActualMatlab()

**Actual/Actual (Matlab)** day count convention, as computed via Matlab `yearfrac` function with the option `0` selected.

The year fraction is
```math
\\frac{\\text{# of days}}{D}
```
where:
* if start date is in a leap year, and before or falling on 29 February, then ``D = 366``,
* otherwise ``D = 365`` 

# Reference
- [Matlab documentation for `yearfrac`](https://www.mathworks.com/help/finance/yearfrac.html#bu_82me-1-Basis)
"""
function yearfrac(startdate::Date, enddate::Date, dc::ActualActualMatlab)
    if startdate > enddate
        return -yearfrac(enddate, startdate, dc)
    end
    y1 = year(startdate)
    m1 = month(startdate)
    d1 = day(startdate)
    if startdate == lastdayofmonth(startdate) && d1 == 29
        return  Dates.value(enddate-startdate)/Dates.value(Date(y1+1, 3, 1)-startdate)
    else
        return Dates.value(enddate-startdate)/Dates.value(Date(y1+1, m1, d1)-startdate)
    end
end
