#=
Date and Time Component (13/13):
OK  today	    - Current date
OK  day	        - Day of month
OK  eomdate	    - Last date of month
OK  hour	    - Hour of date or time
OK  lweekdate	- Date of last occurrence of weekday in month
OK  second	    - Seconds of date or time
OK  minute	    - Minute of date or time
OK  month	    - Month of date
OK  months	    - Number of whole months between dates
OK  nweekdate	- Date of specific occurrence of weekday in month
OK  weeknum	    - Week in year
OK  year	    - Year of date
OK  yeardays	- Number of days in year
=#


# today - Current date
today()

# eomdate	    - Last date of month
eomdate = lastdayofmonth

#= hour	        - Hour of date or time    
    hour(time)
=#

#= second	    - Seconds of date or time
    second(time)
=#

#= minute	    - Minute of date or time
    minute(time)
=#

#= month	        - Month of date
    month(date)
0#

# weeknum	    - Week in year
    weeknum = week

Matlab: weeknum('20 mar 2021') = 12. Return the ISO week date of a Date or DateTime as an Int64. Note that the first week of a year is the week that contains the first Thursday of the year, which can result in dates prior to January 4th being in the last week of the previous year. For example, week(Date(2005, 1, 1)) is the 53rd week of 2004.
week() 
# e.g. week(Date(2021,3,20)) == 11             
=#

weeknum = week

#= year	        - Year of date
    year(date)
=#


#= yeardays      - Number of days in year
    nd = yeardays(year, basis)
=#

function yeardays(year, basis)   # Errors at least for Actual365ICMA(), NL365()
    first   = Date(year,   1, 1)
    last    = Date(year+1, 1, 1)
    if basis isa ActualActualMatlab || basis isa ActualActualICMA || basis isa Actual365ICMA || basis isa ActualActualISDA
        return daysact(first, last)
    elseif basis isa Thirty360SIA || basis isa Actual360 || basis isa Thirty360PSA || basis isa Thirty360 || basis isa ThirtyE360 || basis isa Actual360ICMA || basis isa ThirtyE360ISMA
        return days360(first, last)
    elseif basis isa Actual365Fixed || basis isa NL365
        return days365(first, last)
    elseif basis isa Bus252
        return days252bus(first,last) # We use the default which is calender = NYSE
    end
end

function lweekdate(weekday,year,month)
    date = Date(year, month, daysinmonth(Date(year,month,1)))
    while dayabbr(date) != weekday
        date -= Day(1)
    end
    return date
end

function months(startdate,enddate)
    if startdate > enddate
        return - months(enddate, startdate)
    end
    mmax = (year(enddate)-year(startdate) + 1) * 12     
    for i = 0:mmax   
        if startdate + Month(i) > enddate
            return i-1
        end
    end
end

function nweekdate(n, weekday, year, month)
    date = Date(year, month, 1)
    counter = 0
    while counter < n
        if dayabbr(date) == weekday
            counter += 1
        end
        date += Day(1)
    end
    return date - Day(1)
end