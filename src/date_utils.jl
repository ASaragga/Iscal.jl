# Returns the number of days falling on 29th February between two dates (see conventions NL365 and ActualActualL.  
function num29feb(startdate, enddate)
    if startdate > enddate
        return num29feb(enddate, stardate)
    elseif startdate == enddate 
        return 0
    end 

    y1 = year(startdate)
    y2 = year(enddate)
    dy = y2-y1
    if mod(dy,4) == 0
        if mod(y1,4) == 0
            nleap = 1 + dy/4
        else
            nleap = dy/4
        end
    elseif mod(dy,4) == 1
        if mod(y1,4) == 0 || mod(y1,4) == 3
            nleap = 1 + floor(dy/4)
        else
            nleap = floor(dy/4)
        end
    elseif mod(dy,4) == 2
        if mod(y1,4) == 1
            nleap = floor(dy/4)
        else
            nleap = 1 + floor(dy/4)
        end
    else # case mod(dy,4) == 3
        nleap = 1 + floor(dy/4)
    end
    # The years 1700, 1800, and 1900 were not leap years; neither will 2100, 2200 and 2300. Conversely, the years 1600 and 2000 were leap years as will be 2400 (Wikipedia). 
    c1 = ceil(y1/100)*100
    c2 = floor(y2/100)*100
    notleap = zeros(0)  # notleap is a vector that will contain non-leap years e.g. 1800, 1900, 2100, 2200 
    for x in c1:100:c2
        if mod(x,400) != 0
            append!(notleap,x)
        end
    end
    for x in notleap
        if y1 <= x && x <= y2
            nleap = nleap - 1
        end
    end
    if isleapyear(y1) && startdate > Date(y1,2,29)
        nleap = nleap - 1
    end
    if isleapyear(y2) && enddate < Date(y2,2,29)
        nleap = nleap - 1
    end
    return nleap
end

#= 
A much simpler, but slower function is: 
function num29feb(startdate, enddate)
    if startdate > enddate
        return num29feb(enddate, startdate)
    end
    dr = startdate:Day(1):enddate
    return count(x -> monthday(x) == (2, 29), dr)/365
end
=#

# End-of-month rule. This rule applies only when Maturity is an end-of-month date for a month having 30 or fewer days. 
function eomrule(dt, nmonths, endmonthrule::Bool = true)
    if day(dt) == daysinmonth(Date(year(dt),month(dt),1)) && day(dt) <= 30 && endmonthrule == true 
        newdate = dt + Month(nmonths)
        return Date(year(newdate), month(newdate), daysinmonth(newdate))
    else
        return dt + Month(nmonths)
    end
end
