#= 
Business Days (19/19):
OK  busdate 	- Next or previous business day
OK  datemnth	- Date of day in future or past month
OK  datewrkdy	- Date of future or past workday
OK  days252bus	- Number of business days between dates
OK  days360	    - Days between dates based on 360-day year
OK  days360e	- Days between dates based on 360-day year (European)
OK  days360isda	- Days between dates based on 360-day year (ISDA compliant)
OK  days360psa	- Days between dates based on 360-day year (Public Securities Association (PSA) compliant)
OK  days365	    - Days between dates based on 365-day year
OK  daysact	    - Actual number of days between dates
OK  daysadd	    - Date away from starting date for any day-count basis
OK  daysdif	    - Days between dates for any day-count basis
OK  fbusdate	- First business date of month
OK  isbusday	- True for dates that are business days
OK  lbusdate	- Last business date of month
OK  quarter	    - Returns the quarter of given date
OK  thirdwednesday  - Find third Wednesday of month
OK  wrkdydif	- Number of working days between dates
OK  yearfrac	- Fraction of year between dates
=#

busdate(d::Date, ::Unadjusted, c::BusinessDays.HolidayCalendar) = d
function busdate(dt::Date, ::ModifiedFollowing, c::BusinessDays.HolidayCalendar)
    newDt = busdate(dt, Following(), c)
    if month(newDt) != month(dt)
        return busdate(dt, Preceding(), c);
    end
    return newDt
end
function busdate(dt::Date, ::Following, c::BusinessDays.HolidayCalendar)
    newDt = dt
    while (!isbday(c, newDt))
            newDt = newDt + Day(1)
    end
    return newDt
end
function busdate(dt::Date, ::ModifiedPreceding, c::BusinessDays.HolidayCalendar)
    newDt = busdate(dt, Preceding(), c)
    if month(newDt) != month(dt)
        return busdate(dt, Following(), c);
    end
    return newDt
end
function busdate(dt::Date, ::Preceding, c::BusinessDays.HolidayCalendar)
    newDt = dt
    while (!isbday(c, newDt))
            newDt = newDt - Day(1)
    end
    return newDt
end



function datemnth(d::Date, NumberMonths, dayflag, endmonthrule)
    newdate = d + Month(NumberMonths)
    if dayflag == 0
        if endmonthrule == 1 && daysinmonth(d) <= 30 && day(d) == daysinmonth(d)
            newdate = Date(year(newdate), month(newdate), daysinmonth(newdate))
        end
        return newdate
    elseif dayflag == 1
        return Date(year(newdate),month(newdate), 1)
    elseif dayflag == 2
        return Date(year(newdate), month(newdate), daysinmonth(newdate))
    end
end
datemnth(d, NumberMonths) = datemnth(d::Date, NumberMonths, 0, 0)



function datewrkdy(startdate, nweekdays, c::BusinessDays.HolidayCalendar)
    date = startdate
    nwd = 0
    while nwd < nweekdays
        date +=  Day(1)
        if BusinessDays.isbday(c, date) == true
            nwd += 1
        end
    end
    return date
end



function daysadd(date, numdays, basis)
    newdate = date + Day(numdays)
    if daysdif(date,newdate,basis) == numdays
        return newdate
    else 
        newdate = newdate + Day(numdays - daysdif(date,newdate,basis))
        return newdate
    end
end



function fbusdate(year,month, c::BusinessDays.HolidayCalendar)
    date = Date(year, month, 1)
    while BusinessDays.isbday(c, date) == false
        date += Day(1)
    end
    return date
end



isbusday(date, c::BusinessDays.HolidayCalendar) = BusinessDays.isbday(c, date)

function lbusdate(year,month, c::BusinessDays.HolidayCalendar)
    date = Date(year,month,daysinmonth(Date(year,month,1)))
    while BusinessDays.isbday(c, date) == false
        date -= Day(1)
    end
    return date
end





function quarter(date,month1)
    quarterofyear(date-Month(month1-1))
end
quarter(date) = quarterofyear(date)



function thirdwednesday(month,year)
    d = Date(year, month, 1)
    wed = 0
    while wed < 3
        if dayofweek(d) == 3
            wed +=1
            if wed == 3 
                return d, d + Month(3)
            end
        end 
        d += Day(1)
    end
end



wrkdydif(startdate, enddate, c::BusinessDays.HolidayCalendar) = Dates.value(bdays(c::BusinessDays.HolidayCalendar, startdate, enddate))



function daysdif(startdate, enddate, basis::ActualActualMatlab)
    daysact(startdate,enddate)
end
function daysdif(startdate, enddate, basis::Actual360)
    daysact(startdate,enddate)
end
function daysdif(startdate, enddate, basis::Actual365Fixed)
    daysact(startdate,enddate)
end
function daysdif(startdate, enddate, basis::ActualActualICMA)
    daysact(startdate,enddate)
end
function daysdif(startdate, enddate, basis::Actual360ICMA)
    daysact(startdate,enddate)
end
function daysdif(startdate, enddate, basis::Actual365ICMA)
    daysact(startdate,enddate)
end
function daysdif(startdate, enddate, basis::ActualActualISDA)
    daysact(startdate,enddate)
end
function daysdif(startdate, enddate, basis::Thirty360SIA)
    days360(startdate,enddate)
end
function daysdif(startdate, enddate, basis::Thirty360PSA) 
    days360psa(startdate,enddate)
end
function daysdif(startdate, enddate, basis::Thirty360)
    days360isda(startdate,enddate)
end
function daysdif(startdate, enddate, basis::ThirtyE360)
    days360e(startdate,enddate)
end
function daysdif(startdate, enddate, basis::Thirty360ICMA)
    days360e(startdate,enddate)
end    
function daysdif(startdate, enddate, basis::NL365)
    days365(startdate,enddate)
end
function daysdif(startdate, enddate, basis::Bus252)
    days252bus(startdate,enddate)
end


# helper function
thirty(dy,dm,dd) = (360*dy + 30*dm + dd)

function daysact(startdate,enddate)     # Matlab basis: 0, 2, 3, 8, 9, 10, 12
    return Dates.value(enddate-startdate)
end

function days360(startdate,enddate)     # Matlab basis: 1
    dy = year(enddate)-year(startdate)
    m1 = month(startdate)
    m2 = month(enddate)
    d1 = day(startdate)
    d2 = day(enddate)
    #if dc.eom == true <- We assume that the end of month rule is true
        if startdate == lastdayofmonth(startdate) && m1 == 2
            d1 = 30
            if enddate == lastdayofmonth(enddate) && m2 == 2
                d2 = 30
            end
        end
    #end
    if d1 >= 30
        d1 = 30
        if d2 >= 30
            d2 = 30
        end
    end
    return thirty(dy,m2-m1,d2-d1)
end


function days360psa(startdate,enddate)  # Matlab basis: 4
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
    return thirty(dy,m2-m1,d2-d1)
end

function days360isda(startdate,enddate) # Matlab basis: 5
    dy = year(enddate)-year(startdate)
    dm = month(enddate)-month(startdate)

    d1 = day(startdate)
    d2 = day(enddate)
    if d1 >= 30
        d1 = 30
        if d2 >= 30
            d2 = 30
        end
    end
    return thirty(dy,dm,d2-d1)
end

function days360e(startdate,enddate)    # Matlab basis: 6, 11
    dy = year(enddate)-year(startdate)
    dm = month(enddate)-month(startdate)
    d1 = min(day(startdate),30)
    d2 = min(day(enddate),30)
    return thirty(dy,dm,d2-d1)
end

function days365(startdate,enddate)     # Matlab basis: 7
    if month(enddate) == 2 && day(enddate) == 29
        # Since the startdate is included and the enddate is excluded, when the enddate is 29th February one should't consider this ocorrence in num29feb(startdate,enddate).
        if enddate == startdate
            return 0
        else 
            return (Dates.value(enddate - startdate) - num29feb(startdate, enddate) + 1)
        end
    else
        return (Dates.value(enddate - startdate) - num29feb(startdate, enddate))
    end
end


function days252bus(startdate,enddate,calender)  # Matlab basis: 13
    bdayscount(calender, startdate, enddate)        
end
days252bus(startdate,enddate) = days252bus(startdate,enddate,USNYSE()) # Matlab default is USNYSE()
