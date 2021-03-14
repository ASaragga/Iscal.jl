#= 
Interoperability between Julia and Excel and MATLAB (3/3):
OK  datenum     - Convert date and time to Unix, Excel and MATLAB based serial date numbers
OK  datestr     - Converts serial date numbers based on Unix, Excel and Matlab epochs to date/time
OK  transbasis	- Correpondence between day count basis in Julia versus Excel and MATLAB
=#

function transbasis(n::Int, app::Symbol = :MATLAB)
    if app == :MATLAB
        if n == 0
            return ActualActualMatlab()
        elseif n == 1
            return Thirty360SIA()
        elseif n == 2
            return Actual360()
        elseif n == 3
            return Actual365Fixed()
        elseif n == 4
            return Thirty360PSA()
        elseif n == 5
            return Thirty360()
        elseif n == 6
            return ThirtyE360()
        elseif n == 7
            return NL365()
        elseif n == 8
            return ActualActualICMA() ##
        elseif n == 9
            return Actual360ICMA()
        elseif n == 10
            return Actual365ICMA()
        elseif n == 11
            return ThirtyE360ICMA()
        elseif n == 12
            return ActualActualISDA()
        elseif n == 13
            return Bus252()
        end
    elseif app == :Excel
        if n == 0
            return Thirty360Excel()
        elseif n == 1
            return ActualActualExcel()
        elseif n == 2
            return Actual360()
        elseif n == 3
            return Actual365Fixed()
        elseif n == 4
            return ThirtyE360()
        end
    end
end


#=
datenum convert date and time to Unix, Excel and MATLAB based serial date numbers

    datenum(dt, epochtype)

Input arguments:
    epochtype - Possible values are:
        :unix - The :unix serial date number represents the whole and fractional number of days from the fixed, preset date of January 1, 1970 (midnight UTC/GMT), not counting leap seconds (in ISO 8601: 1970-01-01T00:00:00Z)
        :matlab - The :matlab serial date number represents the whole and fractional number of days from the fixed, preset date of January 0, 0000 in the proleptic ISO calendar
        :excel - Serial date numbers using the default 1900 date system. Excel erroneously believes 1900 was a leap year, so after February 28, 1900, serial numbers consider an extra day. If year is between 0 (zero) and 1899 (inclusive), Excel adds that value to 1900 to calculate the year. For example, DATE(108,1,2) returns January 2, 2008 (1900+108).
=#
function datenum(dt::Union{DateTime, Date}, epochtype::Symbol = :matlab)
    mininday = 86400
    zerounix2matlab = 719529
    zerounix2excel1900 = 25568
    if isa(dt, DateTime) == true
        if epochtype == :unix
            return Dates.datetime2unix(dt)/mininday
        elseif epochtype == :matlab 
            return Dates.datetime2unix(dt)/mininday + zerounix2matlab
        elseif epochtype == :excel
            if dt > DateTime(1900,2,28,23,59,59)
                return Dates.datetime2unix(dt)/mininday + zerounix2excel1900 + 1
            elseif dt >= DateTime(1900,1,1,0,0,0) && dt <= DateTime(1900,2,28,23,59,59)
                return Dates.datetime2unix(dt)/mininday + zerounix2excel1900
            elseif year(dt) >= 0 && year(dt) <= 1899 # If year is between 0 (zero) and 1899 (inclusive), Excel adds that value to 1900 to calculate the year. 
                dt2 = dt + Year(1900)
                return datenum(dt2, :excel)
            end
        end
    else isa(dt, Date) == true
        dt = DateTime(year(dt), month(dt), day(dt),0,0,1)
        if epochtype == :unix
            return floor(Dates.datetime2unix(dt)/mininday)
        elseif epochtype == :matlab
            return floor(Dates.datetime2unix(dt)/mininday) + zerounix2matlab
        elseif epochtype == :excel
            if dt > DateTime(1900,2,28, 24, 0, 0)
                return floor(Dates.datetime2unix(dt)/mininday) + zerounix2excel1900 + 1
            elseif dt >= DateTime(1900,1,1,0,0,0) && dt <= DateTime(1900,2,28,24,0,0)
                return floor(Dates.datetime2unix(dt)/mininday) + zerounix2excel1900
            elseif year(dt) >= 0 && year(dt) <= 1899 # If year is between 0 (zero) and 1899 (inclusive), Excel adds that value to 1900 to calculate the year.
                dt2 = dt + Year(1900)
                return datenum(dt2, :excel)
            end
        end
    end
end


# datestr - Converts serial date numbers based on Unix, Excel and Matlab epoch to date and time
function datestr(serialdate::Real, epochtype)
    mininday = 86400
    zerounix2matlab = 719529
    zerounix2excel1900 = 25568
    if epochtype == :unix
        dttm = Dates.unix2datetime(serialdate * mininday)
        if floor(serialdate) == serialdate
            return Date(year(dttm), month(dttm), day(dttm))
        else
            return dttm
        end
    elseif epochtype == :matlab
        dttm = Dates.unix2datetime((serialdate - zerounix2matlab) * mininday)
        if floor(serialdate) == serialdate
            return Date(year(dttm), month(dttm), day(dttm))
        else
            return dttm
        end
    elseif epochtype == :excel
        if serialdate >= 60
            dttm = Dates.unix2datetime((serialdate - zerounix2excel1900 - 1) * mininday)
        else
            dttm = Dates.unix2datetime((serialdate - zerounix2excel1900) * mininday)
        end
        if floor(serialdate) == serialdate
            return Date(year(dttm), month(dttm), day(dttm))
        else
            return dttm
        end
    end
end    

#=
m2xdate converts serial date numbers from the MATLAB serial date number format to the Excel serial date number format.

    ExcelDateNumber = m2xdate(MATLABDateNumber, Convention)

Input argument: 
    MATLABDateNumber    - Dates in MATLAB serial date number or datetime form
    Convention          - Scalar or an array of flags indicating which Excel serial date number convention should be used when converting from MATLAB serial date numbers; possible values are:
        0 : 1900 date system in which a serial date number of one corresponds to the date 1-Jan-1900 {default}.
        1 : 1904 date system in which a serial date number of zero corresponds to the date 1-Jan-1904.
Output argument: 
        ExcelDateNumber - Serial date numbers in Excel serial date number form.
=#
