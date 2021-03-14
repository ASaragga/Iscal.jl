#=
Timetables in Finance

2.1 - Create Timearrays:
       timetable	  - Timearray array with time-stamped rows and variables of different types
       table2timetable - Convert table to timearray
       array2timetable - Convert homogeneous array to array

2.2 - Size, shape and time ranges:
       istimetable    - True for timearray.
       size           - Size of a timearray.
       width          - Number of variables in a timearray.
       height         - Number of rows in a timearray.
       ndims          - Number of dimensions of a timearray.
       numel          - Number of elements in a timearray.
       summary	 - Print summary of table, timearray, or categorical array
       timerange	 - Time range for timearray row subscripting
       overlapsrange  - TRUE for a timearray whose row times overlap the specified time range.
       withinrange    - TRUE for a timearray whose row times are entirely within the specified time range.
       containsrange  - TRUE for a timearray whose row times fully contain the specified time range.

2.3 - Data reorganization:
       addvars       - Insert new variables at a specified location in a table.
       movevars      - Move table variables to a specified location.
       removevars    - Delete the specified table variables.
       renamevars    - Rename variables in timearray.
       sortrows      - Sort rows of a timearray.
       issorted      - TRUE for a sorted timearray.
       stack         - Stack up data from multiple variables into a single variable.
       unstack       - Unstack data from a single variable into multiple variables.

2.4 - Clean Data
       ismissing     - Find elements in a timetable that contain missing values
       rmmissing     - Remove missing entries
       fillmissing	- Fill missing values

2.5 - Manipulate Data in Timearray: 
       retime        - Resample or aggregate data in timearray, and resolve duplicate or irregular times
       rowfun        - Apply function to table or timearray rows
       varfun        - Apply a function to variables in a timearray
       lag           - Lag data in a timearray
       lead          - Lead data in a timearray
       findgroups    - Find groups and return group numbers
       splitapply    - Split data into groups and apply function

2.61 - Combine Timetables and Synchronize Data:
       synchronize	- Synchronize timetables to common time vector, and resample or aggregate data from input timearrays
       horzcat	- Concatenate arrays horizontally
       vertcat	- Concatenate arrays vertically
       intersect     - Find rows common to two timearrays.
       ismember      - Find rows in one timearray that occur in another timearray.
       setdiff       - Find rows that occur in one timearray but not in another.
       setxor        - Find rows that occur in one or the other of two timearrays, but not both.
       unique        - Find unique rows in a timearray.
       union         - Find rows that occur in either of two timearrays.
       join          - Merge two timearrays by matching up rows using key variables.
       innerjoin     - Inner join between two timearrays.
       outerjoin     - Outer join between two timearrays.
=#

# Load local file
function load(localfile::String)
       return readtimearray(localfile)
end

# Load a csv file to a TimeArray via CSV.jl
function load(filename::String, timestamp)
       return TimeArray(CSV.File(filename), timestamp = timestamp)
end

# Save a TimeArray via CSV.jl
function save(filename::String, ta)
       return CSV.write(filename, ta)
end

# Convert TimeArray to DataFrame
function ta2df(ta)
       return DataFrame(ta)
end

# Convert DataFrame to TimeArray. In this case, one needs to point out the column of time index via the timestamp keyword argument.
function df2ta(df, timestamp)
       return TimeArray(df, timestamp = timestamp)
end


# Convert DataFrame to Matrix
# function Matrix(df)

# Convert TimeArray to Matrix
# function values(ta)

              
# Modify a DataFrame, keeping columns of a certain type and discarding columns of all other types. The default type of columns to keep is Union{Real, Missing}.
function concord(df::DataFrames.DataFrame, T = Union{Real, Missing})
       return select(df, findall(col -> eltype(col) <: T, eachcol(df))) # Select only columns of type T
end

