ta = stock(:GOOG, "1d", Date(2020,1,1), Date(2020,6,30))

# Multiple series plots. The recipe allows TimeArray objects to be passed as input to plot. The recipe will plot each variable as an individual line, aligning all variables to the same y axis. backend).
display(plot(ta[:Open, :High, :Low, :Close]))

# Candlestick plots. We have seriestype = :candlestick support that requires four columns exist in the input. They are open, high, low and close (case-insensitive).
display(plot(ta, seriestype = :candlestick, xticks = 10))

# Heikin-Ashi plots. The Heikin-Ashi chart is constructed like a regular candlestick chart, except the formula for calculating each bar is different.
display(plot(ta, seriestype = :heikinashi, xticks = 10))