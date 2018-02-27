library(PerformanceAnalytics)

# load my source code here
#IMPORT_PATH <- 
source('RiskParityBenchmark.R')

SortFundsByIR <- function(available.funds.code.list, benchmark.ts, end.date, data.path){
  # return a data.frame with order, fund code, IR score
  for (code in available.funds.code.list[,1]) {
    #net.value.file <- paste(data.path, paste(code, 'csv', sep = '.'), sep = '/')
    net.value.ts <- load.all.prices(code, root.path = data.path)
    rt.3.years <- GetReturnsPair(net.value.ts, benchmark.ts, end.date, years.span = 3)
    IR.3.years <- InformationRatio(rt.3.years[,1], rt.3.years[,2])
    if (IR.3.years <= 0) {
      next
    }
    rt.5.years <- GetReturnsPair(net.value.ts, benchmark.ts, end.date, years.span = 5)
    IR.5.years <- InformationRatio(rt.5.years[,1], rt.5.years[,2])
    if (IR.5.years <= 0) {
      next
    }
    # now both 3 & 5 year IR is positive.
    # ??? 2 list??? how to sort ???
  }
}

GetReturnsPair <- function(ts, benchmark.ts, end.date, years.span=3){
  combined.ts <- cbind(ts, benchmark.ts)
  combined.ts <- na.trim(combined.ts)
  combined.ts <- na.omit(combined.ts)
  window.ts <- n.year.window.by.end.date(combined.ts, end.date, n.years=years.span)
  rt <- get.normal.rt(window.ts)
  rt
}