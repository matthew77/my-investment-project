library(PerformanceAnalytics)

# load my source code here
#IMPORT_PATH <- 
source('RiskParityBenchmark.R')

SortFundsByIR <- function(available.funds.code.list, benchmark.ts, end.date, data.path){
  # return a data.frame with order, fund code, IR score
  rs.IR <- data.frame()
  for (code in available.funds.code.list[,1]) {
    print(paste('processing ... ', code))
    net.value.ts <- load.all.prices(code, root.path = data.path)
    # The fund maybe closed, so there maybe a very large 
    # gap between the final ts date and end.date. 2 weeks can be a threshold.
    if(!IsValidateTs(net.value.ts, end.date)) {
      print (paste('###', code, 'may be out of date'))
      next
    }
    rt.3.years <- GetReturnsPair(net.value.ts, benchmark.ts, end.date, years.span = 3)
    IR.3.years <- InformationRatio(rt.3.years[,1], rt.3.years[,2])
    if (IR.3.years <= 0) {
      print (paste('###', code, 'has negtive IR value::: ', IR.3.years))
      next
    }
    rt.5.years <- GetReturnsPair(net.value.ts, benchmark.ts, end.date, years.span = 5)
    IR.5.years <- InformationRatio(rt.5.years[,1], rt.5.years[,2])
    if (IR.5.years <= 0) {
      print (paste('###', code, 'has negtive IR value::: ', IR.3.years))
      next
    }
    # should also check whether the information ratio is significant
    # it seems that it's not important. since from sqrt(52*3) = 12.5. this is a very big number
    # it may be useful when the IR is near zero.
    if (!IsIRSignificant(IR.3.years, nrow(rt.3.years))) {
      print (paste('###', code, ' does not have significant IR value'))
      next
    }
    if (!IsIRSignificant(IR.5.years, nrow(rt.5.years))) {
      print (paste('###', code, ' does not have significant IR value'))
      next
    }
    # now both 3 & 5 year IR is positive.
    # store the IR.3 value 
    value.IR <- data.frame(list(code=code, IR3=IR.3.years, IR5=IR.5.years))
    rs.IR <- rbind(rs.IR, value.IR)
  }
  # now I have the IR available. then sort the data frame to get the score (position)
  rank.by.IR3 <- order(rs.IR[,'IR3'], decreasing = TRUE)
  rank.by.IR5 <- order(rs.IR[,'IR5'], decreasing = TRUE)
  score.IR <- rank.by.IR3 + rank.by.IR5 # the higher the IR, the lower the score(rank).
  rs.IR <- data.frame(rs.IR, rank.IR3=rank.by.IR3, rank.IR5=rank.by.IR5, score=score.IR)
  score.IR.order <- order(rs.IR[,'score'])
  rs.IR <- rs.IR[score.IR.order,]
  rs.IR
}

IsIRSignificant <- function (value.IR, n, significant.level=0.95){
  # n -- number of observations
  is.significant <- FALSE
  t.value <- sqrt(n) * value.IR
  df <- n - 1 #degree of freedom
  critical.value <- qt(significant.level, dof)
  if (t.value > critical.value) {
    is.significant <- TRUE
  }
  is.significant
}

IsValidateTs <- function(ts, end.date, gap.allowed = 10) {
  # if the the end date of ts is very different of end.date. maybe, the fund has been closed.
  # gap.allowed = 10 days, if exceed that number, the ts should not be valid.
  is.ok <- TRUE
  end.date.obj <- as.POSIXct(end.date, tz='GMT')
  # to get when the ts ends
  hist.days <- index(hist.pf.ts)
  ts.end.pos <- hist.days[length(hist.days)] 
  gap <- end.date.obj - ts.end.pos
  if(gap > gap.allowed) {
    is.ok <- FALSE  
  }
  is.ok
}

GetReturnsPair <- function(ts, benchmark.ts, end.date, years.span=3){
  combined.ts <- cbind(ts, benchmark.ts)
  combined.ts <- na.trim(combined.ts)
  combined.ts <- na.omit(combined.ts)
  window.ts <- n.year.window.by.end.date(combined.ts, end.date, n.years=years.span)
  rt <- get.normal.rt(window.ts)
  rt
}