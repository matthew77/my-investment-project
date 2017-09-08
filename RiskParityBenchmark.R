library(BB)
library(biotools)
library(xts)
library(stringr)
#library(futile.logger) #all print() should be replaced with logging.
DATA.ROOT <- 'D:/MyProject/R/my-investment-project/history data'
#DATA.ROOT <- 'E:/projects/rp/R/my-investment-project/history data'
OUTPUT.ROOT <- 'D:/MyProject/R/my-investment-project/output'
#OUTPUT.ROOT <- 'E:/projects/rp/R/my-investment-project/output'
CONFIG.ROOT <- 'D:/MyProject/R/my-investment-project/cfg'
BIG.ASSET.TIME.WINDOW=5 #year
SUB.ASSET.TIME.WINDOW=3 #year
TRADING.DAYS=252
COV.COMP.THRESHOLD=1-0.95  #95%
INIT.PORTFOLIO.MONEY=100000
INIT.INDEX=10000
MIN.WINDOW.TO.START=100 #weeks. if the ts data is less than this threashold, then the data should not be included.
TRANSACTION.COST = TRUE
LOAN.COST = TRUE
#DATA.FILES <- new.env()
#assign('COMM', paste(DATA.ROOT, sep = '\\', 'benchMarkCOMM.csv'))
#assign('SH50', paste(DATA.ROOT, sep = '\\', 'benchMarkSH50.csv'))
#assign('CYB', paste(DATA.ROOT, sep = '\\', 'benchMarkCYB.csv'))

findCovMat <- function(return_matrix) {
  meanv <- apply(return_matrix,2,mean)
  cov_mat <- cov(return_matrix)
  diag_cov_mat <- diag(cov_mat)
  sdevv <- sqrt(diag(cov_mat))
  list(meanv,cov_mat,diag_cov_mat,sdevv)
}

#get.log.rt <- function(prices_matrix) {#Find R: logrets:
#  len <- dim(prices_matrix)[1]
#  D <- dim(prices_matrix)[2]
#  R = matrix(nrow=(len-1),ncol=D)
#  for(i in 1:D) {
#    R[,i] = 100*diff(log(prices_matrix[,i])) ###log rets
#  }
#  R
#}

get.log.rt <- function(prices) {#Find R: logrets:
  mtx <- 100*diff(log(prices))
  mtx <- na.omit(mtx)
}

###
load.all.prices <- function (label = 'all') {
  #load all the prices from the data source folder. 
  files <- NULL
  names <- NULL
  if ('all' %in% label) {
    files <- list.files(path = DATA.ROOT, full.names = TRUE, 
             recursive = FALSE, include.dirs = FALSE, pattern = "*.csv")
    names <- list.files(path = DATA.ROOT, full.names = FALSE, 
            recursive = FALSE, include.dirs = FALSE, pattern = "*.csv")
  } else {
    for(tmplab in label) {
      tmplab <- remove.label.level(tmplab)
      files <- append(files, paste(DATA.ROOT, '/', tmplab, '.csv', sep=''))
      names <- append(names, paste(label, '.csv', sep=''))
    }
  }
  
  ts <- NULL
  labs <- NULL
  for (i in 1:length(files)) {
    tmp.zoo <- read.csv.zoo(files[i], format="%Y/%m/%d", tz='GMT')
    tmp.xts <- as.xts(tmp.zoo)
    if (i == 1) {
      ts <- tmp.xts
    }else {
      ts <- cbind(ts, tmp.xts)
    }
    #remove the surfix -- '*.csv'
    labs <- append(labs, substr(names[i], 1, (nchar(names[i])-4)))
  }
  colnames(ts) <- labs
  #trim and interpolation NA data
  ts <- na.trim(ts)
  ts <- na.approx(ts)
}

remove.label.level <- function(leveled.label) {
  # remove the level. e.g. comm_gold.usd --> gold.usd
  # the purpose for removing level is because, in prices file, no level is presented in the file name.
  search.str <- '_{1}[A-Za-z0-9\\.]+$'
  m <- regexpr(search.str, leveled.label)
  if (m == -1) {
    #no match, so no leveled information
    return (leveled.label)
  }
  lab <- regmatches(leveled.label, m)
  substr(lab, 2, nchar(lab))
}

###
get.pre.n.years.rt <- function(ts, end, n=SUB.ASSET.TIME.WINDOW, period='weeks') {
  # default n=3, 3 years of data
  # period can be:"seconds", "minutes", "hours", "days", "weeks", "months", "quarters", and "years".
  # end= matrix_xts['2006-01-07/2007-01-07']
  tmp.ts <- window.by.end.date(ts, end, time.window = n)
  rt <- get.normal.rt(tmp.ts)
  # get cor(), cov()
  cov <- cov(rt)
  list(rt=rt, cov=cov)
}
# cov2cor(log_rt_cov); if the correlation will be needed.
# how to simulate the excel solver to get the asset weight
# First, the cov matrix must be ready
#aa <- c(0.0025,0.0015,-0.0025,0.0015,0.0225,0.0075,-0.0025,0.0075,0.0625)
#a3cov <- matrix(aa, nrow=3, ncol=3, byrow = TRUE)
#colnames(a3cov) <-  c('aa','bb','cc')
#rownames(a3cov) <- c('aa','bb','cc')
# target function:
# first calcuate the 3 assets' risk contribution
# then compute the variance of the 3 risk contribution
# the target should be that variance of the 3 risk contribution equals to ZERO
#risk.target <- function (w, cov_mtx) {  # inputs are the 3 assets weights
#  risk_bond <- (w[1]*w[1]*cov_mtx[1,1]+w[1]*w[2]*cov_mtx[1,2]+w[1]*w[3]*cov_mtx[1,3])*10000
#  risk_stock <- (w[2]*w[2]*cov_mtx[2,2]+w[1]*w[2]*cov_mtx[1,2]+w[2]*w[3]*cov_mtx[2,3])*10000
#  risk_comm <- (w[3]*w[3]*cov_mtx[3,3]+w[1]*w[3]*cov_mtx[1,3]+w[2]*w[3]*cov_mtx[2,3])*10000
#  var(c(risk_bond, risk_stock, risk_comm))
#}
# set boundry: w1+w2+w3 = 1
#Amat <- matrix(c(1,1,1,1,0,0,0,1,0,0,0,1), nrow = 4, ncol = 3, byrow = TRUE)
#[,1] [,2] [,3]
#[1,]    1    1    1
#[2,]    1    0    0
#[3,]    0    1    0
#[4,]    0    0    1
#b <- c(1,0,0,0)
#meq <- 1 # the first line is equal condition: w1+w2+w3 = 1
#p0 <- c(1/3, 1/3, 1/3) # the init value for w1,w2,w3

#opt_rs <- spg(par=p0, fn=risk.target, cov_mtx= a3cov, project="projectLinear", 
#      projectArgs=list(A=Amat, b=b, meq=meq))

#opt_rs <- spg(par=p0, fn=optim.target, cov.mtx= a3cov, project="projectLinear", 
#              projectArgs=list(A=Amat, b=b, meq=meq))

exe.optim <- function (cov.mtx) {
  num.asset <- nrow(cov.mtx)
  labs <- rownames(cov.mtx)
  # set the equal condition: w1+w2+w3 = 1
  Amat.seq <- rep(1, num.asset)
  b <- 1
  # set the unequal condition, w1,w2,w3 > 0
  for (i in 1:num.asset) {
    for (j in 1:num.asset) {
      if (i == j) {
        Amat.seq <- append(Amat.seq, 1)
      } else {
        Amat.seq <- append(Amat.seq, 0)
      }
    }
    #[,1] [,2] [,3]
    #[1,]    1    1    1
    #[2,]    1    0    0
    #[3,]    0    1    0
    #[4,]    0    0    1    
    b <- append(b, 0)
  }
  Amat <- matrix(Amat.seq, nrow = (num.asset + 1), ncol = num.asset, byrow = TRUE)
  meq <- 1 # the first line is equal condition: w1+w2+w3 = 1
  p0 <-  rep(1/num.asset, num.asset)# the init value for w1,w2,w3...
  optim.rs <- spg(par=p0, fn=optim.target, cov.mtx= cov.mtx, project="projectLinear", 
                projectArgs=list(A=Amat, b=b, meq=meq))
  print(optim.rs)
  weights <- optim.rs$par
  names(weights) <- labs
  return(weights)
}

optim.target <- function (w, cov.mtx) {
  risk.contrib <- c()
  num.asset <- nrow(cov.mtx)
  for(i in 1:num.asset) {
    tmpsum <- 0
    for (j in 1:num.asset) {
      tmpsum <- tmpsum + w[i]*w[j]*cov.mtx[i,j]
    }
    risk.contrib[i] <- tmpsum*10000
  }
  var(risk.contrib)
}

LoadSubPortfolioCfg <- function(label) {
  path <- paste(OUTPUT.ROOT, 'sub_portfolio', label, sep = '/')
  path <- paste(path, 'csv', sep = '.')
  pf <- read.csv(path, row.names = 1, sep = ',')
}

get.fx.lab <- function (lab) {
  currency <- grep('\\.[a-z]{3}$', lab, value=TRUE)
  if(length(currency) == 0) {
    'cny'
  } else {
    tolower(substr(lab, (nchar(lab)-2), nchar(lab)))
  }
}

get.exchange.rate <- function (lab, prices) {
  #last 4 characters stands for foreign currency e.g. .usd
  currency <- tolower(substr(lab, (nchar(lab)-2), nchar(lab)))
  # hard code, since as expect there should be only at most 3~5 currency.
  if (currency == 'usd') {
    return(prices$usdcny)
  } 
  if (currency == 'eur') {
    return(prices$eurcny)
  } 
}

get.pf.value <- function(pf, latest.prices){
  # pf -- volumn, it's a data.frame
  # latest.prices -- should include currency pair
  # return a vector with the following item:
  # rmb_total, sub_usd, sub_rmb, sub_other.fx.....
  labs <- rownames(pf)
  currency_class <- c('rmb_total', 'sub_usd', 'sub_rmb', 'sub_eur')
  total <- c(0,0,0,0)
  names(total) <- currency_class
  for(i in 1:length(labs)) {
    lab <- labs[i]
    vol.num <- pf[lab,]
    price <- as.numeric(latest.prices[, lab])
    which.fx <- get.fx.lab(lab)
    cur_key <- paste('sub', which.fx, sep = '_')
    if (which.fx != 'cny') {  # it's a foreign asset, so need exchange rate
      exrate <- get.exchange.rate(lab, latest.prices)
      exrate <- as.numeric(exrate)  #convert from ts to numeric
      total[cur_key] <- total[cur_key] + price*vol.num  # subtotal for fx assets
      price <- price*exrate    # convert to cny price
    } else { # it's a cny asset
      total['sub_rmb'] <- total['sub_rmb'] + price*vol.num    # sub total cny assets
    }
    total['rmb_total'] <- total['rmb_total'] + price*vol.num    # sub total cny assets
  }
  total
}

get.normal.rt <- function(prices) {
  lag.prices <- lag(prices, 1)
  normal.rt <- prices/lag.prices - 1
  normal.rt <- na.omit(normal.rt)
}

LoadCSVWithLabelAsRowName <- function(bchmrk.w.file) {
  # load these value from file, without lever. the content format should be:
  # return a matrix with rownames equal to the first column
  rec <- NULL
  tryCatch({
      rec <- read.csv(bchmrk.w.file, row.names = 1, stringsAsFactors = FALSE)
    },
    error=function(cond){
      print(cond)
    },
    warning=function(cond){
      print(cond)
    }
  )
  return(rec)
}

is.leaf.lable <- function(lab) {
  # search in the structure.csv file to check if current lable has any '_' appended,
  # if so, then it is not a leaf lable.
  # e.g. the input is 'comm_gold.usd'
  if(lab == 'RPROOT') {
    return(FALSE)
  }
  weights.file <- LoadCSVWithLabelAsRowName(paste(CONFIG.ROOT, 'structure.csv', sep = '/'))
  available.weights.labels <- rownames(weights.file)
  search.str <- paste('^',escape.weight.lable(lab), '_{1}', sep = '')
  #print(search.str)  # log
  is.match = grepl(search.str, available.weights.labels)
  if(length(which(is.match == TRUE)) > 0) {
    # has matched succesfully, so still has sub labels
    # thus return FALSE since it it not a leaf
    return(FALSE)
  }
  return(TRUE)
}

escape.weight.lable <- function(lab) {
  #escpae '\\|' to '\\\\|'
  #str_replace_all(lab, '_', '_')
  # do NOTHING now. originally, '|' is used as the level separator, but, later it turns out
  # the '|' cannot be used as file name!!!
  lab
}

get.sub.lables <- function(lab) {
  # the input tag should contains the hirarachy information 
  # e.g. 'test_test1_test2'
  weights.file <- LoadCSVWithLabelAsRowName(paste(CONFIG.ROOT, 'structure.csv', sep = '/'))
  available.weights.labels <- rownames(weights.file)
  if(lab=='RPROOT') {
    return(grep('^[A-Za-z0-9\\.]+$', available.weights.labels, value = TRUE))
  } else {
    # e.g. when input is stock, it should search for the direct descends
    # stock_sp500, stock_hs300
    search.str <- paste('^',escape.weight.lable(lab), '_{1}[A-Za-z0-9\\.]+$', sep = '')
    #print(search.str)
    return(grep(search.str, available.weights.labels, value = TRUE))
  }
  return(NULL)
}

# should be deleted
#override.sub.weights <- function (w.target, w.source) {
  #both are named vector
#  for(lab in names(w.source)) {
#    w.target[lab] <- w.source[lab]
#  }
#  w.target
#}

convert.to.target.currency <- function(from, to, ts){
  #from/to -- 3 characters for currency e.g. CNH, USD
  lab <- paste(from, to, sep = '')
  currency <- NULL
  is.reverse <- FALSE
  tryCatch({
      currency <- load.all.prices(lab)
    },
    error=function(cond){
      print (paste('currency', lab, 'does not exist', sep=' '))
    },
    warning=function(cond) {
      print (paste('currency', lab, 'does not exist', sep=' '))
    }
  )
  if (is.null(currency)) {
    lab <- paste(to, from, sep = '')
    tryCatch(
      {
        currency <- load.all.prices(lab)
        is.reverse <- TRUE
      },
      error=function(cond){
        print (paste('currency', lab, 'does not exist', sep=' '))
      },
      warning=function(cond) {
        print (paste('currency', lab, 'does not exist', sep=' '))
      }
    )
  }
  if(is.null(currency)) {
    e <- simpleError(paste('currency pair', lab, 'is not available', sep = ' '))
    stop(e)
  }
  #to compare whether ts & currency has same end date, or at least, the currency data
  #should be latest enough to convert the ts
  idx.cur <- index(currency)
  idx.asset <- index(ts)
  if(idx.cur[length(idx.cur)] < idx.asset[length(idx.asset)]) {
    # the currency data is not uptodate!!!
    e <- simpleError(paste('currency pair', lab, 'is not up to date', sep = ' '))
    stop(e)
  }
  if(idx.cur[1] > idx.asset[length(1)]) {
    # the currency data is not uptodate!!!
    e <- simpleError(paste('currency pair', lab, 'does not have enough history data', sep = ' '))
    stop(e)
  }
  # convert to target currency
  if(is.reverse) {
    currency <- 1/currency
  }
  combined.ts <- cbind(ts, currency)
  combined.ts <- na.omit(combined.ts)
  return(combined.ts[,1]*combined.ts[,2])
}

window.by.end.date <- function(ts, end.date, time.window=BIG.ASSET.TIME.WINDOW, period='weeks') {
  time.window.day <- TRADING.DAYS*time.window
  window.ts <- ts
  if(!missing(end.date)){
    window.ts <- ts[paste('/', end.date, sep = '')]
  }
  total <- length(index(window.ts))
  if (total > time.window.day) {
    # subset to n years of window 
    window.ts <- window.ts[(total-time.window.day + 1):total]
  }
  window.ts <- to.period(window.ts, period, OHLC=FALSE)
}

cov.changed <- function(ts, end.date.pre, end.date.current, time.window=BIG.ASSET.TIME.WINDOW) {
  # 
  if(end.date.pre == end.date.current) {
    return(FALSE)
  }
  window.new <- window.by.end.date(ts, end.date.current, time.window = time.window)
  rt.new <- get.normal.rt(window.new)
  rt.new <- as.matrix(rt.new) # boxM does not support xts
  flag.new <- rep(1, nrow(rt.new))
  window.pre <- window.by.end.date(ts, end.date.pre, time.window = time.window)
  rt.pre <- get.normal.rt(window.pre)
  rt.pre <- as.matrix(rt.pre)
  flag.pre <- rep(2, nrow(rt.pre))
  data <- rbind(rt.new, rt.pre)
  group <- c(flag.new, flag.pre)
  fgroup <- factor(group)
  result <- boxM(as.matrix(data), fgroup)
  # when 2 identical matrix were compared, we will get the following results:
  # Chi-Sq (approx.) = 0, df = 15, p-value = 1, so:
  # 'p value is large that we cannot reject the cov is equal'
  if(result$p.value < COV.COMP.THRESHOLD) {
    #so we can reject that the covs are equal
    #which means the cov has been changed
    return (TRUE)
  }
  return(FALSE)
}

load.sub.pf <- function(lab) {
  # the information for sub value including 1)history net value, 2) weights
  # 3) volumn, 4)std, should be already stored on the disk. or it is a init run!
  out <- tryCatch({
      filename <- paste(lab, 'csv', sep = '.')
      path <- paste(OUTPUT.ROOT, 'sub_netvalue', filename, sep = '/')
      sub.pf.value <- read.csv.zoo(path, format="%Y/%m/%d", tz='GMT')
      sub.pf.value <- as.xts(sub.pf.value)
      path <- paste(OUTPUT.ROOT, 'sub_portfolio', filename, sep = '/')
      sub.pf.cfg <- read.csv(path, row.names = 1, sep = ',')
      rs <- list(sub.pf.value, sub.pf.cfg)
      names(rs) <- c('value', 'cfg')
    },
    error=function(cond){
      print(paste('sub portfolio', lab, 'does not exist', sep=' '))
      print(cond)
      return(NULL)
    },
    warning=function(cond) {
      print(paste('sub portfolio', lab, 'does not exist', sep=' '))
      print(cond)
      return(NULL)
    }
  )
  return(out)
}

InitRPPF <- function(label, ts, end.date, hist.window = BIG.ASSET.TIME.WINDOW, period='weeks') {
  rts <- get.pre.n.years.rt(ts, end.date, n = hist.window, period=period)
  #check the ts contains enough data.
  idx <- index(rts$rt)
  pos1 <- idx[1]
  pos2 <- idx[length(idx)]
  span.weeks <- as.numeric(pos2 - pos1, units='weeks')
  if(span.weeks < MIN.WINDOW.TO.START) {
    e <- simpleError(paste('Not enough data to start the process :::', colnames(ts), sep = ' '))
    stop(e)
  }
  cov.mtx <- rts$cov
  weight <- exe.optim(cov.mtx)
  money <- INIT.PORTFOLIO.MONEY * weight
  volumn <- as.numeric(money / ts[pos1])
  std <- sqrt(diag(cov.mtx))
  pf.name <- colnames(rts$rt)
  #w.low, w.high should not be need during init, because the init cfg will be 
  #written to the disk immediately.
  cfg <- data.frame(pf.name, weight, volumn, std, row.names = 1)
  cfg <- SaveRPPFCfg(label, cfg, format(pos1))    #save weight, volumn, std to disk, back up previous cfg if any 
  UpdateCovChangeDate(label, end.date)  #save the init cov change date. 
  # consturct the ts before the run day. 
  first.day.pf <- zoo(INIT.PORTFOLIO.MONEY, pos1) #it's the initial money invested, defined as a constant.
  first.day.pf <- as.xts(first.day.pf)
  init.pf.ts <- CalcuRPTS(label, ts, first.day.pf, cfg, end.date, init.run=TRUE)
  results <- list(cfg, init.pf.ts)
  names(results) <- c('cfg', 'ts')
  return(results)
} 

SaveRPPFCfg <- function(label, cfg, end.date) {
  #write a dataframe/matrix to disk
  #cfg is a list of a list. outter list contains w, volumn, std
  #inner list contains sub assets.
  #end.date is used to log, when renaming file. the format is 2010-06-08
  pf.name <- rownames(cfg)
  weight <- as.numeric(cfg$weight)
  volumn <- as.numeric(cfg$volumn)
  std <- as.numeric(cfg$std)
  # calculate weight high & low. 
  #[standard weight1 +/- one std（return）* w1]/ total weight, 
  #if the above scope is broken, then rebalance
  w.std <- weight*std
  w.low <- weight - w.std
  w.high <- weight + w.std
  df.cfg <- data.frame(pf.name, weight, volumn, std, w.low, w.high, row.names = 1)
  # check if file exists, if so rename it
  cfg.file <- paste(OUTPUT.ROOT, 'sub_portfolio', label, sep = '/')
  cfg.file <- paste(cfg.file, 'csv', sep = '.')
  if(file.exists(cfg.file)) {
    # rename the existing file
    rename.to <- paste(cfg.file, end.date, sep = '.')
    file.rename(cfg.file, rename.to)
  } 
  #write the cfg to a new file
  write.csv(df.cfg, file = cfg.file)
  return(df.cfg)
}

update.sub.pf.value <- function(label, ts) {
  #write the ts to disk
  # check if file exists, if so rename it
  ts.file <- paste(OUTPUT.ROOT, '/sub_netvalue', label, sep = '/')
  ts.file <- paste(ts.file, 'csv', sep = '.')
  colnames(ts) <- c('value')
  #just wirte the ts value directly to the disk. the ts contains the full values.
  #so every time just overwrite is required. 
  write.zoo(ts, ts.file, sep = ',')
} 

need.rebalance <- function(cfg, current.w){
  # cfg is a data.frame
  # current.w is a named vector
  labs <- rownames(cfg)
  for(lab in labs){
    rec <- cfg[lab,]
    w.low <- rec[,'w.low']
    w.high <- rec[,'w.high']
    w <- current.w[lab]
    if(w>w.high || w<w.low){
      return(TRUE)
    }
  }
  return(FALSE)
}

rebalance <- function(cfg, sub.pf.value, current.date.ts) {
  #cfg -- was loaded from file
  #sub.pf.value -- current money value of the sub portfolio
  #current.date.ts -- the market prices of different assets on A specifc date.
  labs <- rownames(cfg)
  for(lab in labs) {
    w <- cfg[lab, 'weight']
    money <- sub.pf.value * w
    vol <- money/as.numeric(current.date.ts[, lab])
    #update the cfg
    print(paste('rebalance on', format(index(current.date.ts)), '--', lab, 'volumn changed by [', vol-cfg[lab, 'volumn'], ']'))
    cfg[lab, 'volumn'] <- vol
  }
  cfg
}

CalcuRPTS <- function (parent.lab, current.sub.ts, hist.pf.ts, cfg, end.date, init.run=FALSE) {
  # current.sub.ts -- contains the ts of each sub asset which consist hist.pf.ts.
  # hist.pf.ts -- history net value of current level (a level higher then current.sub.ts) ts.
  # cfg -- current level portfolio config info including: weight, volumn, std
  # end.date -- the end date for construting the current level ts. 
  ## First check the validities of the input data:
  # 1. the last date of hist.pf.ts should be in current.sub.ts, or the current.sub.ts 
  # was not correctly constructed. 
  end.date.obj <- as.POSIXct(end.date, tz='GMT')
  hist.days <- index(hist.pf.ts)
  start.pos <- hist.days[length(hist.days)] #it's a date
  ts.days <- index(current.sub.ts)
  if (!(start.pos %in% ts.days)) {
    e <- simpleError(paste(hist.days[length(hist.days)], 'is not in local ts file, please update local ts file first and then try again.', sep = ' '))
    stop(e)
  }
  # 2. if end.date is not the last day of current.sub.ts, then an error should be raised. I use this
  # strategy for simplify the process. So each time I run the procedure, the end.date should be set 
  # explicitly so that to make sure the ts data are uptodate.
  if(end.date.obj > ts.days[length(ts.days)]) {
    e <- simpleError(paste('The end date (', end.date, ') exceeds the ts range(', ts.days[length(ts.days)], ')', sep = ' '))
    stop(e)
  }
  # 3. hist.pf.ts end date should not be later then the end.date
  if(start.pos > end.date.obj) {
    e <- simpleError('The end date of the history sub portfolio ts (', start.pos, ') is later than current end date (', end.date, ')')
    stop(e)
  }
  ## loop from hist.pf.ts end to end.date to construct the ts.
  one.day <- as.difftime(1, units = "days")
  start.pos <- start.pos+one.day # should start the process at least one day after the previous run.
  ts.in.range <- current.sub.ts[paste(start.pos, end.date, sep = '/')]
  current.pf.ts = NULL
  labs <- rownames(cfg) # get the assets' names in this sub portfolio
  for(i in 1:nrow(ts.in.range)) {
    p <- ts.in.range[i]
    date.obj <- index(p)  # get the date obj for current date.
    sub.total <- c() #store the value of each asset with names set by labs
    #calculate the total value of each asset in the sub portfolio
    for(lab in labs) {
      #get volumn
      vol <- cfg[lab, 'volumn']
      #get price
      p.asset <- as.numeric(p[,lab])
      value <- vol*p.asset
      sub.total <- append(sub.total, value)
    }
    names(sub.total) <- labs
    sub.sum <- sum(sub.total)
    # Normal Rebalance:::::::::
    # [standard weight1 +/- one std（return）* w1]/ total weight, if the above scope is broken, 
    # then rebalance.
    if(need.rebalance(cfg, sub.total/sub.sum)){
      print(paste(parent.lab, ' --- rebalance happened while processing portfolio :::'))
      # recalculate the volumn
      cfg <- rebalance(cfg, sub.sum, p)
      #update the cfg, and write the update to disk
      cfg <- SaveRPPFCfg(parent.lab, cfg, date.obj) 
    }
    # COV Change Rebalance:::::::::
    # during the init run, the cov will not be changed in the init ts construction.
    if(!init.run) {
      # check if the cov has been changed, if so, produced the new cfg file 
      #     for the use of next time
      pre.date <- GetPreCovChangeDate(parent.lab)
      # for big category such as bond, stock and commodity, use long term cycle to check the change 
      # of cov. for sub category such as sp500 and nasdaq, use short term cycle to check the change
      window <- SUB.ASSET.TIME.WINDOW
      if (parent.lab == 'RPROOT') {
        window=BIG.ASSET.TIME.WINDOW
      }
      if(cov.changed(current.sub.ts, pre.date, format(date.obj), time.window = window)) {
        print(paste(parent.lab, ' --- covariance matrix has been changed since last time'))
        # need a rebalance, since the weights of each asset has been changed.
        rts <- get.pre.n.years.rt(ts, format(date.obj))
        # run the excel solver like program to get the new weights.
        new.weights <- exe.optim(rts$cov)
        #update the new weights in cfg
        for(lab in labs) {
          cfg[lab, 'weight'] <- new.weights[lab]
        }
        
        cfg <- rebalance(cfg, sub.sum, p)
        SaveRPPFCfg(parent.lab, cfg, date.obj)    #save weight, volumn, std to disk, back up previous cfg if any 
        UpdateCovChangeDate(parent.lab, format(date.obj))
      }
    }
    #store the net value of the sub portfolio into a ts
    tmp.zoo <- zoo(sub.sum, index(p))
    current.pf.ts = rbind(current.pf.ts, as.xts(tmp.zoo)) #convert to xts to combine. 
  }
  #combine the history ts and new ts and return. 
  total.pf.ts <- rbind(hist.pf.ts, current.pf.ts)
}

UpdateCovChangeDate <- function(lab, date.str){
  cov.date.file <- paste(OUTPUT.ROOT, 'sub_portfolio', 'cov_change.csv', sep = '/')
  if(!file.exists(cov.date.file)){
    # it is the first record in the file, so just create the file and write the single record.
    rec <- data.frame(date = date.str, row.names = lab, stringsAsFactors = FALSE)
    write.csv(rec, cov.date.file)
  } else {
    #rename file to keep the cov change records.
    cov.date <- read.csv(cov.date.file, row.names = 1, stringsAsFactors = FALSE)
    #bug fix. is.null not correct;;; when only one record in df, the rowname is very wired.better to us subset subset(rec, rownames(rec)=='stock')
    rec.date <- subset(cov.date, colnames(cov.date)==lab)
    if (nrow(rec.date)==0) {
      # this record is new, so just append it into the data.frame.
      rec <- data.frame(date = date.str, row.names = lab, stringsAsFactors = FALSE)
      cov.date <- rbind(cov.date, rec)
      write.csv(cov.date, cov.date.file)
    } else {
      # update the record 
      cov.date[lab, 'date'] <- date.str
      write.csv(cov.date, cov.date.file)
    }
  }
}

GetPreCovChangeDate <- function(lab){
  cov.date.file <- paste(OUTPUT.ROOT, 'sub_portfolio', 'cov_change.csv', sep = '/')
  tryCatch(
    {
      cov.date <- read.csv(cov.date.file, row.names = 1, stringsAsFactors = FALSE)
    },
    error=function(cond){
      stop(cond)
    },
    warning=function(cond) {
      stop(cond)
    }
  )
  # return a date str
  cov.date[lab, 'date']
}

CalcuRPAllocation <- function(lab = 'RPROOT', weight = 1, rec = NULL) {
  if (is.leaf.lable(lab)) {
    # update w.low & w.high before return.
    alloc.cfg <- rec[,c('weight', 'std', 'w.low', 'w.high')]
    std <- alloc.cfg[,'std']
    w <- alloc.cfg[,'weight']
    alloc.cfg[,'w.low'] <- w - w*std
    alloc.cfg[,'w.high'] <- w + w*std
  } else {
    alloc.cfg <- NULL
    # load sub pf file
    sub.pf <- LoadSubPortfolioCfg(lab)
    labs <- rownames(sub.pf)
    for (sublab in labs) {
      # get the sub pf weight and pass into the next lower level.
      rec <- sub.pf[sublab, ]
      relative.weight <- sub.pf[sublab, 'weight']
      abs.weight <- weight * relative.weight
      rec[, 'weight'] <- abs.weight
      abs.cfg <- CalcuRPAllocation(sublab, abs.weight, rec)  #recursive call.
      # rbind the returns into a data.frame and return to next upper level
      alloc.cfg <- rbind(alloc.cfg, abs.cfg)
    }
  }
  return(alloc.cfg)
}

AllocateRPAssetWeight <- function (end.date, lab='RPROOT', period='weeks') {
  # bottom up strategy. first run risk parity in sub category,
  # then up to higher level of the category.
  # !!! the hist.prices should be windowed !!!
  # return type:
  #   - net value of a sub portfolio e.g. stock which contains china, us etc. 
  # end.date format : 2007-01-07
  if (is.leaf.lable(lab)) {
    # load prices according to the lab.
    sub.pf.ts <- load.all.prices(remove.label.level(lab)) 
    print(paste('      getting leaf asset :::', lab))
  } else {
    print(paste('constructing the middle/root layer portfolio for :::', lab))
    convert.to <- get.fx.lab(lab)
    # the current lab descendents
    labs <- get.sub.lables(lab)
    ts <- NULL
    ### construct the ts matrix from sub portfolios. e.g. stock including china stock, us stock etc.
    for (sublab in labs) {
      convert.from <- get.fx.lab(sublab)
      ############## recursive call
      tmp.ts <- AllocateRPAssetWeight(end.date, sublab)
      #convert to target currency
      if(convert.from != convert.to) {
        tmp.ts <- convert.to.target.currency(convert.from, convert.to, tmp.ts)
      }
      # combine the different assets' ts in the same category, e.g. stock: sp500, hs300
      if (is.null(ts)) {
        ts <- tmp.ts
      } else {
        ts <- cbind(ts, tmp.ts)
      }
    }
    colnames(ts) <- labs
    #trim and interpolation NA data
    ts <- na.trim(ts)
    ts <- na.approx(ts)
    ### use the sub portfolio ts matrix to construct the current level portfolio net value. 
    sub.pf <- load.sub.pf(lab)
    if(is.null(sub.pf)) {
      # it's the first time run for this sub portfolio, so should initialized this portfolio
      init.pf <- InitRPPF(lab, ts, end.date, hist.window = BIG.ASSET.TIME.WINDOW, period=period)
      sub.pf.ts <- init.pf$ts
    } else {
      # cfg and net value already saved on disk
      # 1. load sub portfolio ts (previous created).
      #       sub.pf$value   sub.pf$cfg
      # 2. calculate the net value start from previous end [including rebalance], and save to disk
      # actually, the sub pf ts is full ts, so just need to overwrite the file on the disk
      sub.pf.ts <- CalcuRPTS(lab, ts, sub.pf$value, sub.pf$cfg, end.date)
    }
    update.sub.pf.value(lab, sub.pf.ts)  #save sub portfolio net value to disk
    if(lab == 'RPROOT') {
      #If it's the root, then output the standard asset (leaf) allocation information to disk. 
      pf.alloc <- CalcuRPAllocation()
      # compare pf.alloc with the one stored on disk.  use::: all.equal()
      # if changed, then rename the old one and save the new one.
      pre.pf.alloc.file <- paste(OUTPUT.ROOT, 'portfolio.csv', sep = '/')
      pre.pf.alloc <- LoadCSVWithLabelAsRowName(pre.pf.alloc.file)
      if(!identical(pre.pf.alloc, pf.alloc)) {
        #rename the previous created csv and save a new CSV.
        if(file.exists(pre.pf.alloc.file)) {
          # rename the existing file
          rename.to <- paste(pre.pf.alloc.file, end.date, sep = '.')
          file.rename(pre.pf.alloc.file, rename.to)
        } 
        #write the cfg to a new file
        write.csv(pf.alloc, file = pre.pf.alloc.file)
      }
    }
  }
  return(sub.pf.ts)
}

GetRPAllocForIndex <- function(end.date, lever) {
  AllocateRPAssetWeight(end.date) #calculate the reference weights allocation
  path <- paste(OUTPUT.ROOT, 'portfolio.csv', sep = '/')
  w <- read.csv(path, row.names = 1, sep = ',')
  w[,c('weight', 'w.low', 'w.high')] <- w[,c('weight', 'w.low', 'w.high')] * lever
  # --- No need to return unleveled. every assets should be leveled, because if the assets
  # are in different level then they should be treated as different assets [All Weather Strategy]
  #unleveled.rowname <- GetUnleveledRowName(rownames(w))   
  #w <- data.frame(w, unleveled.lab=unleveled.rowname, stringsAsFactors = FALSE)
  w
}

GetUnleveledRowName <- function(leveled.rn) {
  unleveled.rn <- sapply(leveled.rn, remove.label.level)
  #return a named vector: key is leveled, value is unleveled. that's very cool!!!
  unleveled.rn
}

#TODO: constant lever (2X) benchmark coding. 
# - 10000rmb as equity. initial market value should be 20000. and the init benchmark index will be 10000
# - what is the rule to keep the 2X constant lever? within 10%, so within 210%~190%? 
# - I will suggest weekly run (using weekly data)
# !!!How to put the momentum factor into the process??? A obvious question for momentum is that:
# 1 million bond and 1 millon stock has very different risk. you cannot simply short the 1 million 
# in bond and then simple long 1 million in stock. that will put too much in stock. so then answer may be 
# run the exe.optim with specified risk allocation. !!!???

CalcuPRIndex <- function(end.date, lever=2){
  #load init data and parameters
  end.date.obj <- as.POSIXct(end.date, tz='GMT')
  index.cfg.file <- paste(CONFIG.ROOT, 'index_cfg.csv', sep = '/')
  index.cfg <- read.csv(index.cfg.file, row.names = 1)
  # get previous risk parity history ts from file. If the file does not exist, 
  # the it is the first time the RP Index is created.
  path.index <- paste(OUTPUT.ROOT, 'rp_index', sep = '/')
  rp.index.file <- paste(path.index, lever, sep = '/')
  rp.index.file <- paste(rp.index.file, 'X.csv', sep = '')  #..../2X.csv
  vol.file <- paste(path.index, 'index_alloc.csv', sep = '/')
  ts.all <- load.all.prices()
  if(!file.exists(rp.index.file)) {
    # init the index. end.date is the first day the index is created.
    print('First time to create Risk Parity index...')
    index.w <- GetRPAllocForIndex(end.date, lever)  # get the weights according to lever, momentum etc...
    unleveled.labs <- GetUnleveledRowName(rownames(index.w))
    ts <- ts.all[, unleveled.labs]
    equity <- index.cfg['init.equity',]
    loan <- equity*lever - equity   # loan for 2X lever
    price <- ts[end.date,]
    market.value <- equity * index.w[,'weight']
    if(TRANSACTION.COST) {
      #buy fee.
      buy.fee <- market.value * index.cfg['buy.commission',]
      equity <- equity - buy.fee
      
    }
    if(LOAN.COST) {
      #TODO:::: minus the loan interest
    }
    volumn <- as.numeric(market.value / price)
    #save the volumn info to disk
    alloc.vol.info <- data.frame(index.w, volumn=volumn)
    write.csv(alloc.vol.info, vol.file)
    #save the ts -- index(equity) value, market value, loan.
    equity.ts <- zoo(equity, end.date.obj)
    equity.ts <- as.xts(equity.ts)
    market.value.ts <- zoo(market.value, end.date.obj)
    market.value.ts <- as.xts(market.value.ts)
    loan.ts <- zoo(loan, end.date.obj)
    loan.ts <- as.xts(loan.ts)
    index.ts <- cbind(equity.ts, market.value.ts, loan.ts)
    colnames(index.ts) <- c('equity.value', 'market.value', 'loan')
    write.zoo(init.ts, rp.index.file)
  } else {
    # load previous created index value ts.
    index.ts <- read.csv.zoo(rp.index.file, tz='GMT')
    index.ts <- as.xts(index.ts)
    # get the intervals (days or weeks) between the last date in the index file and end.date
    hist.days <- index(index.ts)
    start.pos <- hist.days[length(hist.days)] + as.difftime(1, units = "days")
    ts.in.range <- ts.all[paste(start.pos, end.date, sep = '/')]
    index.ts.to.append <- NULL
  # loop through the interval. on each day/week:
    for(i in 1:nrows(ts.in.range)){
      is.ref.w.changed <- FALSE
      is.rebalanced <- FALSE
      # 1. get the uptodate weight allocation
      print(paste('calculating index on :::::::::::', date.obj))
      index.w <- GetRPAllocForIndex(format(date.obj), lever)  # get the weights according to lever, momentum etc...
      unleveled.labs <- GetUnleveledRowName(rownames(index.w))  #here the labs may contains some new assets, e.g. some assets are removed and others are added.
      p <- ts.in.range[i][, unleveled.labs] # filter out those unused prices.
      date.obj <- index(p)
      # 2. rebalance according to the lever(2X), detailed rebalance information will be recorded.
      #     the total portfolio should be around 200% +/- 10% (it's redundant).
      # load volumn information
      alloc.vol.info <- read.csv(vol.file, row.names = 1)
      unleveled.labs.pre <- GetUnleveledRowName(rownames(alloc.vol.info))   # the labs from last time.
      leveled.labs.pre <- names(unleveled.labs.pre)
      p.with.labs.pre <- ts.in.range[i][, unleveled.labs.pre]               # asset prices before adjustment (if any)
      # check if the weight has been changed this time. if so the weight/volumn file needs to be updated.
      if(!identical(index.w[, 'weight'], alloc.vol.info[, 'weight'])) {
        print('-- reference weights have been changed compared with last time!')
        is.ref.w.changed <- TRUE
      } 
      # -- get current market value 
      market.value.before.rb <- as.numeric(p.with.labs.pre * alloc.vol.info[, 'volumn'])
      names(market.value.before.rb) <- names(unleveled.labs.pre)  # need leveled asset name for each asset's weight
      market.value.before.rb.sum <- sum(market.value.before.rb)
      # -- get current equity value
      pre.index <- index.ts[nrow(index.ts)]
      pre.equity <- as.numeric(pre.index[,'equity.value'])
      pre.market.value.sum <- as.numeric(pre.index[,'market.value'])
      pre.loan <- as.numeric(pre.index[,'loan'])
      if(LOAN.COST) {
        #TODO:::: minus the loan interest
        #charge the interest from last time to current time, with the loan number as the principle.
      }
      current.equity <- pre.equity + (current.market.value.sum - pre.market.value.sum)
      # -- rebalance, calculating current weight.
      w.before.rebalance <- current.market.value/current.equity
      print(paste('current allocation (weights) is ::::: '), w.before.rebalance)
      # -- check with asset has exceeded the upper/lower limit.
      loan.rebalance <- 0
      commission <- 0
      for(lab in names(unleveled.labs)) {   # names(lab) returns the leveled lab.
        rec <- index.w[lab,] 
        ref.w.high <- rec$w.high
        ref.w.low <- rec$w.low
        ref.w <- rec$weight
        asset.price <- as.numeric(p[, unleveled.labs[lab]])
        money.for.gap <- 0
        # here the lab is uptodate, may be the lab (asset) was not included last time run (because 
        # the risk parity assets were adjusted)
        asset.w <- as.numeric(w.before.rebalance[lab])
        if(length(asset.w)==0) {
          #this asset is newly added, so no rebalance need for this asset. just buy in directly
          money.for.gap <-  ref.w * current.equity
          volumn.change <- money.for.gap/asset.price
          print(paste('----', lab, 'is newly added in this time! weight gap :::', ref.w, '| money gap :::', money.for.gap, 
                      '| volumn change :::', volumn.change))
          loan.rebalance <- loan.rebalance + money.for.gap
          #update volumn info
          alloc.vol.info <- rbind(alloc.vol.info, data.frame(rec, volumn=volumn.change))
          # is.ref.w.changed <- TRUE    ---MUST be ref weight changed!!!
        } else {
          # this asset is in the portfolio last time, rebalance maybe needed.
          if(asset.w > ref.w.high || asset.w < ref.w.low) {
            # exceed the upper/lower bond, rebalance required.
            gap <- ref.w - asset.w  # POSITIVE = buy asset required; NEGATIVE = sell asset required.
            money.for.gap <- gap * current.equity  # POSITIVE = borrow loan; NEGATIVE=return loan
            # calculate the volumn for the rebalance adjustment.
            volumn.change <- money.for.gap/asset.price
            print(paste('----', lab, 'exceeds the allocation limit! weight gap :::', gap, '| money gap :::', money.for.gap, 
                        '| volumn change :::', volumn.change))
            # update the loan info. if gap is positive, then borrow loan,
            # if gap is negative, then return loan.
            loan.rebalance <- loan.rebalance + money.for.gap
            # update volumn info
            alloc.vol.info[lab, 'volumn'] <- alloc.vol.info[lab, 'volumn'] + volumn.change
            is.rebalanced <- TRUE
          } else if(is.ref.w.changed){
            # why need this check? because there may be the ref weight has been changed, but 
            # the current asset weight can still be in the new range. in this case, the ref
            # weight should be updated into alloc.vol.info and save onto disk
            alloc.vol.info[lab, c('weight', 'std', 'w.low', 'w.high')] <- 
              rec[, c('weight', 'std', 'w.low', 'w.high')]
          }
          # at last, the remaining labs in labs.pre, are those assets that have been removed from 
          # Risk Parity portfolio. so, those assets should be sold. 
          leveled.labs.pre <- leveled.labs.pre[leveled.labs.pre!=lab]
        }
        if(TRANSACTION.COST && money.for.gap != 0) {
          if(money.for.gap > 0) {
            # buy commission
            commission <- commission + money.for.gap * index.cfg['buy.commission',]
          } else {
            # sell commission
            commission <- commission + abs(money.for.gap) * index.cfg['sell.commission',]
          }
        }
      }
      # sell the remaining assets in labs.pre, since those assets have been removed from portfolio.
      for(lab in leveled.labs.pre) {
        money.for.gap <- as.numeric(current.market.value[lab])
        loan.rebalance <- loan.rebalance - money.for.gap  # return money
        if(TRANSACTION.COST) {
          # sell commission
          commission <- commission + abs(money.for.gap) * index.cfg['sell.commission',]
          print(paste('sell all asset:::', lab, 'since it has been removed from Risk Parity Portfolio'))
        }
        # remove the asset in alloc.vol.info
        alloc.vol.info <- subset(alloc.vol.info, rownames(alloc.vol.info) !=lab )
      }
      # update loan number
      loan.after.rb <- pre.loan + loan.rebalance
      if(TRANSACTION.COST) {
        current.equity <- current.equity - commission
        print(paste('total transaction cost on', date.obj, 'is:::::', commission))
      }
      # 3. save the updated volumn file.
      if(is.ref.w.changed || is.rebalanced) {
        # save the volumn & allocation info onto disk since it has been changed since last time.
        rename.to <- paste(vol.file, end.date, sep = '.')
        file.rename(vol.file, rename.to)
        write.csv(alloc.vol.info, vol.file)
      }
      # calculate the market after rebalance. reorder the sequence of the assets to match
      # the asset price ts (the column sequence.)
      volumn.after.rb <- alloc.vol.info[names(unleveled.labs), 'volumn']
      market.value.after.rb <- as.numeric(p * volumn.after.rb)
      market.value.after.rb.sum <- sum(market.value.after.rb)
      print(paste('The lever now is :::', market.value.after.rb.sum/current.equity))
      # 4. assemble the ts.
      equity.ts <- zoo(current.equity, date.obj)
      equity.ts <- as.xts(equity.ts)
      market.value.ts <- zoo(market.value.after.rb.sum, date.obj)
      market.value.ts <- as.xts(market.value.ts)
      loan.ts <- zoo(loan.after.rb, date.obj)
      loan.ts <- as.xts(loan.ts)
      tmp.index.ts <- cbind(equity.ts, market.value.ts, loan.ts)
      colnames(index.ts) <- c('equity.value', 'market.value', 'loan')
      index.ts.to.append <- rbind(index.ts.to.append, tmp.index.ts)
    }# end for
    # write the results onto the file
    write.zoo(index.ts.to.append, rp.index.file, append = TRUE, col.names = FALSE)
  }# end processing.
}

############## MAIN #####################

############## TEST #####################
#TODO: testing:
# 1. very simple senario test. Only 2 asset (CYB, SH50) make sure the whole process is OK.
  #end.date <- '2016-06-01'
  #rp.ts <- AllocateRPAssetWeight(end.date)
  #DONE!
# 2. 3 asset with same currency
#end.date <- '2016-06-01'
#rp.ts <- AllocateRPAssetWeight(end.date)
#DONE!
# 3. 3 asset with different currency
#end.date <- '2016-06-01'
#rp.ts <- AllocateRPAssetWeight(end.date)
#DONE!
# 4. 3 level tree e.g. stock contains us stock (sp500, nasdaq), china stock(CYB, SH50).
#end.date <- '2016-06-01'
#rp.ts <- AllocateRPAssetWeight(end.date)