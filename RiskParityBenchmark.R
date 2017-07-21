library(BB)
library(biotools)
library(xts)
library(stringr)
#DATA.ROOT <- 'D:\\Nutstore\\Nutstore\\investment\\src\\R\\history data'
DATA.ROOT <- 'E:/projects/rp/R/my-investment-project/history data'
#OUTPUT.ROOT <- 'D:\\Nutstore\\Nutstore\\investment\\src\\R\\output'
OUTPUT.ROOT <- 'E:/projects/rp/R/my-investment-project/output'
BIG.ASSET.TIME.WINDOW=5 #year
SUB.ASSET.TIME.WINDOW=3 #year
TRADING.DAYS=252
COV.COMP.THRESHOLD=1-0.95  #95%
INIT.PORTFOLIO.MONEY=100000

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
  if (label == 'all') {
    files <- list.files(path = DATA.ROOT, full.names = TRUE, 
             recursive = FALSE, include.dirs = FALSE, pattern = "*.csv")
    names <- list.files(path = DATA.ROOT, full.names = FALSE, 
            recursive = FALSE, include.dirs = FALSE, pattern = "*.csv")
  } else {
    for(tmplab in label) {
      tmplab <- remove.label.level(tmplab)
      files <- paste(DATA.ROOT, '\\', tmplab, '.csv', sep='')
      names <- paste(label, '.csv', sep='')
    }
  }
  
  rts <- NULL
  labs <- NULL
  for (i in 1:length(files)) {
    tmp.zoo <- read.csv.zoo(files[i], format="%Y/%m/%d", tz='GMT')
    tmp.xts <- as.xts(tmp.zoo)
    if (i == 1) {
      rts <- tmp.xts
    }else {
      rts <- cbind(rts, tmp.xts)
    }
    #remove the surfix -- '*.csv'
    labs <- append(labs, substr(names[i], 1, (nchar(names[i])-4)))
  }
  colnames(rts) <- labs
  #interpolation NA data
  rts <- na.trim(rts)
  rts <- na.approx(rts)
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
aa <- c(0.0025,0.0015,-0.0025,0.0015,0.0225,0.0075,-0.0025,0.0075,0.0625)
a3cov <- matrix(aa, nrow=3, ncol=3, byrow = TRUE)
colnames(a3cov) <-  c('aa','bb','cc')
rownames(a3cov) <- c('aa','bb','cc')
# target function:
# first calcuate the 3 assets' risk contribution
# then compute the variance of the 3 risk contribution
# the target should be that variance of the 3 risk contribution equals to ZERO
risk.target <- function (w, cov_mtx) {  # inputs are the 3 assets weights
  risk_bond <- (w[1]*w[1]*cov_mtx[1,1]+w[1]*w[2]*cov_mtx[1,2]+w[1]*w[3]*cov_mtx[1,3])*10000
  risk_stock <- (w[2]*w[2]*cov_mtx[2,2]+w[1]*w[2]*cov_mtx[1,2]+w[2]*w[3]*cov_mtx[2,3])*10000
  risk_comm <- (w[3]*w[3]*cov_mtx[3,3]+w[1]*w[3]*cov_mtx[1,3]+w[2]*w[3]*cov_mtx[2,3])*10000
  var(c(risk_bond, risk_stock, risk_comm))
}
# set boundry: w1+w2+w3 = 1
Amat <- matrix(c(1,1,1,1,0,0,0,1,0,0,0,1), nrow = 4, ncol = 3, byrow = TRUE)
#[,1] [,2] [,3]
#[1,]    1    1    1
#[2,]    1    0    0
#[3,]    0    1    0
#[4,]    0    0    1
b <- c(1,0,0,0)
meq <- 1 # the first line is equal condition: w1+w2+w3 = 1
p0 <- c(1/3, 1/3, 1/3) # the init value for w1,w2,w3

opt_rs <- spg(par=p0, fn=risk.target, cov_mtx= a3cov, project="projectLinear", 
      projectArgs=list(A=Amat, b=b, meq=meq))

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

########## the main procedure ###########
#every time (at least each weekend), the procedure should be run in order to get the following
#information:
#-1. load the previous reference covariance.
#1. current net value
#2. current assets weights
#3. current covariance
#4. compare current covariance with previous reference covariance, if the covariance has been changed
#   then the a Total rebalance will be performed.
#5. TOTAL relalance:
#   use the new covariance to run the optim.
#   get the new assets weight
#   rebalance the assets according to the new weights
#6. normal rebalance:
#   each month end, run a normal rebalance. 

load.curt.pf <- function() {
  path <- paste(OUTPUT.ROOT, sep = '\\', 'portfolio.csv')
  pf <- read.csv(path, row.names = 1, sep = ',')
}

get.fx.lab <- function (lab) {
  currency <- grep('\\.[a-z]{3}$', lab, value=TRUE)
  if(length(currency) == 0) {
    currency='cny'
  }
  tolower(substr(lab, (nchar(lab)-2), nchar(lab)))
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

# not sure this function should be valid any more. 
#get.current.weight<- function (pf) {#?????? should be updated with currency, all should be convert to RMB
  #pf is a matrix, 1st column is volumn, 2nd column is close prices
  # rownames are labs for asset
#  sub.total <- pf[,1]*pf[,2]
#  total <- sum(sub.total)
#  sub.total/total
#} 

output.bchmrk.weights <- function(bchmrk.w.file, w) {
  # rename previous file name by appending the timestamp. 
  # create a new file with new weights.
  write.csv(w, file = bchmrk.w.file)
}

load.currt.base.bchmrk.weights <- function(bchmrk.w.file) {
  # load these value from file, without lever. the content format should be:
  # return a matrix with rownames equal to the first column
  w <- read.csv(bchmrk.w.file, row.names = 1)
}

is.leaf.lable <- function(lab) {
  # search in the weights.csv file to check if current lable has any '_' appended,
  # if so, then it is not a leaf lable.
  # e.g. the input is 'comm_gold.usd'
  weights.file <- load.currt.base.bchmrk.weights(paste(OUTPUT.ROOT, 'weights.csv', sep = '\\'))
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
  weights.file <- load.currt.base.bchmrk.weights(paste(OUTPUT.ROOT, 'weights.csv', sep = '\\'))
  available.weights.labels <- rownames(weights.file)
  if(lab=='root') {
    return(grep('^\\w+$', available.weights.labels, value = TRUE))
  } else {
    # e.g. when input is stock, it should search for the direct descends
    # stock_sp500, stock_hs300
    search.str <- paste('^',escape.weight.lable(lab), '_{1}[A-Za-z0-9\\.]+$', sep = '')
    #print(search.str)
    return(grep(search.str, available.weights.labels, value = TRUE))
  }
  return(NULL)
}

override.sub.weights <- function (w.target, w.source) {
  #both are named vector
  for(lab in names(w.source)) {
    w.target[lab] <- w.source[lab]
  }
  w.target
}

convert.to.target.currency <- function(from, to, ts){
  #from/to -- 3 characters for currency e.g. CNH, USD
  lab <- paste(from, to, sep = '')
  currency <- NULL
  is.reverse <- FALSE
  tryCatch(
    {
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
  # e.g. sub.pf.lab='stock', then should load sp500, hs300 from file.
  # and compare with the corrent sp500, hs300 ts. 
  # ts is uptodate, the 2 data window(new and old) just come from this ts
  # first window uptodate data. 5 years of data. !!!should be return!!!
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
  tryCatch(
    {
      filename <- paste(lab, '.csv', sep = '')
      path <- paste(OUTPUT.ROOT, 'sub_netvalue', filename, sep = '\\')
      sub.pf.value <- read.csv.zoo(filename, format="%Y/%m/%d", tz='GMT')
      sub.pf.value <- as.xts(sub.pf.value)
      path <- paste(OUTPUT.ROOT, 'sub_portfolio', filename, sep = '\\')
      sub.pf.cfg <- read.csv(path, row.names = 1, sep = ',')
    },
    error=function(cond){
      print(paste('sub portfolio', lab, 'does not exist', sep=' '))
      return(NULL)
    },
    warning=function(cond) {
      print(paste('sub portfolio', lab, 'does not exist', sep=' '))
      return(NULL)
    }
  )
  rs <- list(sub.pf.value, sub.pf.cfg)
  names(rs) <- c('value', 'cfg')
  return(rs)
}

init.pf <- function(ts, end.date, hist.window = BIG.ASSET.TIME.WINDOW, period='weeks') {
  rts <- get.pre.n.years.rt(ts, end.date, n = hist.window, period=period)
  cov.mtx <- rts$cov
  rp.weights <- exe.optim(cov.mtx)
  money <- INIT.PORTFOLIO.MONEY * rp.weights
  volumn <- money / ts[end.date]
  std <- sqrt(diag(cov.mtx))
  results <- list(rp.weights, volumn, std, INIT.PORTFOLIO.MONEY)
  names(results) <- c('w','vol','std', 'net')
  return(results)
} 

update.sub.pf.cfg <- function(label, cfg, end.date) {
  #write a dataframe/matrix to disk
  #cfg is a list of a list. outter list contains w, volumn, std
  #inner list contains sub assets.
  #end.date is used to log, when renaming file. the format is 2010-06-08
  w <- cfg$w
  pf.name <- names(w)
  weight <- as.numeric(w)
  volumn <- as.numeric(cfg$vol)
  std <- as.numeric(cfg$std)
  df.cfg <- data.frame(pf.name, weight, volumn, std, row.names = 1)
  # check if file exists, if so rename it
  cfg.file <- paste(OUTPUT.ROOT, 'sub_portfolio', label, sep = '\\')
  cfg.file <- paste(cfg.file, 'csv', sep = '.')
  if(file.exists(cfg.file)) {
    # rename the existing file
    rename.to <- paste(cfg.file, end.date, sep = '.')
    file.rename(cfg.file, rename.to)
  } 
  #write the cfg to a new file
  write.csv(df.cfg, file = cfg.file)
}

update.sub.pf.value <- function(label, ts) {
  #write the ts to disk
  # check if file exists, if so rename it
  ts.file <- paste(OUTPUT.ROOT, '\sub_netvalue', label, sep = '\\')
  ts.file <- paste(cfg.file, 'csv', sep = '.')
  colnames(ts) <- c('value')
  #just wirte the ts value directly to the disk. the ts contains the full values.
  #so every time just overwrite is required. 
  write.zoo(ts, ts.file, sep = ',')
} 

calc.rp.pf.value <- function (current.sub.ts, hist.pf.ts, cfg, end.date) {
  # current.sub.ts -- contains the ts of each sub asset which consist hist.pf.ts.
  # hist.pf.ts -- history net value of current level (a level higher then current.sub.ts) ts.
  # cfg -- current level portfolio config info including: weight, volumn, std
  # end.date -- the end date for construting the current level ts. 
  ## First check the validities of the input data:
  # 1. the last date of hist.pf.ts should be in current.sub.ts, or the current.sub.ts 
  # was not correctly constructed. 
  end.date <- as.POSIXct(end.date, tz='GMT')
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
  if(end.date > ts.days[length(ts.days)]) {
    e <- simpleError(paste('The end date (', end.date, ') exceeds the ts range(', ts.days[length(ts.days)], ')', sep = ' '))
    stop(e)
  }
  # 3. hist.pf.ts end date should not be later then the end.date
  if(start.pos > end.date) {
    e <- simpleError('The end date of the history sub portfolio ts (', start.pos, ') is later than current end date (', end.date, ')')
    stop(e)
  }
  ## loop from hist.pf.ts end to end.date
  #######TODO::::::::::
}

allocate.asset.weight <- function (lab='root', end.date, period='weeks') {
  # bottom up strategy. first run risk parity in sub category,
  # then up to higher level of the category.
  # !!! the hist.prices should be windowed !!!
  # return type: a list contains:
  #   - net value of a sub portfolio e.g. stock which contains china, us etc. 
  #   - weight list of the sub portfolio, the weight need to be updated every time 
  #     the level goest up. e.g. from sp500 to stock -- from 1 to 0.3.
  # end.date format : 2007-01-07
  if (is.leaf.lable(lab)) {
    # return price and weight, obviously weight = 1
    # load prices according to  lab and end.date
    ts <- load.all.prices(remove.label.level(lab))
    w <- c(1)
    names(w) <- c(lab)
    result <- list(ts, w)
    return (result)
  } else {
    convert.to <- get.fx.lab(lab)
    # the current lab descendents
    labs <- get.sub.lables(lab)
    ts <- NULL
    w <- c()
    ### construct the ts matrix from sub portfolios. e.g. stock including china stock, us stock etc.
    for (sublab in labs) {
      convert.from <- get.sub.lables(sublab)
      ############## recursive call
      tmp <- allocate.asset.weight(sublab, end.date)
      #convert to target currency
      if(convert.from != convert.to) {
        tmp$ts <- convert.to.target.currency(convert.from, convert.to, tmp$ts)
      }
      # combine the different assets' ts in the same category, e.g. stock: sp500, hs300
      if (is.null(ts)) {
        ts <- tmp$ts
      } else {
        ts <- cbind(ts, tmp$ts)
      }
      # combine the weights of sub items, and they will be updated when optim is done
      w <- append(w, tmp$w)
    }
    colnames(ts) <- labs
    ### use the sub portfolio ts matrix to construct the current level portfolio net value. 
    sub.pf <- load.sub.pf(lab)
    if(is.null(sub.pf)) {
      # it's the first time run for this sub portfolio, so should initialized this portfolio
      init.param <- init.pf(ts, end.date, hist.window = BIG.ASSET.TIME.WINDOW, period=period)
      update.sub.pf.cfg(lab, init.param, end.date)    #save weight, volumn, std to disk, back up previous cfg if any 
      #as it is the first time run to construct the sub portfolio, e.g. stock sub portfolio
      #it's easy to understand the net value of the sub portfolio will just equal the intial 
      #money invested in the sub portfolio. 
      tmp.date <- as.Date(end.date, tz='GMT')
      sub.pf.ts <- zoo(init.param$net, tmp.date) #it's the initial money invested, defined as a constant.
    } else {
      # cfg and net value already saved on disk
      # 1. load sub portfolio ts (previous created).
      #       sub.pf$value   sub.pf$cfg
      # 2. calculate the net value start from previous end [including rebalance], and save to disk
      # actually, the sub pf ts is full ts, so just need to overwrite the file on the disk
      sub.pf.ts <- calc.rp.pf.value(sub.pf$cfg, sub.pf$value, end.date)
      # 3. check if the cov has been changed, if so, produced the new cfg file 
      #     for the use of next time
      
    }
    update.sub.pf.value(lab, as.xts(sub.pf.ts))  #save sub portfolio net value to disk 
    
    if(cov.changed(ts, end.date)){
      # !!! the changed weight can only be used next time of execution!!!
      # it's understandable that you get the signal that the weight needs to be changed,
      # but you can only changed the weight start from now, not from previous time run. 
      # what should be done this time is that to change the portfolio volumn at the end 
      # of the current time window.
      rts <- get.pre.n.years.rt(ts, end.date)
      # run the excel solver like program to get the weight.
      cov.mtx <- cov(rts)
      rp.weights <- exe.optim(cov.mtx)
    } else {
      rp.weights <- load.sub.pf(lab) #??????????????
    }
    # update the calculated weights for sub category. initially the sub category weights
    # are set to 1
    w <- override.sub.weights (w, rp.weights)
    # construt the risk parity sub portfolio net value. 
    rp.ts <- construct.rp.pf(ts, rp.weights, cov.mtx)
    
    #return to the upper level.
    w.current <- c(1)
    names(w.current) <- c(lab)
    w <- append(w, w.current)
    result <- list(ts, w)
    return(result)
  }
  return(NULL)
}

############## MAIN #####################
#load latest prices, including everyting in the datasource folder
all.prices <- load.all.prices()

# check whether it's the first time run? if so, the initialization will be needed
pf.output <- paste(OUTPUT.ROOT, sep = '\\', 'portfolio.csv')
if (!file.exists(pf.output)) { #first time run
  init.rp.strategy(all.prices)
} else {  #regular weekly run

  #calculate current portfolio net value
  pf <- load.curt.pf()
  latest.prices <- all.prices[length(index(all.prices)),]
  net.value <- get.pf.value(pf, latest.prices)
  # print/log the net value information 
  #
  # need rebalance ?????????
  ## if it's first time run, then 
  # --- normal rebalance
  # --- cov change rebalance
}

#aa2 <- as.vector(aa)
#names(aa2) <-  c('qq','ee','rr')
#aa2['qq']