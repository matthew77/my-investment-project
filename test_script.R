#test script
rp.ts <- read.csv.zoo('D:\\MyProject\\R\\my-investment-project\\output\\sub_netvalue\\RPROOT.csv', format='%Y-%m-%d', tz='GMT')
stock.ts <- read.csv.zoo('D:\\MyProject\\R\\my-investment-project\\output\\sub_netvalue\\stock.csv', format='%Y-%m-%d', tz='GMT')
comm.ts <-  read.csv.zoo('D:\\MyProject\\R\\my-investment-project\\output\\sub_netvalue\\comm.usd.csv', format='%Y-%m-%d', tz='GMT')
bond.ts <-  read.csv.zoo('D:\\MyProject\\R\\my-investment-project\\output\\sub_netvalue\\bond.usd.csv', format='%Y-%m-%d', tz='GMT')