## README ####
### V 2.2.1
##############

#!/usr/bin/env Rscript
# ChiDist_Grouped.R
# Usage: Rscript ChiDist_Grouped.R <numCells> <data.tsv> <data_filtered.tsv> <out.tsv> <params.RData>

suppressWarnings(library(stringr))
suppressWarnings(library(splines))

# helper for ANY of these “[x, y, z]” columns
safe_split <- function(x) {
  txt <- gsub("^\\[|\\]$", "", as.character(x))
  if (is.na(txt) || txt == "") return(numeric(0))
  as.numeric(unlist(strsplit(txt, ",\\s*")))
}

# ---- 1. parse command-line args ----
args      <- commandArgs(trailingOnly=TRUE)
numcell   <- as.integer(args[1])
data      <- read.table(args[2], sep="\t", stringsAsFactors=FALSE)
data_filtered <- read.table(args[3], sep="\t", stringsAsFactors=FALSE)
out_file  <- args[4]
param_file<- if (length(args)>=5) args[5] else file.path(tempdir(),"ChiDist_pars.RData")

## 设分布都是ZINB分布，先用极大似然估计参数
CalcZINB = function(par)
{
  sum = 0
  for (i in 1:dim(data_filtered)[1])
  {
    s <- as.character(data_filtered$V4[i])
    cc_str <- gsub("^\\[|\\]$", "", s)
    if (is.na(cc_str) || cc_str == "") {
      cc <- numeric(0)
    } else {
      # now we know cc_str is a non-empty "[…]" content string
      pieces <- strsplit(cc_str, ",\\s*")[[1]]
      cc     <- as.numeric(pieces)
    }
    expr1 <- safe_split(data_filtered$V15[i])
    expr2 <- safe_split(data_filtered$V16[i])
    if (sum(is.na(expr1))>0 | sum(is.na(expr2))>0 | sum(is.na(cc))>0){
      next()
    }
    aveexpr1 = mean(expr1)
    aveexpr2 = mean(expr2)
    gc = (data_filtered$V13[i] + data_filtered$V14[i]) / 2
    gcpre = log(gcprelist[round(gc * 1000) + 1])
    lambda = exp(par[5])
    linpart = par[6] + par[7] * gcpre + par[8] * aveexpr1 + par[9] * aveexpr2
    zeroprob = min(exp(linpart) / (1 + exp(linpart)), 0.9999999999999)
    for (j in 1:length(cc)) {
      mu = par[1] + par[2] * gcpre + exp(par[3]) * expr1[j] + exp(par[4]) * expr2[j]
      mu = exp(mu)
      sum = sum + dnbinom(as.numeric(cc[j]),
                          1 / lambda,
                          prob = 1 / (1 + mu * lambda),
                          log = TRUE) + log(1 - zeroprob)
      sum = sum - log(1 - dnbinom(0, 1 / lambda, prob = 1 / (1 + mu * lambda))) # 因为是truncated NB
    }
    sum = sum + (numcell - length(cc)) * log(zeroprob)
    sum = sum - log(1 - zeroprob ^ numcell - numcell * (1 - zeroprob) *
                      zeroprob ^ (numcell - 1)) #条件在至少两个细胞支持
    if (is.nan(sum)){
      break
    }
  }
  - sum
}

CalcPValueZINB_Fusion_2.0.0_3 = function(par)
{
  Pv = rep(0, dim(data)[1])
  for (i in 1:dim(data)[1]){
    s <- as.character(data$V4[i])
    cc_str <- gsub("^\\[|\\]$", "", s)
    if (is.na(cc_str) || cc_str == "") {
      cc <- numeric(0)
    } else {
      # now we know cc_str is a non-empty "[…]" content string
      pieces <- strsplit(cc_str, ",\\s*")[[1]]
      cc     <- as.numeric(pieces)
    }
    if (sum(cc) <= 3 | data$V3[i] <= 1){
      Pv[i] = 1
      next()
    }
    expr1 <- safe_split(data$V15[i])
    expr2 <- safe_split(data$V16[i])
    if (sum(is.na(expr1))>0 | sum(is.na(expr2))>0 | sum(is.na(cc))>0){
      Pv[i] = 1
      next()
    }
    aveexpr1 = mean(expr1) / 2
    aveexpr2 = mean(expr2) / 2
    gc = (data$V13[i] + data$V14[i]) / 2
    gcpre = log(gcprelist[round(gc * 1000) + 1])
    numcellsup = length(cc)
    linpart = par[6] + par[7] * gcpre + par[8] * aveexpr1 + par[9] * aveexpr2
    zeroprob = min(exp(linpart) / (1 + exp(linpart)), 0.9995)
    betamu = 1 - zeroprob
    betavar = betamu * (1 - betamu) / numcell
    myalpha = (betamu * (1 - betamu) / betavar - 1) * betamu
    mybeta = (1 / betamu - 1) * myalpha
    for (j in 1:length(cc)) {
      cc[j] = min(15, cc[j])
    }
    totalread = sum(cc)
    mu = par[1] + par[2] * gcpre + exp(par[3]) * 2 * aveexpr1 + exp(par[4]) * 2 * aveexpr2
    truemu = exp(mu)
    mu = max(truemu, 0.1)
    lambda = exp(par[5])
    truevar = truemu + truemu ^ 2 * lambda
    lambda = max(0.001, (truevar - mu) / (mu ^ 2))
    simuset = rnbinom(10000, 1 / lambda, prob = 1 / (1 + mu * lambda))
    simuset = simuset[simuset > 0]
    whilecount = 0
    while (length(simuset) < 2000){
      newsimuset = rnbinom(10000, 1 / lambda, prob = 1 / (1 + mu * lambda))
      newsimuset = newsimuset[newsimuset>0]
      simuset = c(simuset, newsimuset)
      whilecount = whilecount + 1
      if (whilecount > 100){
        simuset = c(simuset, 1)
      }
    }
    simuset[simuset >= 15] = 15
    totalreadsample = c()
    for (j in 1:100){
      nonzeronum = rbinom(1, numcell, 1-zeroprob)
      usedreadnum = sample(simuset, nonzeronum, replace = T)
      totalreadsample = c(totalreadsample, sum(usedreadnum))
    }
    totalreadmean = mean(totalreadsample)
    totalreadsd = sd(totalreadsample)
    thispvalue = pnorm(totalread, totalreadmean, totalreadsd, lower.tail = F)
    Pv[i] = thispvalue
  }
  #Pv = p.adjust(Pv, method = 'fdr')
  Pv
}

### filter bad gc
neiflag = c()
for (i in 1:dim(data)[1]) {
  if (max(data$V14[i], data$V13[i]) <= 0.05||min(data$V14[i], data$V13[i]) >= 0.95) {
    neiflag[i] = 1
  } else{
    neiflag[i] = 0
  }
}
data <- data[data$neiflag == 0,]
if (nrow(data)>0) {
  row.names(data) <- seq_len(nrow(data))
} else {
  row.names(data) <- integer(0)
}

# fit a spline for gc
gc = round((data$V14 + data$V13) / 2, digits = 3)
gclist = sort(unique(gc))
gclist = round(gclist, digits = 3)
gclist = unique(gclist)
chimexpr = gclist
for (i in 1:length(gclist)) {
  usedata = data[gc == gclist[i],]
  exprtotal = 0
  count = 0
  for (j in 1:dim(usedata)[1]) {
    s <- as.character(usedata$V4[j])
    cc_str <- gsub("^\\[|\\]$", "", s)
    if (is.na(cc_str) || cc_str == "") {
      cc <- numeric(0)
    } else {
      # now we know cc_str is a non-empty "[…]" content string
      pieces <- strsplit(cc_str, ",\\s*")[[1]]
      cc     <- as.numeric(pieces)
    }
    if (length(cc) >= 1) {
      exprtotal = exprtotal + sum(cc)
      count = count + length(cc)
    }
  }
  chimexpr[i] = exprtotal / count
}
lmres = lm(chimexpr ~ bs(gclist, df = 5, intercept = T))
basepre = predict(lmres, data.frame(gclist = gclist))
diff = abs(chimexpr - basepre)
thred = mean(diff) + 1.5 * sqrt(var(diff))
for (i in 1:length(gclist)) {
  pre = basepre[i]
  if (chimexpr[i] - pre > thred) {
    chimexpr[i] = pre
  }
}
lmres = lm(chimexpr ~ bs(gclist, df = 5, intercept = T))
gcprelist = predict(lmres, data.frame(gclist = seq(0, 1, length.out = 1001)))
gcprelist[gcprelist <= 0] = 0.01

print('Start Reading ChiDist File!')
neiflag = rep(0, dim(data)[1])
totalcount = rep(0, dim(data)[1])
datacc = rep('', dim(data)[1])
dataexpmin = rep('', dim(data)[1])
dataexpmax = rep('', dim(data)[1])
datagcmin = rep(0, dim(data)[1])
datagcmax = rep(0, dim(data)[1])
for (i in 1:dim(data)[1]) {
  s <- as.character(data$V4[i])
  cc_str <- gsub("^\\[|\\]$", "", s)
  if (is.na(cc_str) || cc_str == "") {
    cc <- numeric(0)
  } else {
    # now we know cc_str is a non-empty "[…]" content string
    pieces <- strsplit(cc_str, ",\\s*")[[1]]
    cc     <- as.numeric(pieces)
  }
  totalcount[i] = sum(cc)
  smallindex = which(cc <= 2)
  bigindex = which(cc > 2)
  if (max(cc) < 2 | length(cc) < 2) {
    neiflag[i] = 1
  }
  count1 = sum(cc<=2)
  usesmallcellnum = floor(data$V5[i] * 1.1 + 2)
  usesmallcell = c()
  if (count1 <= usesmallcellnum){
    usesmallcell = smallindex
  }else{
    usesmallcell = sample(smallindex, usesmallcellnum)
  }
  cc = cc[union(bigindex, usesmallcell)]
  if (neiflag[i] == 1) {
    next()
  }
  expmin = c()
  expmax = c()
  gene1expr <- safe_split(data$V15[i])
  gene2expr <- safe_split(data$V16[i])
  # if either is empty or any value below threshold, skip
  if (length(gene1expr)==0 ||
      length(gene2expr)==0 ||
      min(c(gene1expr, gene2expr)) < 0.0001) {
    neiflag[i] <- 1
    next()
  }
  datagcmin[i] = min(data$V14[i], data$V13[i])
  datagcmax[i] = max(data$V14[i], data$V13[i])
  for (j in 1:length(gene1expr)) {
    expmin = c(expmin, log(1 + min(
      as.numeric(gene1expr[j]), as.numeric(gene2expr[j])
    ) / 100))
    expmax = c(expmax, log(1 + max(
      as.numeric(gene1expr[j]), as.numeric(gene2expr[j])
    ) / 100))
  }
  expmin = expmin[union(bigindex, usesmallcell)]
  expmax = expmax[union(bigindex, usesmallcell)]
  datacc[i] = paste0('[', str_c(as.character(cc), collapse = ', '), ']')
  dataexpmin[i] = paste0('[', str_c(as.character(expmin), collapse = ', '), ']')
  dataexpmax[i] = paste0('[', str_c(as.character(expmax), collapse = ', '), ']')
}
data$V4 = datacc
data$V15 = dataexpmin
data$V16 = dataexpmax
data$V13 = datagcmin
data$V14 = datagcmax
data$neiflag = neiflag
data$totalreadcount = totalcount
data <- data[data$neiflag == 0,]
row.names(data) <- seq_len(nrow(data))

#  Fit or load ZINB parameters
if (file.exists(param_file)) {
  message("Loading existing parameters from ", param_file)
  load(param_file)            # brings `optres` into workspace
  par <- optres$par
} else {
  message("Fitting ZINB model (this may take a while)...")
  optres <- optim(
    rnorm(9, mean = -0.4, sd = 0.1),
    CalcZINB,
    control = list(maxit = 10000)
  )
  par <- optres$par
  message("Saving fitted parameters to ", param_file)
  save(optres, file = param_file)
}

# ---- 5. now do the grouping + p-value calculation ----

# 5.1 annotate fusion name & cell type
data$FusionName <- paste(data$V1, data$V2, sep="--")
data$CellType   <- data[[26]]    # 26th column holds your cell type

# 5.2 find every unique (FusionName, CellType) pair
groups <- unique(data[, c("FusionName","CellType")])

# 5.3 for each group, subset the data, rerun p-value fn, collect results
result_list <- vector("list", nrow(groups))
for (i in seq_len(nrow(groups))) {
  fg <- groups$FusionName[i]
  ct <- groups$CellType[i]
  
  idx <- data$FusionName == fg & data$CellType == ct
  sub <- data[idx, , drop=FALSE]
  if (nrow(sub)==0) next
  
  # hack: overwrite global `data` so CalcPValue sees it
  data <<- sub
  pv <- CalcPValueZINB_Fusion_2.0.0_3(par)
  
  # append p-values and group IDs
  sub$PvalueFusion <- pv
  sub$FusionName   <- fg
  sub$CellType     <- ct
  
  result_list[[i]] <- sub
}

# glue them back together
grouped_data <- do.call(rbind, result_list)

# 6) Write out grouped results
write.table(grouped_data, file=out_file,
            sep="\t", quote=FALSE, row.names=FALSE)
message("Done! grouped p-values written to ", out_file)