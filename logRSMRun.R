# install.packages('getopt')
library(getopt)

#' http://stackoverflow.com/questions/1815606/rscript-determine-path-of-the-executing-script
thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(sub(needle, "", cmdArgs[match]))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}

subexperiment = function(subsets, i, b, shifted, optWeights) {
  print(paste('calculating', i, 'for', subsets, 'and', b, 'and weights', optWeights))
  
  inputpath = paste('dane_', i, '.csv', sep = '')
  dane = read.csv(inputpath, header = FALSE)
  class = dane[, 1]
  data = dane[, 2:ncol(dane)]
  
  if (optWeights == 't') {
    myweights = calculate_weigths_with_t(class, data)
    fileName = 'log_weights_t_'
  } else if (optWeights == 'univ') {
    myweights = calculate_weigths_with_univ(class, data)
    fileName = 'log_weights_univ_'
  } else {
    myweights = NULL
    fileName = 'log_no_weights_'
  }
  
  colnames(data) <- c(paste('V_S_' , 1:shifted, sep=''), paste('V_N_' , (shifted+1):ncol(data), sep=''))
  header = c('lp',
             paste('positions_V_S_' , 1:shifted, sep=''),  paste('positions_V_N_',  (shifted+1):ncol(data), sep=''),
             paste('scores_V_S_' , 1:shifted, sep=''),     paste('scores_V_N_',     (shifted+1):ncol(data), sep=''), 
             paste('selections_V_S_' , 1:shifted, sep=''), paste('selections_V_N_', (shifted+1):ncol(data), sep=''))

  reg = logRSM(class, data, m = subsets, initial_weights = myweights, stopControl = logRSMStop(B = b))
  result = rev(colnames(data)[order(reg$scores)])
  print (result[1:10])

  
  allpath    = paste(fileName, b, '-', subsets, '-', 'all',   '.csv', sep = '')
  resultpath    = paste(fileName, b, '-', subsets, '-', 'varnames',   '.csv', sep = '')
  headerpath    = paste(fileName, b, '-', subsets, '-', 'weights',    '.csv', sep = '')
  selectionpath = paste(fileName, b, '-', subsets, '-', 'selections', '.csv', sep = '')
  if (i == 1) {
    write.table(t(colnames(data)), file = headerpath, col.names=FALSE, row.names=FALSE, sep = ",", quote = FALSE, append = FALSE)
    write.table(t(header), file = allpath, col.names=FALSE, row.names=FALSE, sep = ",", quote = FALSE, append = FALSE)
  }
  print(rank(reg$scores))
  print(order(reg$scores))
  write.table(t(c(i, rank(-reg$scores, ties.method = 'first'), reg$scores, reg$ns)), file= allpath, col.names=FALSE, row.names=FALSE, sep = ",", quote = FALSE, append = TRUE)
  write.table(t(reg$scores), file = headerpath,    col.names=FALSE, row.names=FALSE, sep = ",", quote = FALSE, append = TRUE)
  write.table(t(result),     file = resultpath,    col.names=FALSE, row.names=FALSE, sep = ",", quote = FALSE, append = TRUE)
  write.table(t(reg$ns),     file = selectionpath, col.names=FALSE, row.names=FALSE, sep = ",", quote = FALSE, append = TRUE)
}

source(paste(dirname(thisFile()), '/../logrsm/R/logRSM.R', sep = ''))

spec <- matrix(c(
  'subset',       's', 1, 'integer',    'size of subset (required)',
  'iterations',   'b', 1, 'integer',    'number of iterations (required)',
  'significant',  'g', 1, 'integer',    'number of significant variables',
  'weights',      'w', 1, 'character',  'weights (possible values: none, t, cor, univ)',
  'help',         'h', 0, 'logical',    'this help'
),ncol=5,byrow=T)

docli = function () {
  opt = getopt(spec);

  if ( !is.null(opt$help) || is.null(opt$subset) || is.null(opt$iterations) || is.null(opt$significant) || is.null(opt$weights)) {  
    cat(getopt(spec, usage=TRUE))
    q(status=1)
  }

  if (!is.element(opt$weights, c('none', 'cor', 't', 'univ'))) {
    cat(paste('unsupported weight calculation:', paste(opt$weights, ',', sep = ''), 'possible values: none, t, cor, univ'))
    q(status=1)
  }

  for (i in 1:1) {
    print(system.time(subexperiment(opt$subset, i, opt$iterations, opt$significant, opt$weights)))
  }
}

run450 = function() {
  for (i in 1:100) {
    print(System.time(subexperiment(5, i, 450, 2, 'univ')))
  }
}

run4500 = function() {
  for (i in 1:100) {
    print(System.time(subexperiment(5, i, 4500, 2, 'univ')))
  }
}

runRFInternal = function(i, shifted) {
  print(paste('calculating', i))
  
  inputpath = paste('dane_', i, '.csv', sep = '')
  dane = read.csv(inputpath, header = FALSE)
  class = dane[, 1]
  data = dane[, 2:ncol(dane)]
  
  colnames(data) <- c(paste('V_S_' , 1:shifted, sep=''), paste('V_N_' , (shifted+1):ncol(data), sep=''))
  header = c('lp',
             paste('IncMSE V_S_' , 1:shifted, sep=''), paste('IncMSE V_N_',  (shifted+1):ncol(data), sep=''),
             paste('IncNodePurity V_S_' , 1:shifted, sep=''), paste('IncNodePurity V_N_',     (shifted+1):ncol(data), sep='')) 
  
  rf = randomForest(y = t(class), x = data, keep.forest = FALSE, importance = TRUE)
  
  if (i == 1) {
    write.table(t(colnames(data)), file = 'result.csv', col.names=FALSE, row.names=FALSE, sep = ",", quote = FALSE, append = FALSE)
  }
  imp = importance(rf)
  write.table(c(i, t(imp[,1]), t(imp[,2])), file = 'result.csv', col.names=FALSE, row.names=FALSE, sep = ",", quote = FALSE, append = TRUE)
}

runRF = function() {
  for (i in 1:100) {
    print(system.time(runRFInternal(i, 20)))
  }
}

calculate_y = function(dane, args) {
  vec = numeric(ncol(dane))
  for (k in 1:length(args)) {
    vec[k] = args[k]
  }
  y = numeric(nrow(dane))
  for (j in 1:nrow(dane)) {
    y[j] = ifelse (sum(vec * dane[j,1:ncol(dane)]) > 0, 1, 0)
  }
  return (y)  
}

generate_y = function(i, args, expname) {
  datapath = paste(getwd(), '/data/', sep = '')
  inputpath = paste(datapath, 'normal/data_', i, '.csv', sep = '')
  dane = read.csv(inputpath, header = FALSE)
  yfilename = paste(getwd(), '/experiments/', expname, '/y_', i, '.csv', sep = '')
  y = calculate_y(dane, args)
  write.table(y, file = yfilename, col.names = FALSE, row.names=FALSE, sep = ",", quote = FALSE)
}

generate_yy = function(args, expname) {
  exppath = paste(getwd(), '/experiments/', expname, sep = '')
  dir.create(exppath)
  for (i in 1:100) {
    generate_y(i, args, expname)
  }
}

generate_yyy = function() {
  generate_yy(c(0.7, 0.5), '07_05')
  generate_yy(c(0.5, 0.3), '05_03')
  generate_yy(c(0.7, 0.6, 0.5), '07_06_05')
  generate_yy(c(0.5, 0.4, 0.3), '05_04_03')
  generate_yy(c(0.7, 0.6, 0.5, 0.4, 0.3), '07_06_05_04_03')
  generate_yy(c(0.7, 0.6, 0.5, 0.4, 0.3, 0.2), '07_06_05_04_03_02')
}
#myexp(1, c(3,5,6), '03_05_06')
