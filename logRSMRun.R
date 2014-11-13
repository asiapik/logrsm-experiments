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
  } else if (optWeights == 'cor') {
    myweights = calculate_weigths_with_cor(class, data)
    fileName = 'log_weights_cor_'
  } else {
    myweights = NULL
    fileName = 'log_no_weights_'
  }
  
  colnames(data) <- c(paste('V_S_' , 1:shifted, sep=''), paste('V_N_' , (shifted+1):ncol(data), sep='') )
  
  reg = logRSM(class, data, m = subsets, initial_weights = myweights, stopControl = logRSMStop(B = b))
  result = rev(colnames(data)[order(reg$scores)])
  print (result[1:10])
  
  resultpath    = paste(fileName, b, '-', subsets, '-', 'varnames',   '.csv', sep = '')
  headerpath    = paste(fileName, b, '-', subsets, '-', 'weights',    '.csv', sep = '')
  selectionpath = paste(fileName, b, '-', subsets, '-', 'selections', '.csv', sep = '')
  if (i == 1) {
    write.table(t(colnames(data)), file = headerpath, col.names=FALSE, row.names=FALSE, sep = ",", quote = FALSE, append = TRUE)
  }
  write.table(t(reg$scores), file = headerpath,    col.names=FALSE, row.names=FALSE, sep = ",", quote = FALSE, append = TRUE)
  write.table(t(result),     file = resultpath,    col.names=FALSE, row.names=FALSE, sep = ",", quote = FALSE, append = TRUE)
  write.table(t(reg$ns),     file = selectionpath, col.names=FALSE, row.names=FALSE, sep = ",", quote = FALSE, append = TRUE)
}

source(paste(dirname(thisFile()), '/../logrsm/R/logRSM.R', sep = ''))

spec <- matrix(c(
  'subset',       's', 1, 'integer',    'size of subset (required)',
  'iterations',   'b', 1, 'integer',    'number of iterations (required)',
  'significant',  'g', 1, 'integer',    'number of significant variables',
  'weights',      'w', 1, 'character',  'weights (possible values: none, t, cor)',
  'help',         'h', 0, 'logical',    'this help'
),ncol=5,byrow=T)

opt = getopt(spec);

if ( !is.null(opt$help) || is.null(opt$subset) || is.null(opt$iterations) || is.null(opt$significant) || is.null(opt$weights)) {
  cat(getopt(spec, usage=TRUE))
  q(status=1)
}

if (!is.element(opt$weights, c('none', 'cor', 't'))) {
  cat(paste('unsupported weight calculation:', paste(opt$weights, ',', sep = ''), 'possible values: none, t, cor'))
  q(status=1)
}

for (i in 1:2) {
  print(system.time(subexperiment(opt$subset, i, opt$iterations, opt$significant, opt$weights)))
}
