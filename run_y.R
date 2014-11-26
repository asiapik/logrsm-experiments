source('calculate_y.R')

subexperiment_y = function(args, expname, stopDesc, i, subsets, stopControl, shifted, optWeights = NULL) {
  datapath = paste(getwd(), '/data/', sep = '')
  inputpath = paste(datapath, 'normal/data_', i, '.csv', sep = '')
    
  print(paste('calculating', i, 'for', subsets, 'experiment' , expname))
  
  data = read.csv(inputpath, header = FALSE)
  class = calculate_y(data, args)
  
  if (optWeights == 't') {
    myweights = calculate_weigths_with_t(class, data)
    fileName = paste(expname, '_t_', sep = '')
  } else if (optWeights == 'cor') {
    myweights = calculate_weigths_with_cor(class, data)
    fileName = paste(expname, '_cor_', sep = '')
  } else if (optWeights == 'univ') {
    myweights = calculate_weigths_with_univ(class, data)
    fileName = paste(expname, '_cor_', sep = '')
  } else {
    myweights = NULL
    fileName = paste(expname, '_nw_', sep = '')
  }
  
  colnames(data) <- c(paste('V_S_', 1:shifted, sep=''), paste('V_N_' , (shifted+1):ncol(data), sep=''))
  header = c('lp',
             paste('positions_V_S_' ,  1:shifted, sep=''), paste('positions_V_N_',  (shifted+1):ncol(data), sep=''),
             paste('scores_V_S_' ,     1:shifted, sep=''), paste('scores_V_N_',     (shifted+1):ncol(data), sep=''), 
             paste('selections_V_S_' , 1:shifted, sep=''), paste('selections_V_N_', (shifted+1):ncol(data), sep=''))
  
  reg = logRSM(class, data, m = subsets, initial_weights = myweights, stopControl = stopControl)
  result = rev(colnames(data)[order(reg$scores)])
  print (result[1:10])
  
  allpath       = paste(fileName, stopDesc, '-', subsets, '-', 'all',   '.csv', sep = '')
  resultpath    = paste(fileName, stopDesc, '-', subsets, '-', 'varnames',   '.csv', sep = '')
  headerpath    = paste(fileName, stopDesc, '-', subsets, '-', 'weights',    '.csv', sep = '')
  selectionpath = paste(fileName, stopDesc, '-', subsets, '-', 'selections', '.csv', sep = '')
  if (i == 1) {
    write.table(t(colnames(data)), file = headerpath, col.names=FALSE, row.names=FALSE, sep = ",", quote = FALSE, append = FALSE)
    write.table(t(header), file = allpath, col.names=FALSE, row.names=FALSE, sep = ",", quote = FALSE, append = FALSE)
  }
  write.table(t(c(i, rank(-reg$scores, ties.method = 'first'), reg$scores, reg$ns)), file= allpath, col.names=FALSE, row.names=FALSE, sep = ",", quote = FALSE, append = TRUE)
  write.table(t(reg$scores),             file = headerpath,    col.names=FALSE, row.names=FALSE, sep = ",", quote = FALSE, append = TRUE)
  write.table(t(result),                 file = resultpath,    col.names=FALSE, row.names=FALSE, sep = ",", quote = FALSE, append = TRUE)
  write.table(t(c(sum(reg$ns), reg$ns)), file = selectionpath, col.names=FALSE, row.names=FALSE, sep = ",", quote = FALSE, append = TRUE)
}

experiment_y = function(args, expname, stopDesc, subsets, stopControl, shifted, optWeights) {
  for (i in 1:100) {
    subexperiment_y(args, expname, stopDesc, i, subsets, stopControl, shifted, optWeights)
  }
}

experiment_with_min_ns = function(args, expname, shifted) {
  experiment_y(args, expname, 'min_ns-30', 4, stopControl = logRSMStop(min_ns = 30), shifted, optWeights = 'NONE')
}

experiment_with_b = function(args, expname, shifted) {
  experiment_y(args, expname, 'b-600', 4, stopControl = logRSMStop(B = 600), shifted, optWeights = 'univ')
  experiment_y(args, expname, 'b-600', 4, stopControl = logRSMStop(B = 600), shifted, optWeights = 't')
}

experiment_all_with_min_ns = function() {
  experiment_with_min_ns(c(0.7, 0.5), '07_05', 2)
  experiment_with_min_ns(c(0.5, 0.3), '05_03', 2)
  experiment_with_min_ns(c(0.7, 0.6, 0.5), '07_06_05', 3)
  experiment_with_min_ns(c(0.5, 0.4, 0.3), '05_04_03', 3)
  experiment_with_min_ns(c(0.7, 0.6, 0.5, 0.4, 0.3), '07_06_05_04_03', 5)
  experiment_with_min_ns(c(0.7, 0.6, 0.5, 0.4, 0.3, 0.2), '07_06_05_04_03_02', 6)
}

experiment_all_with_b = function() {
  experiment_with_b(c(0.7, 0.5), '07_05', 2)
  experiment_with_b(c(0.5, 0.3), '05_03', 2)
  experiment_with_b(c(0.7, 0.6, 0.5), '07_06_05', 3)
  experiment_with_b(c(0.5, 0.4, 0.3), '05_04_03', 3)
  experiment_with_b(c(0.7, 0.6, 0.5, 0.4, 0.3), '07_06_05_04_03', 5)
  experiment_with_b(c(0.7, 0.6, 0.5, 0.4, 0.3, 0.2), '07_06_05_04_03_02', 6)
}
