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

generate_all = function() {
  generate_yy(c(0.7, 0.5), '07_05')
  generate_yy(c(0.5, 0.3), '05_03')
  generate_yy(c(0.7, 0.6, 0.5), '07_06_05')
  generate_yy(c(0.5, 0.4, 0.3), '05_04_03')
  generate_yy(c(0.7, 0.6, 0.5, 0.4, 0.3), '07_06_05_04_03')
  generate_yy(c(0.7, 0.6, 0.5, 0.4, 0.3, 0.2), '07_06_05_04_03_02')
}
