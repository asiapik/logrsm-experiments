generate_data = function () {  
  for (i in 1:100) {
    cols = 500
    rows = 100
    data = matrix(rnorm(rows * cols, mean = 0, sd = 1), rows, cols)
    filename = paste('data_', i, '.csv', sep= '')
    write.table(data, file = filename, col.names = FALSE, row.names=FALSE, sep = ",", quote = FALSE)
  }
}