

#### General Functions

finite.diff = function(matrix, degree = 2, dx = 1, dy = 1) {
  matrix = as.array(matrix)
  res = array(0, dim = dim(matrix))
  u = array(0, dim = dim(matrix))
  v = array(0, dim = dim(matrix))
  
  
  if (degree == 1) {
    for (i in 3:(nrow(res)-2)) {
      for (j in 3:(ncol(res)-2)) {
        a = 1/12 * matrix[i-2,j] - 2/3 * matrix[i-1,j] + 2/3 * matrix[i+1,j] - 1/12 * matrix[i+2,j]
        b = 1/12 * matrix[i,j-2] - 2/3 * matrix[i,j-1] + 2/3 * matrix[i,j+1] - 1/12 * matrix[i,j+2]
        res[i,j] = (a/dx + b/dy)
        u[i,j] = a / dx
        v[i,j] = b / dy
      }
    }
  }
  
  if (degree == 2) {
    for (i in 3:(nrow(res)-2)) {
      for (j in 3:(ncol(res)-2)) {
        a = -1/12 * matrix[i-2,j] + 4/3 * matrix[i-1,j] - 5/2 * matrix[i,j] + 4/3 * matrix[i+1,j] - 1/12 * matrix[i+2,j]
        b = -1/12 * matrix[i,j-2] + 4/3 * matrix[i,j-1] - 5/2 * matrix[i,j] + 4/3 * matrix[i,j+1] - 1/12 * matrix[i,j+2]
        res[i,j] = a/dx^2 + b/dy^2
        u[i,j] = a / dx^2
        v[i,j] = b / dy^2
      }
    }
  }
  list(res = res, u = u, v = v)
}




get.year = function (x) { as.POSIXlt(x)$year + 1900 }
get.month = function (x) { as.numeric(format(x, '%m')) }
get.day = function (x) { as.numeric(format(x, '%d')) }
