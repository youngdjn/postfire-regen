extract.single <- function(layer, plots) {
  vals <- extract(layer,plots,method="bilinear")
  return(vals)
}
