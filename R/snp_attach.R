#' Attach rds_file
#'
#' For ease of use this function is added
#' It is the exact same as bigsnpr::snp_attach see ?bigsnpr::snp_attach for further documentation
#'
#' @param rds_file An rds file
#' @return The values from the rds file.
#' @export
snp_attach <- function(rds_file){
  return(bigsnpr::snp_attach(rds_file))
}
