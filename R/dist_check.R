#' Check distribution of liabilities
#'
#' @param info A rds file made by gen_sim
#' @return qqplot, mean and sample variance for the liability of each family
#' @examples
#' dist_chec(snp_attach("test.rds"))
#' @export
dist_check = function(info){
  l = list()
  c = 0
  for (i in G2$fam){
    c = c + 1
    if (c != 1){
      if (c %% 2 == 0){
        l[[c-(c/2)]] = i
      }
      else{
        l[[c - ((c-1)/2)-1]] = l[[c - ((c-1)/2)-1]] + i
      }
    }
  }
  par(mfrow=c(1,length(l)))
  for (j in l){
    qqnorm(j)
  }
  c = 0
  for (j in l){
    out = c(mean(j), var(j))
    names(out) = c(sprintf("mean_l%s", c), sprintf("variance_l%s", c))
    print(out)
    c = c + 1
  }
}
