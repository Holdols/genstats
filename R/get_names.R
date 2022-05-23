#' Adds a postfix to given names for each family member
#' @param names Names of the the values to add a postfix to e.g. l_g.
#' @param n_sib Amount of siblings.
#' @param id If TRUE add a postfix for the subject.
#' @param parents If TRUE add postfix for parents.
#' @return A vector containing the names for each family member.
#' @examples
#' get_names(c("l_g", "l_e"), n_sib = 2)
#' @export
get_names <- function(names, n_sib=0, id= TRUE, parents=TRUE){

  len = length(names)
  prefix = c(rep("", len*id), rep("p",len*2*parents), rep("s",n_sib*len))

  n_peop = 1*id + 2*parents + n_sib
  index = c(rep(0,id), rep(c(1:2),parents), rep(c(3:(2+n_sib)),n_sib!=0))

  name = rep("", len*(n_peop))
  j = 1
  for (i in 1:n_peop){
    name[j:(j+len-1)] = paste0((paste(names,prefix[j:(j+len-1)], sep="_")), rep(index[i], len))
    j = j + len
  }
  return(name)
}
