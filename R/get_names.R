# prefix = ting før ordet,
# suffix = ting efter ordet. Jeg tror aldrig jeg har hørt om postfix før, men måske det også er en ting.
#' Adds a postfix to given names for each family member
#' @param names Names of the the values to add a postfix to e.g. l_g
#' @param n_sib Amount of siblings
#' @param id If TRUE add a postfix for the subject
#' @param parents If TRUE add postfix for parents
#' @return A vector containing the names for each family member
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

# jeg vil ikke mene, at det er nødvendigt at have l_e(environment liability) med i jeres covariance matrix, men i stedet den *fulde* liability.
# Derudover, så er vi ikke rigtig interesseret i andet end den fulde liability fra familie medlemmerne. Vi har ikke noget genetisk data på dem
# og hvis vi inkluderer flere indgange end nødvendigt i covariance matricen, så tager det bare længere tid til at konvergere.
# I praksis er den tid måske ikke overvældende i dette tilfælde, men det er alligevel værd at overveje.

# som udgangspunkt kan jeg godt lide idéen med denne funktion. Men overvej om i ikke kan lave den på en måde i stil med get_cov..
# evnetuelt bare sæt række og søjle navne på jeres covmat, og brug dem.
