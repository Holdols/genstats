
#'  Normalize fbm
#'
#' @param Gx Fbm containing snp data
#' @param probs Probabilities of mutaion
#' @param beta_cut Causal SNPs
#' @return Normalized values
normalized_prod = function(Gx, beta, MAF){
  out = c(sweep(sweep(Gx[,], FUN = '-', STATS=2*MAF, MARGIN = 2), FUN='/', STATS=sqrt(2*MAF*(1-MAF)), MARGIN = 2)  %*% beta)
  return(out)
}


#'  Find start and and end index of fbm
#'
#' @param i Iteration from for loop
#' @param block_size Size of fbm to be processed
#' @return List of start and and end index
get_index = function(i, block_size){
  start = (i-1) * block_size + 1
  end = start + block_size - 1
  return(list('start'=start, 'end'=end))
}


#'  Simulate block of snps for each family member
#'
#' @param i Iteration from for loop
#' @param beta Causal snps
#' @param MAF a vector containing minor allele frequencies
#' @param N Number of subjects
#' @param block_size Size of fbm to be processed
#' @return List containing block of simulated snps and genetic liabilities
get_member = function(i, beta, MAF, N=1e5, block_size=1000){

  index = get_index(i, block_size)
  MAF_cut = MAF[index$start:index$end]
  beta_cut = beta[index$start:index$end]

  Gx = matrix(rbinom(block_size*N, 2, MAF_cut), ncol=block_size, nrow=N, byrow = T)

  l_g_x = normalized_prod(Gx, beta_cut, MAF_cut)

  return(list('Gx'=Gx, 'l_g_x'=l_g_x))
}


#'  Assign snp for the subject
#'
#' @param G1 SNPs for parent1
#' @param G2 SNPs for parent2
#' @param N Number of subjects
#' @param block_size Size of fbm to be processed
#' @return Block of simulated SNPs
assign_snp = function(G1, G2, N=1e5, block_size=1000){
  k = matrix(rnorm(block_size*N, 0, 0.0000001), ncol=block_size, nrow=N, byrow = T)
  G_input = round((G1 + G2)/2 - k)
  return(G_input)
}



#'  Simulate block of snps for whole family
#'
#' @param i Iteration from for loop
#' @param beta Causal snps
#' @param MAF a vector containing minor allele frequencies
#' @param n_sibs Number of siblings
#' @param N Number of subjects
#' @param block_size Size of fbm to be processed
#' @return Block of simulated SNPs for subject and family's liabilities
sim_fam = function(i, G, beta, MAF, N=1e5, n_sib = 0, block_size=1000){

  p1 = get_member(i, beta, MAF, N, block_size)
  p2 = get_member(i, beta, MAF, N, block_size)

  index = get_index(i, block_size)
  G[,index$start:index$end] = assign_snp(p1$Gx, p2$Gx, N, block_size)

  fam_tibble = tibble('l_g_partial_p1'=p1$l_g_x, 'l_g_partial_p2'=p2$l_g_x)

  if (n_sib >0) {
    s_tibble = lapply(1:n_sib, function(j){
      Gs = assign_snp(G1, G2, N, block_size)
      probs = MAF[index$start:index$end]
      beta_cut = beta[index$start:index$end]
      normalized_prod(Gs, beta_cut, probs)
      }) %>% do.call("cbind", .)

    colnames(s_tibble) = get_names(c("l_g_partial"), id=FALSE, parents = FALSE, n_sib)
    fam_tibble = bind_cols(fam_tibble, s_tibble)
  }

  return(fam_tibble)

}




#'  Simulate SNP for all subjects and liabilities for family
#'
#' @param filename Filename for fbm
#' @param beta Causal snps
#' @param MAF Probabilities of mutaion
#' @param n_sibs Number of siblings
#' @param N Number of subjects
#' @param M Number of SNPs
#' @param block_size Size of fbm to be processed
#' @return Simulated SNPs for all subjects and family's liabilities
#' @export
G_func_fam = function(filename, beta, MAF, N=1e5, M=1e5, n_sib = 0, block_size=1000){
  G = bigstatsr::FBM.code256(nrow = N,
                             ncol = M,
                             code = c(0L, 1L, 2L, rep(NA_integer_, 256 - 3)),
                             backingfile = filename)

  liabil = future.apply::future_lapply(1:(M/block_size), function(i){
    fam_tibble = sim_fam(i, G, beta, MAF, N, n_sib, block_size)
    return(fam_tibble)
  },future.seed = TRUE) %>%
    do.call(bind_cols, .) %>%
    collapse_data(., n_sib)


  return(list("G"=G, "liabil"=liabil))
}


#'  Simulate SNP for all subjects
#'
#' @param filename Filename for fbm
#' @param MAF Probabilities of mutaion
#' @param n_sibs Number of siblings
#' @param N number of subjects
#' @param block_size Size of fbm to be processed
#' @return Simulated SNPs for all subjects
#' @export
G_func_simple = function(filename, MAF, N=1e5, M=1e5, block_size=1000){
  G = bigstatsr::FBM.code256(nrow = N,
                             ncol = M,
                             code = c(0L, 1L, 2L, rep(NA_integer_, 256 - 3)),
                             backingfile = filename)


  iterations = M/block_size
  null_catcher = future.apply::future_lapply(1:iterations, function(i){
    index = get_index(i, block_size)

    MAF_cut = MAF[current_start:current_end]

    G[,current_start:current_end] = matrix(rbinom(block_size*N, 2, MAF_cut), ncol=block_size, nrow=N)

    NULL

  }, future.seed = TRUE) %>% do.call('cbind', .)

  return(G)
}


#'  Calculates the sum of colunms in a tibble for each person and value
#'
#' @param data Tibble with colunms to sum over,
#' @param n_sibs Number of siblings
#' @return Collapsed tibble
collapse_data = function(data, n_sib){
  search = get_names("", id=FALSE, n_sib)
  expect_names = get_names("l_g", id=FALSE, n_sib)

  iterations = length(search)

  out = lapply(1:iterations, function(i) {
    new_col = rowSums(select(data, contains(search[i])))
    return(new_col)

  }) %>% do.call(bind_cols, .)

  colnames(out) = expect_names

  return(out)

}


#' Simulate beta
#'
#' @param C Amount of causal SNPs
#' @param h2 The heritability
#' @param M size of the beta vector
#' @return A vector of size m containing a value from rnorm(1, 0, sqrt(h_sq/C)) \cr
#' a C random places.
#' @export
beta_func = function(M=1e5, h_sq=0.5, C=1000){
  beta = rep(0, M)
  beta[sample(M,C)] = rnorm(C, 0, sqrt(h_sq/C))
  return (beta)
}



#' Calculate genetic liabilities and simulate enviromental liabilities for subject and family
#'
#' @param G A FBM containing SNP data
#' @param MAF a vector containing minor allele frequencies
#' @param beta Causal snps
#' @param n_sib Amount of siblings
#' @param N Amount of individuals
#' @param h_sq heritability
#' @param block_size the size of each iteration
#' @return A tibble containing genetic and enviromental liability for each family member
#'
#' @export
liabilities_func_fam = function(G, beta, MAF, liab, N=1e5, n_sib = 0, K=0.05, h_sq=0.5, block_size = 1000){
  l_g_0 = lapply(1:N/block_size, function(i) {
    index = get_index(i, block_size)
    normalized_prod(G[index$start:index$end,],beta, MAF)
     }) %>% do.call("c", .)


  l_out = tibble("l_g_0" = l_g_0,
                   "l_0"   =  l_g_0 + rnorm(N, 0, sqrt(1-h_sq)),
                   "l_p1"  = liab$l_g_p1 + rnorm(N, 0, sqrt(1-h_sq)),
                   "l_p2"  = liab$l_g_p2 + rnorm(N, 0, sqrt(1-h_sq)))


  if (n_sib != 0){
    search = get_names(c("l_g"), id=FALSE, parents = FALSE, n_sib)
    l_sibs = lapply(1:n_sib, function(i){
      select(liabil, search[i]) + rnorm(N, 0, sqrt(1-h_sq))
    }) %>% bind_cols(.) %>%
      setNames(., get_names(c("l"), id=FALSE, parents = FALSE, n_sib)) %>%
      bind_cols(l_out)

  }
  order = get_names(c("l","l_g", "pheno"), n_sib)
  T_ =  qnorm(1-K)
  l_out %>% mutate(across(.cols != contains("g"),
                          .fns = ~ (.x > T_)-0,
                          .names = "{stringr::str_replace(.col, 'l','pheno')}")) %>%
    relocate(order)


  return(l_out)

}


#' Calculate genetic liabilities and simulate enviromental liabilities for subject
#'
#' @param G A FBM containing SNP data
#' @param MAF a vector containing minor allele frequencies
#' @param beta Causal snps
#' @param N Amount of individuals
#' @param h_sq heritability
#' @param block_size the size of each iteration
#' @return A tibble containing genetic and enviromental liability for each subject
#' @export
liabilities_func_simple = function(G, beta, MAF, N=1e5,  K=0.05, h_sq=0.5, block_size = 1000){

  l_g_0 = lapply(1:N/block_size, function(i) {
    index = get_index(i, block_size)
    normalized_prod(G[index$start:index$end,], beta, MAF)
  }) %>% do.call("c", .)

  T_ =  qnorm(1-K)
  l_out = tibble('l_g' = l_g,
                 'l'   = l_g + rnorm(N, 0, sqrt(1-h_sq))) %>%
    mutate('pheno' = (l > T_)-0)

  return(l_out)
}


#'  Simulate vector containing minor allele frequencies
#'
#' @param M Number of SNPs
#' @return Vector containing minor allele frequencies
#' @export
MAF_func = function(M=1e5){
  out = runif(M, 0.01, 0.49)
  return(out)
}





#'  Simulate genetic data
#'
#' @param filename Filename for fbm
#' @param h2 The heritability
#' @param fam Boolean deciding if simulation should include family structure
#' @param n_sibs Number of siblings
#' @param C Amount of causal SNP's
#' @param K The prevalance of trait
#' @param N number of subjects
#' @param M Number of SNPs
#' @param block_size Size of fbm to be processed
#' @return bigsnp object containing fbm, information about subject and family for liabilities, and generated SNP info
#' @export
gen_sim = function (filename, N=1e5, M=1e5, n_sib = 0, K=0.05, h_sq=0.5, C=1000, block_size=1000, fam = TRUE) {
  # Make MAF
  MAF = MAF_func(M)
  beta = beta_func(M, h_sq, C)

  if (fam){
    G_l = G_func_fam(filename,  beta, MAF, N, M, n_sib, block_size)
    G = G_l$G
    l_out = liabilities_func_fam(G, beta, MAF, liab, N, n_sib, K, h_sq, block_size)

  }
  else{
    G = G_func_simple(filename, MAF, N, M, block_size)
    l_out = liabilities_func_simple(G, beta, MAF, N, K, h_sq, block_size)
  }


  obj.bigsnp = list(genotypes = G, # genotypes, FBM object
                    map = tibble(snp = 1:ncol(G), beta, MAF), # map, i.e. SNP info
                    fam = bind_cols(tibble(FID = 1:nrow(G)), l_out)) # fam, i.e. info on individuals


  #saving the bigsnp object
  bigsnpr::snp_save(obj.bigsnp)

}


