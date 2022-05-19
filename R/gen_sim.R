#'  Create and save FBM
#'
#' @param filename Filename for FBM
#' @param N Amount of individuals
#' @param M Amount of SNPs
#' @return A FBM
create_fbm = function(filename, N, M){
  if (file.exists(paste(filename, '.rds', sep=''))){
    obj.bigsnp = bigsnpr::snp_attach(paste(filename, '.rds', sep=''))
    return(obj.bigsnp$genotypes)
  } else{

    G = bigstatsr::FBM.code256(nrow = N,
                               ncol = M,
                               code = c(0L, 1L, 2L, rep(NA_integer_, 256 - 3)),
                               backingfile = filename)

    obj.bigsnp = list(genotypes = G, # genotypes, FBM object
                      map = dplyr::tibble(snp = 1:ncol(G)), # map, i.e. SNP info
                      fam = dplyr::tibble('FID' = 1:nrow(G))) # fam, i.e. info on individuals

    bigsnpr::snp_save(obj.bigsnp)
    return(obj.bigsnp$genotypes)
  }
}



#'  Normalize FBM
#'
#' @param Gx small matrix containing snp data
#' @param probs Probabilities of mutaion
#' @param beta_cut Causal SNPs
#' @return Matrix containing normalized values
normalized_prod = function(Gx, beta, MAF){
  out = c(sweep(sweep(Gx[,], FUN = '-', STATS=2*MAF, MARGIN = 2), FUN='/', STATS=sqrt(2*MAF*(1-MAF)), MARGIN = 2)  %*% beta)
  return(out)
}


#'  Find start and and end index of FBM
#'
#' @param i Iteration from for-loop
#' @param block_size Size of FBM to be processed
#' @return List of start and and end index
get_index = function(i, block_size){
  start = (i-1) * block_size + 1
  end = start + block_size - 1
  return(list('start'=start, 'end'=end))
}



#'  Assign snp for the subject
#'
#' @param G1 SNPs for parent1
#' @param G2 SNPs for parent2
#' @param N Number of subjects
#' @param block_size Size of FBM to be processed
#' @return Matrix of simulated SNPs
assign_snp = function(G1, G2, N=1e5, block_size=1000){
  k = matrix(rnorm(block_size*N, 0, 0.0000001), ncol=block_size, nrow=N, byrow = T)
  G_input = round((G1 + G2)/2 - k)
  return(G_input)
}


#'  Simulate block of snps for each family member
#'
#' @param i Iteration from for loop
#' @param beta Causal snps
#' @param MAF a vector containing minor allele frequencies
#' @param N Number of subjects
#' @param block_size Size of FBM to be processed
#' @return List containing block of simulated snps and genetic liabilities
get_member = function(i, beta, MAF, N=1e5, block_size=1000){

  index = get_index(i, block_size)
  MAF_cut = MAF[index$start:index$end]
  beta_cut = beta[index$start:index$end]

  Gx = matrix(rbinom(block_size*N, 2, MAF_cut), ncol=block_size, nrow=N, byrow = T)

  l_g_x = normalized_prod(Gx, beta_cut, MAF_cut)

  return(list('Gx'=Gx, 'l_g_x'=l_g_x))
}



#'  Simulate block of snps for whole family
#'
#' @param i Iteration from for loop
#' @param G Empty FBM
#' @param beta Causal snps
#' @param MAF a vector containing minor allele frequencies
#' @param n_sibs Number of siblings
#' @param N Number of subjects
#' @param block_size Size of FBM to be processed
#' @return Tibble of simulated SNPs for subject and family's liabilities
#' @importFrom magrittr "%>%"
sim_fam = function(i, G, beta, MAF, N=1e5, n_sib = 0, block_size=1000){

  p1 = get_member(i, beta, MAF, N, block_size)
  p2 = get_member(i, beta, MAF, N, block_size)

  index = get_index(i, block_size)
  G[,index$start:index$end] = assign_snp(p1$Gx, p2$Gx, N, block_size)

  fam_tibble = dplyr::tibble('l_g_partial_p1'=p1$l_g_x, 'l_g_partial_p2'=p2$l_g_x)

  if (n_sib >0) {
    s_tibble = lapply(1:n_sib, function(j){
      Gs = assign_snp(p1$Gx, p2$Gx, N, block_size)
      probs = MAF[index$start:index$end]
      beta_cut = beta[index$start:index$end]
      normalized_prod(Gs, beta_cut, probs)
      }) %>% do.call("cbind", .)

    colnames(s_tibble) = get_names(c("l_g_partial"), id=FALSE, parents = FALSE, n_sib)
    fam_tibble = dplyr::bind_cols(fam_tibble, s_tibble)
  }

  return(fam_tibble)

}



#'  Calculates the sum of colunms in a tibble for each person and value
#'
#' @param data Tibble with colunms to sum over,
#' @param n_sibs Number of siblings
#' @return Collapsed tibble
#' @importFrom magrittr "%>%"
collapse_data = function(data, n_sib){
  search = get_names("", id=FALSE, n_sib)

  iterations = length(search)

  out = lapply(1:iterations, function(i) {
    new_col = rowSums(dplyr::select(data, dplyr::contains(search[i])))
    return(new_col)

  }) %>% dplyr::bind_cols(.)

  colnames(out) = get_names("l_g", id=FALSE, n_sib)

  return(out)

}





#'  Simulate SNP for all subjects and liabilities for family
#'
#' @param filename Filename for FBM
#' @param beta Causal snps
#' @param MAF Probabilities of mutaion
#' @param n_sibs Amount of siblings
#' @param N Amount of individuals
#' @param M Amount of SNPs
#' @param block_size Size of FBM to be processed
#' @return List containing simulated SNPs for all subjects in a FBM and family's liabilities in a tibble
#' @importFrom magrittr "%>%"
#' @export
G_func_fam = function(filename, beta, MAF, N=1e5, M=1e5, n_sib = 0, block_size=1000){

  # filename check
  stopifnot(is_numeric(beta), length(beta) == M)
  stopifnot(is_numeric(MAF), length(MAF) == M)
  stopifnot(is.double(N), N > 0, N%%1 == 0)
  stopifnot(is.double(M), M > 0, M%%1 == 0)
  stopifnot(is.double(n_sib), n_sib >= 0, n_sib%%1 == 0)
  stopifnot(is.double(block_size), block_size > 0, block_size%%1 == 0)


  G = create_fbm(filename, N, M)

  liabil = future.apply::future_lapply(1:(M/block_size), function(i){
    fam_tibble = sim_fam(i, G, beta, MAF, N, n_sib, block_size)
    return(fam_tibble)
  },future.seed = TRUE) %>%
    dplyr::bind_cols(.) %>%
    collapse_data(., n_sib)


  return(list("G"=G, "liabil"=liabil))
}



#'  Simulate SNP for all subjects
#'
#' @param filename Filename for FBM
#' @param MAF A vector containing minor allele frequencies
#' @param n_sib Amount of siblings
#' @param N Amount of individuals
#' @param block_size Size of FBM to be processed
#' @return FBM containing simulated SNPs for all subjects
#' @importFrom magrittr "%>%"
#' @export
G_func_simple = function(filename, MAF, N=1e5, M=1e5, block_size=1000){

  stopifnot(is_numeric(MAF), length(MAF) == M)
  stopifnot(is.double(N), N > 0, N%%1 == 0)
  stopifnot(is.double(M), M > 0, M%%1 == 0)
  stopifnot(is.double(block_size), block_size > 0, block_size%%1 == 0)


  G = create_fbm(filename, N, M)


  iterations = M/block_size
  null_catcher = future.apply::future_lapply(1:iterations, function(i){
    index = get_index(i, block_size)

    MAF_cut = MAF[index$start:index$end]

    G[,index$start:index$end] = matrix(rbinom(block_size*N, 2, MAF_cut), ncol=block_size, nrow=N)

    NULL

  }, future.seed = TRUE) %>% do.call('cbind', .)

  return(G)
}

#' Calculate genetic liabilities and simulate enviromental liabilities for subject and family
#'
#' @param G A FBM containing SNP data
#' @param beta A vector containing the casual SNPs
#' @param MAF A vector containing minor allele frequencies
#' @param liab A tibble containing genetic liabilities for each family member
#' @param n_sib Amount of siblings
#' @param N Amount of individuals
#' @param h_sq heritability
#' @param block_size the size of each iteration
#' @return A tibble containing genetic and enviromental liability for each family member
#' @importFrom magrittr "%>%"
#'
#' @export
liabilities_func_fam = function(G, beta, MAF, liab, N=1e5, n_sib = 0, K=0.05, h2=0.5, block_size = 1000){

  stopifnot(is_tibble(liab))
  stopifnot(is.double(N), N > 0, N%%1 == 0)
  stopifnot(is.double(n_sib), n_sib >= 0, n_sib%%1 == 0)
  stopifnot(is.double(K), K < 1 || K > 0)
  stopifnot(is.double(h2), h2 < 1 || h2 > 0)
  stopifnot(is.double(block_size), block_size > 0, block_size%%1 == 0)

  l_g_0 = lapply(1:(N/block_size), function(i) {
    index = get_index(i, block_size)
    normalized_prod(G[index$start:index$end,],beta, MAF)
     }) %>% do.call("c", .)


  l_out = dplyr::tibble("l_g_0" = l_g_0,
                   "l_f_0"   =  l_g_0 + rnorm(N, 0, sqrt(1-h2)),
                   "l_f_p1"  = liab$l_g_p1 + rnorm(N, 0, sqrt(1-h2)),
                   "l_f_p2"  = liab$l_g_p2 + rnorm(N, 0, sqrt(1-h2))) %>%
    dplyr::bind_cols(liab)


  if (n_sib != 0){
    search = get_names(c("l_g"), id=FALSE, parents = FALSE, n_sib)
    l_out = lapply(1:n_sib, function(i){
      dplyr::select(liab, search[i]) + rnorm(N, 0, sqrt(1-h2))
    }) %>% dplyr::bind_cols(.) %>%
      stats::setNames(., get_names(c("l_f"), id=FALSE, parents = FALSE, n_sib)) %>%
      dplyr::bind_cols(l_out)

  }
  order = get_names(c("l_f","l_g", "pheno"), n_sib)
  T_ =  qnorm(1-K)
  l_out = l_out %>% dplyr::mutate(dplyr::across(.cols = !dplyr::contains("g"),
                          .fns = ~ (.x > T_)-0,
                          .names = "{stringr::str_replace(.col, 'l_f','pheno')}")) %>%
    dplyr::relocate(order)

  return(l_out)

}


#' Calculate genetic liabilities and simulate enviromental liabilities for subject
#'
#' @param G A FBM containing SNP data
#' @param MAF A vector containing minor allele frequencies
#' @param beta a vector containing the casual SNPs
#' @param N Amount of individuals
#' @param h_sq heritability
#' @param block_size the size of each iteration
#' @return A tibble containing genetic and enviromental liability for each subject
#' @importFrom magrittr "%>%"
#' @export
liabilities_func_simple = function(G, beta, MAF, N=1e5,  K=0.05, h2=0.5, block_size = 1000){

  stopifnot(is.double(N), N > 0, N%%1 == 0)
  stopifnot(is.double(K), K < 1 || K > 0)
  stopifnot(is.double(h2), h2 < 1 || h2 > 0)
  stopifnot(is.double(block_size), block_size > 0, block_size%%1 == 0)

  l_g_0 = lapply(1:(N/block_size), function(i) {
    index = get_index(i, block_size)
    normalized_prod(G[index$start:index$end,], beta, MAF)
  }) %>% do.call("c", .)


  T_ =  qnorm(1-K)
  l_out = dplyr::tibble('l_g_0' = l_g_0,
                 'l_f_0'   = l_g_0 + rnorm(N, 0, sqrt(1-h2))) %>%
    dplyr::mutate('pheno_0' = (l_f_0 > T_)-0)

  return(l_out)
}


#'  Simulate vector containing minor allele frequencies
#'
#' @param M Number of SNPs
#' @return Vector containing minor allele frequencies
#' @export
MAF_func = function(M=1e5){

  stopifnot(is.double(M), M > 0, M%%1 == 0)

  out = runif(M, 0.01, 0.49)
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
beta_func = function(M=1e5, h2=0.5, C=1000){


  stopifnot(is.double(M), M > 0, M%%1 == 0)
  stopifnot(is.double(h2), h2 < 1 || h2 > 0)
  stopifnot(is.double(C), C > 0, C%%1 == 0)

  beta = rep(0, M)
  beta[sample(M,C)] = rnorm(C, 0, sqrt(h2/C))
  return (beta)
}




#'  Simulate genetic data
#'
#' @param filename Filename for FBM
#' @param h2 The heritability
#' @param fam Boolean deciding if simulation should include family structure
#' @param n_sibs Number of siblings
#' @param C Amount of causal SNP's
#' @param K The prevalance of trait
#' @param N number of subjects
#' @param M Number of SNPs
#' @param block_size Size of FBM to be processed
#' @return Bigsnp object containing FBM, information about subject and family for liabilities, and generated SNP info
#' @export
gen_sim = function (filename, N=1e5, M=1e5, n_sib = 0, K=0.05, h2=0.5, C=1000, block_size=1000, fam = TRUE) {
  # Make MAF

  stopifnot(is.double(N), N > 0, N%%1 == 0)
  stopifnot(is.double(M), M > 0, M%%1 == 0)
  stopifnot(is.double(n_sib), n_sib >= 0, n_sib%%1 == 0)
  stopifnot(is.double(K), K < 1 || K > 0)
  stopifnot(is.double(h2), h2 < 1 || h2 > 0)
  stopifnot(is.double(block_size), block_size > 0, block_size%%1 == 0)
  stopifnot(is.double(C), C > 0, C%%1 == 0)
  stopifnot(fam==TRUE || fam==FALSE)


  MAF = MAF_func(M)
  beta = beta_func(M, h2, C)

  if (fam){
    G_l = G_func_fam(filename,  beta, MAF, N, M, n_sib, block_size)
    G = G_l$G
    l_out = liabilities_func_fam(G, beta, MAF, G_l$liabil, N, n_sib, K, h2, block_size)

  }
  else{
    G = G_func_simple(filename, MAF, N, M, block_size)
    l_out = liabilities_func_simple(G, beta, MAF, N, K, h2, block_size)
  }

  obj.bigsnp = bigsnpr::snp_attach(paste('test', '.rds', sep=''))
  obj.bigsnp$genotypes = G # genotypes, FBM object
  obj.bigsnp$map = dplyr::tibble(snp = 1:ncol(G), beta, MAF) # map, i.e. SNP info
  obj.bigsnp$fam = dplyr::bind_cols('FID' = 1:nrow(G), l_out) # fam, i.e. info on individuals


  #saving the bigsnp object
  bigsnpr::snp_save(obj.bigsnp)
  print('Succsesfully simulated and saved the data')
}


