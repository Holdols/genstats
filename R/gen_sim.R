
#'  Normalize fbm
#'
#' @param Gx Fbm containing snp data
#' @param probs Probabilities of mutaion
#' @param beta_cut Causal SNPs
#' @return Normalized values
#' @concept GenSim
#' @noRd
normalized_prod = function(Gx, beta_cut, probs){
  out = c(sweep(sweep(Gx[,], FUN = '-', STATS=2*probs, MARGIN = 2), FUN='/', STATS=sqrt(2*probs*(1-probs)), MARGIN = 2)  %*% beta_cut)
  return(out)
}


#'  Find start and and end index of fbm
#'
#' @param i Iteration from for loop
#' @param block_size Size of fbm to be processed
#' @return List of start and and end index
#' @noRd
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
#' @noRd
get_member = function(i, beta, MAF, N=1e5, block_size=1000){

  index = get_index(i, block_size)
  probs = MAF[index$start:index$end]
  beta_cut = beta[index$start:index$end]

  Gx = matrix(rbinom(block_size*N, 2, probs), ncol=block_size, nrow=N, byrow = T)
  l_g_x = normalized_prod(Gx, beta_cut, probs)

  return(list('Gx'=Gx, 'l_g_x'=l_g_x))
}


#'  Assign snp for the subject
#'
#' @param G1 SNPs for parent1
#' @param G2 SNPs for parent2
#' @param N Number of subjects
#' @param block_size Size of fbm to be processed
#' @return Block of simulated SNPs
#' @noRd
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
#' @noRd
sim_fam = function(i, beta, MAF, n_sib = 0, N=1e5, block_size=1000){
  p1 = get_member(i, beta, MAF, N, block_size)
  p2 = get_member(i, beta, MAF, N, block_size)

  G1 = p1$Gx
  G2 = p2$Gx
  l_g_p1 = p1$l_g_x
  l_g_p2 = p2$l_g_x

  G_input = assign_snp(G1, G2, N, block_size)

  fam_tibble = tibble('l_g_partial_p1'=l_g_p1, 'l_g_partial_p2'=l_g_p2)


  if (n_sib >0) {
    l_g_sibs = rep(list(rep(0,N)),n_sib)
    for (j in 1:n_sib){
      Gs = assign_snp(G1, G2, N, block_size)
      index = get_index(i, block_size)
      probs = MAF[index$start:index$end]
      beta_cut = beta[index$start:index$end]
      temp = normalized_prod(Gs, beta_cut, probs)
      l_g_sibs[[j]] = temp
    }

    s_tibble = do.call(cbind, l_g_sibs)
    col_names_g = get_names(c("l_g_partial"), id=FALSE, parents = FALSE, n_sib)
    colnames(s_tibble) = col_names_g
    fam_tibble = bind_cols(fam_tibble, s_tibble)
  }


  return(list('G_input'= G_input, 'fam_tibble'= fam_tibble))

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
G_func_fam = function(filename, beta, MAF, n_sib = 0, N=1e5, M=1e5, block_size=1000){
  G = FBM.code256(nrow = N,
                  ncol = M,
                  code = c(0L, 1L, 2L, rep(NA_integer_, 256 - 3)),
                  backingfile = filename)


  l_g_p1 = rep(0, M)
  l_g_p2 = rep(0, M)

  iterations = M/block_size

  liabil = future_lapply(1:iterations, function(i){
    fam_sim = sim_fam(i, beta, MAF, n_sib, N, block_size)
    index = get_index(i, block_size)
    G[,index$start:index$end] = fam_sim$G_input
    fam_tibble = fam_sim$fam_tibble


    return(fam_tibble)
  },future.seed = TRUE) %>%
    do.call(bind_cols, .) %>% collapse_data(., n_sib)


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
  G = FBM.code256(nrow = N,
                  ncol = M,
                  code = c(0L, 1L, 2L, rep(NA_integer_, 256 - 3)),
                  backingfile = filename)


  iterations = M/block_size
  future_lapply(1:iterations, function(i){
    current_start = (i-1) * block_size + 1
    current_end = current_start + block_size - 1

    probs = MAF[current_start:current_end]

    G[,current_start:current_end] = matrix(rbinom(block_size*N, 2, probs), ncol=block_size, nrow=N)

  }, future.seed = TRUE) %>% do.call('cbind', .)

  return(G)
}


#'  Calculates the sum of colunms in a tibble for each person and value
#'
#' @param data Tibble with colunms to sum over,
#' @param n_sibs Number of siblings
#' @return Collapsed tibble
#' @noRd
collapse_data = function(data, n_sib){
  search = get_names("", id=FALSE, n_sib)
  expect_names = get_names("l_g", id=FALSE, n_sib)
  col_names = colnames(data)

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
beta_func = function(C=1000, h_sq=0.5, M=1e5){
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
liabilities_func_fam = function(G, MAF, beta, n_sib = 0, N=1e5, h_sq=0.5, block_size = 1000){
  l_g_0 = rep(0, N)
  iterations = N/block_size

  for (i in 1:iterations){
    index = get_index(i, block_size)
    probs = MAF[index$start:index$end]
    beta_cut = beta[index$start:index$end]
    l_g_0[index$start:index$end] = normalized_prod(G[index$start:index$end,], probs, beta_cut)

  }

  l_e_0 = rnorm(N, 0, sqrt(1-h_sq))
  l_e_p1 = rnorm(N, 0, sqrt(1-h_sq))
  l_e_p2 = rnorm(N, 0, sqrt(1-h_sq))

  l_e_out = tibble(l_g_0, l_e_0, l_e_p1, l_e_p2)


  if (n_sib != 0){
    col_names = get_names(c("l_e"), id=FALSE, parents = FALSE, n_sib)
    l_e_sibs = as_tibble(matrix(rnorm(N * n_sib, 0, sqrt(1-h_sq)), ncol = n_sib))
    colnames(l_e_sibs) = col_names



    return(bind_cols(l_e_out,l_e_sibs))
  }
  else{
    return(l_e_out)
  }
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
#'
#' @export
liabilities_func_simple = function(G, MAF, beta, N=1e5, h_sq=0.5, block_size = 1000){
  l_g = rep(0, N)
  iterations = N/block_size

  for (i in 1:iterations){
    current_start = (i-1) * block_size + 1
    current_end = current_start + block_size - 1
    l_g[current_start:current_end] = c(sweep(sweep(G[current_start:current_end,], FUN = '-', STATS=2*MAF, MARGIN = 2), FUN='/',
                                             STATS=sqrt(2*MAF*(1-MAF)), MARGIN = 2)  %*% beta)

  }

  l_e = rnorm(N, 0, sqrt(1-h_sq))

  return(tibble(l_g, l_e))
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
gen_sim = function (filename, h_sq=0.5, fam = TRUE, n_sib = 0, C=1000, K=0.05, N=1e5, M=1e5, block_size=1000) {
  # Make MAF
  MAF = MAF_func(M)
  beta = beta_func(C, h_sq, M)

  if (fam == TRUE){
    G_l = G_func_fam(filename, beta, MAF, n_sib, N, M, block_size)
    G = G_l$G
    liabil = G_l$liabil
    l_g_e = liabilities_func(G, MAF, beta, n_sib, N, h_sq, block_size)
    order = get_names(c("l_g","l_e"), n_sib)
    out = bind_cols(liabil, l_g_e) %>% relocate(order)
  }
  else{
    G = G_func_simple(filename, MAF, b, N, M, block_size)
    out = liabilities_func_simple(G, MAF, beta, N, h_sq, block_size)
  }


  obj.bigsnp = list(genotypes = G, # genotypes, FBM object
                    map = tibble(snp = 1:ncol(G), MAF, beta), # map, i.e. SNP info
                    fam = bind_cols(tibble(FID = 1:nrow(G)), out)) # fam, i.e. info on individuals


  #saving the bigsnp object
  snp_save(obj.bigsnp)

}



