#' Local version of High-definition likelihood inference of genetic correlations (HDL-L)
#'
#' The function returns the estimate and standard error of the genetic correlation between two traits based on GWAS summary statistics.
#'
#' @param gwas1.df A data frame including GWAS summary statistics of genetic variants for trait 1.
#' The input data frame should include following columns: SNP, SNP ID; A1, effect allele; A2, reference allele;
#' N, sample size; Z, z-score; If Z is not given, alternatively, you may provide: b, estimate of marginal effect in GWAS; se, standard error of the estimates of marginal effects in GWAS. If the GWAS is based on logistic regression, `b` should be the logarithm of OR (odds ratio) and `se` is the standard error of log(OR). Notice: SNPs with missing N will be removed.
#' @param gwas2.df A data frame including GWAS summary statistics of genetic variants for trait 2.
#' The input data frame should include following columns: SNP, SNP ID; A1, effect allele; A2, reference allele;
#' N, sample size; Z, z-score; If Z is not given, alternatively, you may provide: b, estimate of marginal effect in GWAS; se, standard error of the estimates of marginal effects in GWAS. If the GWAS is based on logistic regression, `b` should be the logarithm of OR (odds ratio) and `se` is the standard error of log(OR). Notice: SNPs with missing N will be removed.
#' @param Trait1name The trait name for gwas1.df
#' @param Trait2name The trait name for gwas2.df
#' @param LD.path Path to the directory where linkage disequilibrium (LD) information is stored.
#' @param bim.path Path to the directory where the decompressed bim.path.zip file is stored.
#' @param Nref Sample size of the reference sample where LD is computed. If the default UK Biobank reference sample is used, Nref = 335272
#' @param N0 Number of individuals included in both cohorts.
#' @param output.file Where the log and results should be written. If you do not specify a file, the log will be printed on the console.
#' @param eigen.cut Specifies which eigenvalues and eigenvectors in each LD score matrix should be used for HDL. The default value is 0.99. Users are allowed to specify a numeric value between 0 and 1 for eigen.cut.
#' @param intercept.output Logical, FALSE by default. Should the intercept terms be included in estimates.df?
#' @param fill.missing.N If NULL (default), the SNPs with missing N are removed. One of "median", "min" or "max" can be given so that the missing N will be filled accordingly. For example, "median" means the missing N are filled with the median N of the SNPs with available N.
#' @param lim Tolerance limitation, default lim = exp(-18).
#' @param chr The chromosome to which the region belongs.
#' @param piece The piece of the genome to which the region belongs. The whole genome is divided into 2,468 smaller, semi-independent blocks, each defined by LD calculated by Plink. The SNP information in each local region is included in these two data sets: "UKB_snp_counter_imputed.RData" and "UKB_snp_list_imputed_vector.RData".
#' @param alpha The significance level used to compute the cut-off value 'c' for the likelihood-based confidence interval (LCI). The cut-off 'c' is calculated using the formula: c = exp(-qchisq(1 - alpha, 1) / 2). For example, when alpha = 0.05, it corresponds to a 95% LCI.
#' @note Users can download the precomputed eigenvalues and eigenvectors of LD correlation matrices for European ancestry population. The download link can be found at https://doi.org/10.5281/zenodo.11001214
#' These are the LD matrices and their eigen-decomposition from 335,272 genomic British UK Biobank individuals.
#'
#' @return A dataframe is returned with:
#'
#' \itemize{
#' \item{Trait1 }{The name of trait 1.}
#' \item{Trait2 }{The name of trait 2.}
#' \item{chr }{Which chromosome.}
#' \item{piece }{which piece of the chromosome}
#' \item{estimates.df }{A detailed matrix includes the estimates, P-value or likelihood-based intervals of heritability, genetic covariance and genetic correlation.}
#' \item{eigen.use }{The eigen.cut used in computation.}
#' }
#'
#' @author Yuying Li, Xia Shen
#'
#' @references
#' Li, Y., Pawitan, Y. & Shen, X. An enhanced framework for local genetic correlation analysis. Nat Genet 57, 1053–1058 (2025). https://doi.org/10.1038/s41588-025-02123-3
#' Li Y, Zhai R, Yang Z, Li T, Pawitan Y, Shen X (2025). High-definition likelihood inference of colocalization reveals protein biomarkers for human complex diseases
#' @seealso
#' HDL tutorial: https://github.com/zhenin/HDL
#' @export
#'


HDL.C <-function(gwas1.df, gwas2.df, Trait1name, Trait2name, LD.path, bim.path, Nref = 335272, N0, output.file = "", eigen.cut = 0.99, intercept.output = FALSE, fill.missing.N = NULL, lim = exp(-18), chr, proname, alpha = 0.05){
    piece = proname
    if(output.file != ""){
      if(file.exists(output.file) == T){
        system(paste0("rm ",output.file))
      }
    }
    time.start <- date()
    cat("Analysis starts on",time.start,"\n")
    if(output.file != ""){
      cat("Analysis starts on",time.start,"\n", file = output.file, append = T)
    }

    if(is.na(as.numeric(eigen.cut))){
      error.message <- "The input of eigen.cut has to be a number between 0 and 1. \n"
      if(output.file != ""){
        cat(error.message, file = output.file, append = T)
      }
      stop(error.message)
    }

    #Likelihood ratio function for heritability estimation
    llfun <-function(param, N, M, Nref, lam, bstar, lim){
  h2 = param[1]
  int= param[2]
  lamh2 = h2/M*lam^2 - h2*lam/Nref + int*lam/N
  lamh2 = ifelse(lamh2<lim, lim,lamh2)
  ll =-1/2*( sum(log(lamh2)) + sum(bstar^2/(lamh2)))
  return(ll)
}

  llfun0 <-function(int, N, M, Nref, lam, bstar, lim){
  h2 = 0
  lamh2 = h2/M*lam^2 - h2*lam/Nref + int*lam/N
  lamh2 = ifelse(lamh2<lim, lim,lamh2)
  ll =-1/2*( sum(log(lamh2)) + sum(bstar^2/(lamh2)))
  return(ll)
}


    #Likelihood ratio function for genetic covariance estimation
    llfun.gcov.part.2 = function(param,h11,h22,rho12, M, N1, N2, N0, Nref,
                             lam1, lam2, bstar1, bstar2, lim){
  h12 = param[1]
  int = param[2]
  ## sample fractions
  p1 = N0/N1; p2= N0/N2
  ## must follow the formula for lamh2 used in llfun4
  lam11 =   h11[1]/M*lam1^2 - h11[1]*lam1/Nref + h11[2]*lam1/N1
  lam11 = ifelse(lam11<lim, lim,lam11)
  lam22 = h22[1]/M*lam2^2 - h22[1]*lam2/Nref + h22[2]*lam2/N2
  lam22 = ifelse(lam22<lim, lim,lam22)
  #lam12 = h12/M*lam1*lam2 - p1*p2*h12*lam1/Nref + p1*p2*int*lam1/N0
  if (N0>0) lam12 = h12/M*lam1*lam2 + p1*p2*int*lam1/N0  ## key change here
  if (N0==0) lam12 = h12/M*lam1*lam2
  ##  resid of bstar2 ~bstar1
  ustar = bstar2 - lam12/lam11*bstar1  ## see note
  lam22.1 = lam22 - lam12^2/lam11
  lam22.1 = ifelse(lam22.1<lim, lim,lam22.1)
  ll = -1/2*(sum(log(lam22.1)) + sum(ustar^2/(lam22.1)))
  return(ll)
}

llfun0.gcov.part.2 = function(int,h11,h22,rho12, M, N1, N2, N0, Nref,
                              lam1, lam2, bstar1, bstar2, lim){
  h12 = 0
  ## sample fractions
  p1 = N0/N1; p2= N0/N2
  ## must follow the formula for lamh2 used in llfun4
  lam11 =   h11[1]/M*lam1^2 - h11[1]*lam1/Nref + h11[2]*lam1/N1
  lam11 = ifelse(lam11<lim, lim,lam11)
  lam22 = h22[1]/M*lam2^2 - h22[1]*lam2/Nref + h22[2]*lam2/N2
  lam22 = ifelse(lam22<lim, lim,lam22)
  #lam12 = h12/M*lam1*lam2 - p1*p2*h12*lam1/Nref + p1*p2*int*lam1/N0
  if (N0>0) lam12 = h12/M*lam1*lam2 + p1*p2*int*lam1/N0  ## key change here
  if (N0==0) lam12 = h12/M*lam1*lam2
  ##  resid of bstar2 ~bstar1
  ustar = bstar2 - lam12/lam11*bstar1  ## see note
  lam22.1 = lam22 - lam12^2/lam11
  lam22.1 = ifelse(lam22.1<lim, lim,lam22.1)
  ll = -1/2*(sum(log(lam22.1)) + sum(ustar^2/(lam22.1)))
  return(ll)
}



    # Load SNP list and counter files with error handling
    #load(paste0(LD.path, "HDLL_LOC_snps.RData"))

# Function to find LD RData file based on chromosome (chr) and piece
find_LD_rda_file <- function(LD.files, chr, piece) {
  pattern1 <- paste0("ukb_chr", chr, "\\.", piece, ".*_LDSVD\\.rda$")
  pattern2 <- paste0("ukb_", chr, "\\.", piece, ".*_LDSVD\\.rda$")

  # Try the first pattern
  LD_rda_file <- LD.files[grep(x = LD.files, pattern = pattern1)]

  # If no file is found, try the second pattern
  if (length(LD_rda_file) == 0) {
    LD_rda_file <- LD.files[grep(x = LD.files, pattern = pattern2)]
  }

  # Return the first matching file, or NULL if none found
  if (length(LD_rda_file) > 0) {
    return(LD_rda_file[1])
  } else {
    return(NULL)
  }
}

# List all files in LD.path
LD.files <- list.files(LD.path, full.names = TRUE)

# Use function to find the specific LD file
LD_rda_file <- find_LD_rda_file(LD.files, chr, piece)

# Load the LD file if found
if (!is.null(LD_rda_file)) {
  load(LD_rda_file)
} else {
  stop("No matching LD file found in LD.path.")
}


    bim.files <- list.files(bim.path)

    LD_bim_file <- bim.files[grep(x = bim.files, pattern = paste0("ukb_",chr,".",piece, "[\\._].*bim"))]
if(length(LD_bim_file)==0){
     LD_bim_file <- bim.files[grep(x = bim.files, pattern = paste0("ukb_chr",chr,".",piece, "[\\._].*bim"))]

}

    snps.ref.df <- read.table(paste(bim.path, LD_bim_file, sep = "/"))

    colnames(snps.ref.df) <- c("chr","id","non","pos","A1","A2")
    snps.ref <- snps.name.list  <- snps.ref.df$id
    A2.ref <- snps.ref.df$A2
    names(A2.ref) <- snps.ref
    if(is.data.frame(gwas1.df)!=TRUE){gwas1.df = as.data.frame(gwas1.df)}
    if(is.data.frame(gwas2.df)!=TRUE){gwas2.df = as.data.frame(gwas2.df)}

    gwas1.df <- gwas1.df %>% filter(SNP %in% snps.name.list)
    gwas2.df <- gwas2.df %>% filter(SNP %in% snps.name.list)

    gwas1.df$A1 <- toupper(as.character(gwas1.df$A1))
    gwas1.df$A2 <- toupper(as.character(gwas1.df$A2))
    gwas2.df$A1 <- toupper(as.character(gwas2.df$A1))
    gwas2.df$A2 <- toupper(as.character(gwas2.df$A2))

    ## Check the format of summary statistics ##
    if(!("Z" %in% colnames(gwas1.df))){
      if(("b" %in% colnames(gwas1.df)) && ("se" %in% colnames(gwas1.df))){
        if(abs(median(gwas1.df$b) - 1) < 0.1){
          cat("Taking log(b) in GWAS 1 because b is likely to be OR in stead of log(OR). \n")
          if(output.file != ""){
            cat("Taking log(b) in GWAS 1 because b is likely to be OR in stead of log(OR). \n", file = output.file, append = T)
          }
          gwas1.df$Z <- log(gwas1.df$b) / gwas1.df$se
        } else{
          gwas1.df$Z <- gwas1.df$b / gwas1.df$se
        }
      } else{
        error.message <- "Z is not available in GWAS 1, meanwhile either b or se is missing. Please check."
        if(output.file != ""){
          cat(error.message, file = output.file, append = T)
        }
        stop(error.message)

      }
    }

    if(!("Z" %in% colnames(gwas2.df))){
      if(("b" %in% colnames(gwas2.df)) && ("se" %in% colnames(gwas2.df))){
        if(abs(median(gwas2.df$b) - 1) < 0.1){
          cat("Taking log(b) in GWAS 2 because b is likely to be OR in stead of log(OR). \n")
          if(output.file != ""){
            cat("Taking log(b) in GWAS 2 because b is likely to be OR in stead of log(OR). \n", file = output.file, append = T)
          }
          gwas2.df$Z <- log(gwas2.df$b) / gwas2.df$se
        } else{
          gwas2.df$Z <- gwas2.df$b / gwas2.df$se
        }
      } else{
        error.message <- "Z is not available in GWAS 2, meanwhile either b or se is missing. Please check."
        if(output.file != ""){
          cat(error.message, file = output.file, append = T)
        }
        stop(error.message)

      }
    }

    k1.0 <- length(unique(gwas1.df$SNP))
    k2.0 <- length(unique(gwas2.df$SNP))

    ## Imputed data may have missing N ##
    if(is.null(fill.missing.N)){
      gwas1.df <- gwas1.df %>% filter(!is.na(Z), !is.na(N))
      gwas2.df <- gwas2.df %>% filter(!is.na(Z), !is.na(N))
    } else if(fill.missing.N == "min"){
      gwas1.df <- gwas1.df %>% filter(!is.na(Z))
      gwas1.df$N[is.na(gwas1.df$N)] <- min(gwas1.df$N, na.rm = T)

      gwas2.df <- gwas2.df %>% filter(!is.na(Z))
      gwas2.df$N[is.na(gwas2.df$N)] <- min(gwas2.df$N, na.rm = T)
    } else if(fill.missing.N == "max"){
      gwas1.df <- gwas1.df %>% filter(!is.na(Z))
      gwas1.df$N[is.na(gwas1.df$N)] <- max(gwas1.df$N, na.rm = T)

      gwas2.df <- gwas2.df %>% filter(!is.na(Z))
      gwas2.df$N[is.na(gwas2.df$N)] <- max(gwas2.df$N, na.rm = T)
    } else if(fill.missing.N == "median"){
      gwas1.df <- gwas1.df %>% filter(!is.na(Z))
      gwas1.df$N[is.na(gwas1.df$N)] <- median(gwas1.df$N, na.rm = T)

      gwas2.df <- gwas2.df %>% filter(!is.na(Z))
      gwas2.df$N[is.na(gwas2.df$N)] <- median(gwas2.df$N, na.rm = T)
    } else{
      error.message <- "If given, the argument fill.missing.N can only be one of below: 'min', 'max', 'median'."
      if(output.file != ""){
        cat(error.message, file = output.file, append = T)
      }
      stop(error.message)
    }

    k1 <- length(unique(gwas1.df$SNP))
    k2 <- length(unique(gwas2.df$SNP))
    k1.percent <- paste("(",round(100*k1 / length(snps.name.list), 2), "%)", sep="")
    k2.percent <- paste("(",round(100*k2 / length(snps.name.list), 2), "%)", sep="")

    cat(k1.0 - k1, "SNPs were removed in GWAS 1 due to missing N or missing test statistic."," \n")
    cat(k2.0 - k2, "SNPs were removed in GWAS 2 due to missing N or missing test statistic."," \n")
    if(output.file != ""){
      cat(k1.0 - k1, " SNPs were removed in GWAS 1 due to missing N or missing test statistic."," \n", file = output.file, append = T)
      cat(k2.0 - k2, " SNPs were removed in GWAS 2 due to missing N or missing test statistic."," \n", file = output.file, append = T)
    }

    cat(k1, "out of", length(snps.name.list), k1.percent, "SNPs in reference panel are available in GWAS 1."," \n")
    cat(k2, "out of", length(snps.name.list), k2.percent, "SNPs in reference panel are available in GWAS 2."," \n")
    if(output.file != ""){
      cat(k1, "out of", length(snps.name.list), k1.percent, "SNPs in reference panel are available in GWAS 1."," \n", file = output.file, append = T)
      cat(k2, "out of", length(snps.name.list), k2.percent, "SNPs in reference panel are available in GWAS 2."," \n", file = output.file, append = T)
    }
    if(k1 < length(snps.name.list)*0.99){
      error.message <- "Warning: More than 1% SNPs in reference panel are missed in GWAS 1. This may generate bias in estimation. Please make sure that you are using correct reference panel.  \n"
      if(output.file != ""){
        cat(error.message, file = output.file, append = T)
      }
      cat(error.message)
    }
    if(k2 < length(snps.name.list)*0.99){
      error.message <- "Warning: More than 1% SNPs in reference panel are missed in GWAS 2. This may generate bias in estimation. Please make sure that you are using correct reference panel.  \n"
      if(output.file != ""){
        cat(error.message, file = output.file, append = T)
      }
      cat(error.message)
    }

    ## stats
    N1 <- median(gwas1.df$N)
    N2 <- median(gwas2.df$N)
    N <- sqrt(N1)*sqrt(N2)
    p1 <- N0/N1
    p2 <- N0/N2

  rho12 <- suppressWarnings(inner_join(gwas1.df %>% select(SNP, Z) %>% filter(abs(Z)<8), gwas2.df %>% select(SNP, Z) %>% filter(abs(Z)<8), by = "SNP") %>% summarise(x=cor(Z.x, Z.y, use = "complete.obs")) %>% unlist)

    bstar1.v <- bstar2.v <- lam.v <- list()

    gwas1.df.subset <- gwas1.df %>% filter(SNP %in% snps.ref) %>% distinct(SNP, A1, A2, .keep_all = TRUE)

    ## Check if there are multiallelic or duplicated SNPs
        if(any(duplicated(gwas1.df.subset$SNP)) == TRUE){
          gwas1.df.subset.duplicated <- gwas1.df.subset %>%
            mutate(row.num = 1:n()) %>%
            filter(SNP == SNP[duplicated(SNP)]) %>%
            mutate(SNP_A1_A2 = paste(SNP, A1, A2, sep = "_"))
          snps.ref.df.duplicated <- snps.ref.df %>%
            filter(id %in% gwas1.df.subset.duplicated$SNP)
          SNP_A1_A2.valid <- c(paste(snps.ref.df.duplicated$id, snps.ref.df.duplicated$A1, snps.ref.df.duplicated$A2, sep = "_"),
                               paste(snps.ref.df.duplicated$id, snps.ref.df.duplicated$A2, snps.ref.df.duplicated$A1, sep = "_"))
          row.remove <- gwas1.df.subset.duplicated %>% filter(!(SNP_A1_A2 %in% SNP_A1_A2.valid)) %>% select(row.num) %>% unlist()
          gwas1.df.subset <- gwas1.df.subset[-row.remove,]
        }

        bhat1.raw <- gwas1.df.subset[, "Z"] / sqrt(gwas1.df.subset[, "N"])
        A2.gwas1 <- gwas1.df.subset[, "A2"]
        names(bhat1.raw) <- names(A2.gwas1) <- gwas1.df.subset$SNP
        idx.sign1 <- A2.gwas1 == A2.ref[names(A2.gwas1)]
        bhat1.raw <- bhat1.raw*(2*as.numeric(idx.sign1)-1)


        gwas2.df.subset <- gwas2.df %>% filter(SNP %in% snps.ref) %>% distinct(SNP, A1, A2, .keep_all = TRUE)

        ## Check if there are multiallelic or duplicated SNPs
        if(any(duplicated(gwas2.df.subset$SNP)) == TRUE){
          gwas2.df.subset.duplicated <- gwas2.df.subset %>%
            mutate(row.num = 1:n()) %>%
            filter(SNP == SNP[duplicated(SNP)]) %>%
            mutate(SNP_A1_A2 = paste(SNP, A1, A2, sep = "_"))
          snps.ref.df.duplicated <- snps.ref.df %>%
            filter(id %in% gwas2.df.subset.duplicated$SNP)
          SNP_A1_A2.valid <- c(paste(snps.ref.df.duplicated$id, snps.ref.df.duplicated$A1, snps.ref.df.duplicated$A2, sep = "_"),
                               paste(snps.ref.df.duplicated$id, snps.ref.df.duplicated$A2, snps.ref.df.duplicated$A1, sep = "_"))
          row.remove <- gwas2.df.subset.duplicated %>% filter(!(SNP_A1_A2 %in% SNP_A1_A2.valid)) %>% select(row.num) %>% unlist()
          gwas2.df.subset <- gwas2.df.subset[-row.remove,]
        }
        bhat2.raw <- gwas2.df.subset[, "Z"] / sqrt(gwas2.df.subset[, "N"])
        A2.gwas2 <- gwas2.df.subset[, "A2"]
        names(bhat2.raw) <- names(A2.gwas2) <- gwas2.df.subset$SNP
        idx.sign2 <- A2.gwas2 == A2.ref[names(A2.gwas2)]
        bhat2.raw <- bhat2.raw*(2*as.numeric(idx.sign2)-1)

        M <- length(LDsc)
        bhat1 <- bhat2 <- numeric(M)
        names(bhat1) <- names(bhat2) <- snps.ref
        bhat1[names(bhat1.raw)] <- bhat1.raw
        bhat2[names(bhat2.raw)] <- bhat2.raw


        a11 <- bhat1**2
        a22 <- bhat2**2
        a12 <- bhat1*bhat2

        reg = lm(a11~ LDsc)
        h11.ols <- c(summary(reg)$coef[1:2,1:2]*c(N1,M))

        reg = lm(a22~ LDsc)
        h22.ols <- c(summary(reg)$coef[1:2,1:2]*c(N2,M))

        reg = lm(a12~ LDsc)
        if (N0>0) h12.ols = c(summary(reg)$coef[1:2,1:2]*c((N0/p1/p2),M))
        if (N0==0) h12.ols = c(summary(reg)$coef[1:2,1:2]*c(N,M))

        #  ................................ weighted LS: use estimated h2
        # vars from Bulik-Sullivan

        h11v = (h11.ols[2]*LDsc/M + 1/N1)^2
        h22v = (h22.ols[2]*LDsc/M + 1/N2)^2

        reg = lm(a11~ LDsc, weight=1/h11v)
        h11.wls <- c(summary(reg)$coef[1:2,1:2]*c(N1,M))

        reg = lm(a22~ LDsc, weight=1/h22v)
        h22.wls <- c(summary(reg)$coef[1:2,1:2]*c(N2,M))

        if (N0>0) h12v = sqrt(h11v*h22v) + (h12.ols[2]*LDsc/M + p1*p2*rho12/N0)^2
        if (N0==0) h12v = sqrt(h11v*h22v) + (h12.ols[2]*LDsc/M)^2

        reg = lm(a12~ LDsc, weight=1/h12v)
        if (N0>0) h12.wls = c(summary(reg)$coef[1:2,1:2]*c((N0/p1/p2),M))
        if (N0==0) h12.wls = c(summary(reg)$coef[1:2,1:2]*c(N,M))

        ## .................................  likelihood based
        ## ....  estimate h2s
        bstar1 = crossprod(V,bhat1)  ##
        bstar2 = crossprod(V,bhat2)  ##


    ### Select eigen value and eigen vector
    M.ref <- length(snps.ref)
    eigen_select.fun <- function(x,k){
      return(x[1:k])
    }

    eigen_select_num.fun <- function(eigen.percent, cut.percent){
      if(!(eigen.percent[length(eigen.percent)] > cut.percent)){
        return(length(eigen.percent))
      } else{
        return(which(eigen.percent > cut.percent)[1])
      }
    }
        bstar1.v <- c(bstar1.v, list(bstar1))
        bstar2.v <- c(bstar2.v, list(bstar2))
        lam.v <- c(lam.v, list(lam))

      eigen.cut <- as.numeric(eigen.cut)
      eigen.use <- eigen.cut
      eigen.num.v.cut <- c()
      nsnps.v <- length(snps.ref)
      for(i in 1:length(nsnps.v)){
        lam.i <- lam.v[[i]]
        eigen.percent <- numeric(length(lam.i))
        temp <- 0
        for(j in 1:length(lam.i)){
          temp <- temp + lam.i[j]
          eigen.percent[j] <- temp/nsnps.v[i]
        }
        eigen.num.cut <- eigen_select_num.fun(eigen.percent, eigen.cut)
        eigen.num.v.cut <- c(eigen.num.v.cut, eigen.num.cut)
      }

      lam.v.cut <- mapply(eigen_select.fun,lam.v,eigen.num.v.cut)
      bstar1.v.cut <- mapply(eigen_select.fun,bstar1.v,eigen.num.v.cut)
      bstar2.v.cut <- mapply(eigen_select.fun,bstar2.v,eigen.num.v.cut)

    ##Estimate heritabitlity for trait 1
    starting_values_1 <- c(h11.wls[2],0,0.5)
    starting_values_2 <- c(1,0.5,1.5, 0)
    best_value <- -Inf
    ndeps_values <- c(1e-5, 1e-8, 1e-16)
    converged = FALSE
for (start_val_1 in starting_values_1) {
  for (start_val_2 in starting_values_2) {
    for (n_deps in ndeps_values) {
      # Perform the optimization with the current set of starting values and n_deps
      opt = optim(c(start_val_1,start_val_2), llfun, N=N1, Nref=Nref, lam=unlist(lam.v.cut), bstar=unlist(bstar1.v.cut), M=M.ref, lim=lim, method ='L-BFGS-B', lower=c(0,0), upper=c(1,20),control = list(factr = 1e-8, maxit = 1000,trace=0,ndeps = c(n_deps, n_deps*100),fnscale=-1))
     opt0 <- optim(start_val_2, llfun0, N=N1, Nref=Nref, lam=unlist(lam.v.cut), bstar=unlist(bstar1.v.cut), M=M.ref, lim=lim, method ='L-BFGS-B', lower=0, upper=20,control = list(factr = 1e-8, maxit = 1000,trace=0,ndeps = n_deps*100,fnscale=-1))

      # Check if this result is the best so far, based on value and convergence
      if (opt$convergence == 0 && opt0$convergence== 0 && opt$value > best_value && opt$value>opt0$value) {
        ll_alt.1 = best_value = opt$value
        h11.hdl.cut <- opt$par
        ll_alt.1conv <- opt$convergence

        ll_null.1 = opt0$value
        ll_null.1conv= opt0$convergence
        converged = TRUE
      }
    }
  }
}

error.message.h1 = error.message.h2 = error.message.rg = NA
if (!converged | h11.hdl.cut[1] < 0 | h11.hdl.cut[1] > 1) {
    error.message.h1 <- "Optimization of trait 1 heritability did not converge."
    if(output.file != ""){
      cat(error.message.h1, file = output.file, append = T)
    }
    cat(error.message.h1)
}
    ##Estimate heritabitlity for trait 2
    starting_values_1 <- c(h22.wls[2],0,0.5)
    starting_values_2 <- c(1,0.5,1.5, 0)
    best_value <- -Inf
    ndeps_values <- c(1e-5, 1e-8, 1e-16)
    converged = FALSE
for (start_val_1 in starting_values_1) {
    for (start_val_2 in starting_values_2) {
        for (n_deps in ndeps_values) {
            # Perform the optimization with the current set of starting values and n_deps
            opt = optim(c(start_val_1,start_val_2), llfun, N=N2, Nref=Nref, lam=unlist(lam.v.cut), bstar=unlist(bstar2.v.cut), M=M.ref, lim=lim, method ='L-BFGS-B', lower=c(0,0), upper=c(1,20),control = list(factr = 1e-8, maxit = 1000,trace=0,ndeps = c(n_deps, n_deps*100),fnscale=-1))
            opt0 <- optim(start_val_2, llfun0, N=N2, Nref=Nref, lam=unlist(lam.v.cut), bstar=unlist(bstar2.v.cut), M=M.ref, lim=lim, method ='L-BFGS-B', lower=0, upper=20,control = list(factr = 1e-8, maxit = 1000,trace=0,ndeps = n_deps*100,fnscale=-1))

            # Check if this result is the best so far, based on value and convergence
            if (opt$convergence == 0 && opt0$convergence== 0 && opt$value > best_value && opt$value>opt0$value) {
                ll_alt.2 = best_value = opt$value
                h22.hdl.cut <- opt$par
                ll_alt.2conv <- opt$convergence

                ll_null.2 = opt0$value
                ll_null.2conv= opt0$convergence

                converged = TRUE
            }
        }
    }
}
if (!converged | h22.hdl.cut[1] < 0 | h22.hdl.cut[1] > 1) {
    error.message.h2 <- "Optimization of trait 2 heritability did not converge."
    if(output.file != ""){
        cat(error.message.h2, file = output.file, append = T)
    }
    cat(error.message.h2)
}

      lam.v.use <- lam.v.cut
      bstar1.v.use <- bstar1.v.cut
      bstar2.v.use <- bstar2.v.cut

      h11.hdl.use <- h11.hdl.cut
      h22.hdl.use <- h22.hdl.cut

      h11 <- h11.hdl.use[1]
      h22 <- h22.hdl.use[1]

    ##Estimate genetic covariance

if (h11 <= 0 | h22 <= 0) {
   h12 = NA
   int = NA
   rg = NA
   rg.lower = rg.upper = NA
  ll_alt.h12 = NA
  ll_alt.h12conv = NA
  ll_null.h12 = NA
  ll_null.h12conv = NA
  h12.hdl.use <- NA
}else{
starting_values_1 <- c(h12.wls[2],0,-sqrt(h11*h22)*0.5, sqrt(h11*h22)*0.5)
starting_values_2 <- c(rho12, 1, 0)
best_value <- -Inf
ndeps_values <- c(1e-5, 1e-8, 1e-16)
converged = FALSE
for (start_val_1 in starting_values_1) {
  for (start_val_2 in starting_values_2) {
    for (n_deps in ndeps_values) {
      # Perform the optimization with the current set of starting values and n_deps
     opt <- optim(c(start_val_1, start_val_2), llfun.gcov.part.2, h11 = h11.hdl.use, h22 = h22.hdl.use,rho12 = rho12, M = M.ref, N1 = N1, N2 = N2, N0 = N0, Nref = Nref, lam1 = unlist(lam.v.use), lam2 = unlist(lam.v.use), bstar1 = unlist(bstar1.v.use), bstar2 = unlist(bstar2.v.use), lim = lim, method = 'L-BFGS-B', lower = c(-sqrt(h11* h22), -20), upper = c(sqrt(h11* h22), 20), control = list(factr = 1e-8, maxit = 1000, trace = 0, ndeps = c(n_deps, n_deps * 100), fnscale = -1))
      opt0 <- optim(start_val_2, llfun0.gcov.part.2, h11=h11.hdl.use, h22=h22.hdl.use, rho12=rho12, M=M.ref, N1=N1, N2=N2, N0=N0, Nref=Nref, lam1=unlist(lam.v.use), lam2=unlist(lam.v.use),bstar1=unlist(bstar1.v.use), bstar2=unlist(bstar2.v.use),lim=lim, method ='L-BFGS-B', lower=-20, upper=20,control=list(factr = 1e-8, maxit = 1000,trace=0,ndeps=n_deps*100,fnscale=-1))

      # Check if this result is the best so far, based on value and convergence
      if (opt$convergence == 0 && opt0$convergence== 0 && opt$value > best_value && opt$value>opt0$value) {
         ll_alt.h12 = best_value = opt$value
         h12.hdl.use <- opt$par
         h12 <- h12.hdl.use[1]
         int <- h12.hdl.use[2]
        ll_alt.h12conv <- opt$convergence

        ll_null.h12 = opt0$value
        ll_null.h12conv = opt0$convergence
        converged = TRUE
      }
    }
  }
}
if(converged == FALSE){
  ll_alt.h12 = NA
  ll_alt.h12conv = NA
  ll_null.h12 = NA
  ll_null.h12conv = NA
  h12.hdl.use <- NA
  h12 <- NA
  int <- NA
  error.message.rg <- "Optimization of genetic covariance did not converge. \n"
    if(output.file != ""){
            cat(error.message.rg, file = output.file, append = T)
    }
          stop(error.message.rg)
}

}
    ##Compute LRT p-value
    lr.1 <- -2*(ll_null.1 - ll_alt.1)
    p.h1.lrt <- pchisq(lr.1, 1, lower.tail = FALSE)

    lr.2 <- -2*(ll_null.2 - ll_alt.2)
    p.h2.lrt <- pchisq(lr.2, 1, lower.tail = FALSE)

    lr.12 <- -2*(ll_null.h12 - ll_alt.h12)
    p.h12.lrt <- pchisq(lr.12, 1, lower.tail = FALSE)

    ##Calculate the likelihood-based confidence interval of rg

    if(h11 <= 0){h11.LCI =0}else{h11.LCI = h11}
    if(h22 <= 0){h22.LCI =0}else{h22.LCI = h22}

    if(h11.LCI == 0 | h22.LCI == 0){
      h11 = h11.LCI
      h22 = h22.LCI
      rg = NA
      rg.lower = rg.upper = NA
    }else{
    h12_val = seq(-1*sqrt(h11.LCI* h22.LCI), 1*sqrt(h11.LCI* h22.LCI), length = 10000)

    # Apply the function over each row of the grid
    ll_values = sapply(h12_val, function(p) {
    llfun.gcov.part.2(c(p, int), h11=h11.hdl.use, h22=h22.hdl.use, rho12=rho12, M=M.ref, N1=N1, N2=N2, N0=N0, Nref=Nref, lam1=unlist(lam.v.use), lam2=unlist(lam.v.use), bstar1=unlist(bstar1.v.use), bstar2=unlist(bstar2.v.use), lim=lim)
    })
    h12.LCI = h12_val[which.max(ll_values)]
    likelihoods <- exp(ll_values - max(ll_values))
    ###1. get confidence interval by wilks statistic ratio
    alpha = alpha
    c <- exp(-qchisq(1 - alpha, 1) / 2)
    l_loc = which(likelihoods > c)
    h12_CI = h12_val[l_loc]
    h12.lower= min(h12_CI)
    h12.upper= max(h12_CI)

    rg = h12.LCI/sqrt(h11.LCI*h22.LCI)
    rg.lower = max(-1, h12.lower/sqrt(h11.LCI*h22.LCI))
    rg.upper = min(1, h12.upper/sqrt(h11.LCI*h22.LCI))
    }

    p_score <- case_when(
      is.na(p.h12.lrt)         ~ -999,  # Flag missing values
      p.h12.lrt < 0.0005        ~ 3,
      p.h12.lrt < 0.005         ~ 2,
      p.h12.lrt < 0.05         ~ 1,
      TRUE                     ~ 0
    )

    rg_score <- case_when(
      is.na(rg)                ~ -999,  # Flag missing values
      abs(rg) > 0.9                 ~ 3,
      abs(rg) > 0.7                 ~ 2,
      abs(rg) > 0.5                 ~ 1,
      TRUE                     ~ 0
    )

    total_score <- ifelse(p_score >= 0 & rg_score >= 0,
                          p_score + rg_score,
                          NA_integer_)


    evidence_for_colocalization <- case_when(
      is.na(p.h12.lrt) | is.na(rg)        ~ "NA",
      p.h12.lrt >= 0.05 | abs(rg) <= 0.5  ~ " ",
      total_score >= 6                    ~ "*****",
      total_score == 5                    ~ "****",
      total_score == 4                    ~ "***",
      total_score == 3                    ~ "**",
      total_score == 2                    ~ "*",
      TRUE                                ~ " "
    )


    if(intercept.output == T){
      h11.intercept <- h11.hdl.use[2]
      h22.intercept <- h22.hdl.use[2]
      h12.intercept <- h12.hdl.use[2]
    }

    output <- function(value){
      if(is.na(value)){
        value.out <- NA
      } else if(abs(value) < 1e-4){
        value.out <- formatC(value, format = "e", digits = 2)
      } else {
        value.out <- round(value, digits = 4)
      }
    }

    if(h11 == 0 | h22 == 0 ){
      cat("Warning: Heritability of one trait was estimated to be 0, which may due to:
          1) The true heritability is very small;
          2) The sample size is too small;
          3) Many SNPs in the chosen reference panel misses in the GWAS;
          4) There is a severe mismatch between the GWAS population and the population for computing reference panel")
    }
    cat("\n")

if(intercept.output == FALSE){
      estimates.df <- matrix(c(h11,p.h1.lrt,h22,p.h2.lrt, h12, rg, rg.lower, rg.upper, p.h12.lrt, evidence_for_colocalization), nrow=1,ncol=10)
      colnames(estimates.df) <- c("Heritability_1","P_value_Heritability_1","Heritability_2", "P_value_Heritability_2", "Genetic_Covariance","Genetic_Correlation","Lower_bound_rg","Upper_bound_rg", "P", "Evidence_for_colocalization")
    } else{
      estimates.df <- matrix(c(h11,p.h1.lrt,h22,p.h2.lrt, h12, rg, rg.lower, rg.upper, p.h12.lrt,evidence_for_colocalization, h11.intercept, NA, h22.intercept, NA, h12.intercept,NA,NA,NA,NA,NA), nrow=2,ncol=10)
      colnames(estimates.df) <- c("Heritability_1","P_value_Heritability_1", "Heritability_2", "P_value_Heritability_2", "Genetic_Covariance","Genetic_Correlation","Lower_bound_rg","Upper_bound_rg", "P", "Evidence_for_coloalization")
      rownames(estimates.df) <- c("Estimate", "Intercept")
    }

    end.time <- date()
   if(!is.na(error.message.h1)){
  cat(error.message.h1, "\n")
}
if(!is.na(error.message.h2)){
  cat(error.message.h2, "\n")
}
if(!is.na(error.message.rg)){
  cat(error.message.rg, "\n")
}
    cat("\n")
    cat("Estimates: \n")
    cat(paste0("Heritability of phenotype 1: ",
        output(h11),
        ", P = ", output(p.h1.lrt), " \n"))
    cat(paste0("Heritability of phenotype 2: ",
        output(h22),
        ", P = ", output(p.h2.lrt), " \n"))
    cat(paste0("Genetic Correlation: ",
        output(rg),
        ", 95% confidence interval (", output(rg.lower),", ",output(rg.upper), "), P = ", output(p.h12.lrt), " \n"))
    cat("Evidence for colocalization:",
        evidence_for_colocalization, " \n")
        cat("
Codes:
      |rG| 0 ' ' 0.5 '+' 0.7 '++' 0.9 '+++' 1
P
1
' '          ' '     ' '     ' '      ' '
0.05
'+'          ' '     '*'     '**'     '***'
0.005
'++'         ' '     '**'    '***'    '****'
0.0005
'+++'        ' '     '***'   '****'   '*****'
0  \n")
    if(h11 == 0 | h22 == 0 ){
      cat("Warning: Heritability of one trait was estimated to be 0, which may due to:
          1) The true heritability is very small;
          2) The sample size is too small;
          3) Many SNPs in the chosen reference panel misses in the GWAS.")
    }

    cat("\n")
    cat("Analysis finished at",end.time,"\n")

    if(output.file != ""){
      cat("\n", file = output.file, append = TRUE)
      cat("\n", file = output.file, append = TRUE)
    if(!is.na(error.message.h1)){
        cat(error.message.h1, file = output.file, append = TRUE, "\n")
    }
    if(!is.na(error.message.h2)){
    cat(error.message.h2, file = output.file, append = TRUE, "\n")
    }
    if(!is.na(error.message.rg)){
    cat(error.message.rg, file = output.file, append = TRUE, "\n")
    }

      cat("\n")
      cat("Estimates: \n")
      cat(paste0("Heritability of phenotype 1: ",
          output(h11),
          ", P = ", output(p.h1.lrt), " \n"), file = output.file, append = TRUE)
      cat(paste0("Heritability of phenotype 2: ",
          output(h22),
          ", P = ", output(p.h2.lrt), " \n"), file = output.file, append = TRUE)
      cat(paste0("Genetic Correlation: ",
          output(rg),
          ", 95% confidence interval (", output(rg.lower),", ",output(rg.upper), "), P = ", output(p.h12.lrt), " \n"), file = output.file, append = TRUE)
      cat("Evidence for colocalization:",
          evidence_for_colocalization, " \n", file = output.file, append = TRUE)
    cat("
Codes:
      |rG| 0 ' ' 0.5 '+' 0.7 '++' 0.9 '+++' 1
P
1
' '          ' '     ' '     ' '      ' '
0.05
'+'          ' '     '*'     '**'     '***'
0.005
'++'         ' '     '**'    '***'    '****'
0.0005
'+++'        ' '     '***'   '****'   '*****'
0  \n", file = output.file, append = TRUE)
      if(h11 == 0 | h22 == 0 ){
        cat("Warning: Heritability of one trait was estimated to be 0, which may due to:
            1) The true heritability is very low;
            2) The sample size of the GWAS is too small;
            3) Many SNPs in the chosen reference panel are absent in the GWAS;
            4) There is a severe mismatch between the GWAS population and the population for computing reference panel", file = output.file, append = TRUE)
      }



      cat("\n", file = output.file, append = TRUE)
      cat("Analysis finished at",end.time,"\n", file = output.file, append = TRUE)
      cat("The results were saved to", output.file)
      cat("\n")
      cat("The results were saved to", output.file, file = output.file, append = TRUE)
      cat("\n", file = output.file, append = TRUE)
    }
    result.df <- data.frame(Trait1 = Trait1name, Trait2 = Trait2name, chr = chr, piece = piece, eigen_use = eigen.use, estimates.df
  )
     return(result.df)
    }
