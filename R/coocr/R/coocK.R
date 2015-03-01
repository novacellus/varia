###################################################################################

coocK2 <- function(corp, attr, lm1, lm2, dis=8, lim=200) {

  # affichage d'une ligne de texte brut (20 tokens)
  # pour une paire de lemmes (a - b  ou b - a)
  # avec indication d'un attribut

  library(plyr, quietly=TRUE, warn.conflicts=FALSE)
  library(rcqp, quietly=TRUE, warn.conflicts=FALSE)
  library(MASS, quietly=TRUE, warn.conflicts=FALSE)
  options(warn=-1)

  obj  <- paste(corp, ".lemma", sep="")
  obj2 <- paste(corp, ".word", sep="")
  obj3 <- paste(corp, attr, sep=".")
  a <- cqi_str2id(obj, lm1)
  b <- cqi_str2id(obj, lm2)
  if (a < 0) stop(paste(lm1, "n'existe pas !"), call.=FALSE)
  if (b < 0) stop(paste(lm2, "n'existe pas !"), call.=FALSE)
  listattr <- cqi_attributes(corp, "s")
  if ((attr != "") & ((attr %in% listattr)==FALSE)) {
    stop(paste(attr, "n'est pas un attribut de ce corpus !"), call.=FALSE)
  }
  LL1 <- "Lsub1"
  LL2 <- "Lsub2"
  dem1 <- paste("[lemma = '" , lm1 , "'][]{0," , dis , "}[lemma = '" , lm2 , "']", sep="")
  dem2 <- paste("[lemma = '" , lm2 , "'][]{0," , dis , "}[lemma = '" , lm1 , "']", sep="")

  cqi_query(corp, LL1, dem1)
  cqi_query(corp, LL2, dem2)
  sc1 <- paste(corp, LL1, sep = ":")
  sc2 <- paste(corp, LL2, sep = ":")
  if (cqi_subcorpus_size(sc1) > 0) {
    A1 <- cqi_dump_subcorpus(sc1)}
  else A1 <- NULL
  if (cqi_subcorpus_size(sc2) > 0) {
    A2 <- cqi_dump_subcorpus(sc2)}
  else A2 <- NULL
  AA <- rbind(A1, A2)
  if (is.null(AA)) stop("Aucun résultat !", call. = FALSE)
  if (nrow(AA)> 1) AA <- AA[order(AA[,1]),]
  res <- AA[,1]
  res <- unique(res)
  cat("\n","nombre d'occurrences :", length(res), "\n\n")


  if (attr != "") {
    res2 <- cqi_cpos2struc(obj3, res)
    res3 <- cqi_struc2str(obj3, res2)
  }
  if (attr == "") {
    res3 <- rep("", length(res))
  }

  d <- res - 6
  vec <- rep(0, 20*length(d))
  mat <- matrix(vec, nrow=length(d), ncol=20)
  for (i in 1:20) {
    mat[,i] <- d + i
  }

  ench2 <- function(x, obj2){
    a <- cqi_cpos2str(obj2, x)
    b <- paste(a[1:20], sep="", collapse=" ")
    return(b)
  }

  tabw <- apply(mat, 1, ench2, obj2)
  res4 <- cbind(as.matrix(tabw), as.matrix(res3))
  if (nrow(res4) <= lim) {
    write.matrix(res4, file="", sep=" ")
  }
  if (nrow(res4) > lim){
    tirage <- sort(as.integer(runif(lim,1,nrow(res4))))
    res5 <- res4[tirage,]
    write.matrix(res5, file="", sep=" ")
  }
  return(res4)
}

#############################################################################

coocK3 <- function(corp, attr, lm1, lm2, lm3, dis=8, lim=200) {

  # affichage d'une ligne de texte brut (24 tokens)
  # pour un triplet de lemmes (dans un ordre quelconque)
  # avec indication d'un attribut

  library(plyr, quietly=TRUE, warn.conflicts=FALSE)
  library(rcqp, quietly=TRUE, warn.conflicts=FALSE)
  library(MASS, quietly=TRUE, warn.conflicts=FALSE)
  options(warn=-1)

  obj  <- paste(corp, ".lemma", sep="")
  obj2 <- paste(corp, ".word", sep="")
  obj3 <- paste(corp, attr, sep=".")
  a <- cqi_str2id(obj, lm1)
  b <- cqi_str2id(obj, lm2)
  c <- cqi_str2id(obj, lm3)
  if (a < 0) stop(paste(lm1, "n'existe pas !"), call.=FALSE)
  if (b < 0) stop(paste(lm2, "n'existe pas !"), call.=FALSE)
  if (c < 0) stop(paste(lm3, "n'existe pas !"), call.=FALSE)
  listattr <- cqi_attributes(corp, "s")
  if ((attr != "") & ((attr %in% listattr)==FALSE)) {
    stop(paste(attr, "n'est pas un attribut de ce corpus !"), call.=FALSE)
  }

  LL1 <- "Lsub1"
  LL2 <- "Lsub2"
  LL3 <- "Lsub3"
  LL4 <- paste(LL1, LL2, sep="")
  LL5 <- paste(LL1, LL3, sep="")
  LL6 <- paste(LL2, LL3, sep="")
  dem1 <- paste("[lemma = '" , lm1 , "'][]{0," , dis , "}[lemma = '" , lm2 , "'][]{0," , dis , "}[lemma = '" , lm3 , "']", sep="")
  dem2 <- paste("[lemma = '" , lm1 , "'][]{0," , dis , "}[lemma = '" , lm3 , "'][]{0," , dis , "}[lemma = '" , lm2 , "']", sep="")
  dem3 <- paste("[lemma = '" , lm2 , "'][]{0," , dis , "}[lemma = '" , lm1 , "'][]{0," , dis , "}[lemma = '" , lm3 , "']", sep="")
  dem4 <- paste("[lemma = '" , lm2 , "'][]{0," , dis , "}[lemma = '" , lm3 , "'][]{0," , dis , "}[lemma = '" , lm1 , "']", sep="")
  dem5 <- paste("[lemma = '" , lm3 , "'][]{0," , dis , "}[lemma = '" , lm1 , "'][]{0," , dis , "}[lemma = '" , lm2 , "']", sep="")
  dem6 <- paste("[lemma = '" , lm3 , "'][]{0," , dis , "}[lemma = '" , lm2 , "'][]{0," , dis , "}[lemma = '" , lm1 , "']", sep="")

  cqi_query(corp, LL1, dem1)
  cqi_query(corp, LL2, dem2)
  cqi_query(corp, LL3, dem3)
  cqi_query(corp, LL4, dem4)
  cqi_query(corp, LL5, dem5)
  cqi_query(corp, LL6, dem6)

  sc1 <- paste(corp, LL1, sep = ":")
  sc2 <- paste(corp, LL2, sep = ":")
  sc3 <- paste(corp, LL3, sep = ":")
  sc4 <- paste(corp, LL4, sep = ":")
  sc5 <- paste(corp, LL5, sep = ":")
  sc6 <- paste(corp, LL6, sep = ":")
  if (cqi_subcorpus_size(sc1) > 0) {
    A1 <- cqi_dump_subcorpus(sc1)}
  else A1 <- NULL
  if (cqi_subcorpus_size(sc2) > 0) {
    A2 <- cqi_dump_subcorpus(sc2)}
  else A2 <- NULL
  if (cqi_subcorpus_size(sc3) > 0) {
    A3 <- cqi_dump_subcorpus(sc3)}
  else A3 <- NULL
  if (cqi_subcorpus_size(sc4) > 0) {
    A4 <- cqi_dump_subcorpus(sc4)}
  else A4 <- NULL
  if (cqi_subcorpus_size(sc5) > 0) {
    A5 <- cqi_dump_subcorpus(sc5)}
  else A5 <- NULL
  if (cqi_subcorpus_size(sc6) > 0) {
    A6 <- cqi_dump_subcorpus(sc6)}
  else A6 <- NULL
  AA <- rbind(A1, A2, A3, A4, A5, A6)
  if (is.null(AA)) stop("Aucun résultat !", call. = FALSE)
  if (nrow(AA) > 1)   AA <- AA[order(AA[,1]),]
  res <- AA[,1]
  res <- unique(res)
  cat("\n","nombre d'occurrences :", length(res), "\n\n")

  if (attr != "") {
    res2 <- cqi_cpos2struc(obj3, res)
    res3 <- cqi_struc2str(obj3, res2)
  }
  if (attr == "") {
    res3 <- rep("", length(res))
  }

  d <- res - 6
  vec <- rep(0, 20*length(d))
  mat <- matrix(vec, nrow=length(d), ncol=20)
  for (i in 1:20) {
    mat[,i] <- d + i
  }

  ench2 <- function(x, obj2){
    a <- cqi_cpos2str(obj2, x)
    b <- paste(a[1:20], sep="", collapse=" ")
    return(b)
  }

  tabw <- apply(mat, 1, ench2, obj2)
  res4 <- cbind(as.matrix(tabw), as.matrix(res3))
  if (nrow(res4) <= lim) {
    write.matrix(res4, file="", sep=" ")
  }
  if (nrow(res4) > lim){
    tirage <- sort(as.integer(runif(lim,1,nrow(res4))))
    res5 <- res4[tirage,]
    write.matrix(res5, file="", sep=" ")
  }
  return(res4)
}
