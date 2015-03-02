#########################################################################
# troisième groupe : fonctions utilitaires
#########################################################################
VB0 <- function(corp, attr, num=0){

  # création des objets nécessaires
  # aux analyses globales du vocabulaire.
  # si l'attribut est numérique,
  # mettre num=1

  library(rcqp, quietly=TRUE, warn.conflicts=FALSE)
  options(warn=-1)
  gcinfo(verbose=FALSE)

  cqi_query(corp, "Tt", '[lemma=".*"]')
  quer <- paste(corp, "Tt", sep=":")
  dumprov <- cqi_dump_subcorpus(quer)
  cps <- dumprov[,1]
  rm(dumprov)
  gc()

  cplem <- paste(corp, "lemma", sep=".")
  cpfrm <- paste(corp, "word", sep=".")
  cppos <- paste(corp, "pos", sep=".")
  cpattr <- paste(corp, attr, sep=".")

  tabcp <- as.data.frame(cps)
  tabcp[,2] <- cqi_cpos2id(cplem, cps)
  tabcp[,3] <- cqi_cpos2id(cpfrm, cps)
  tabcp[,4] <- cqi_cpos2id(cppos, cps)
  tabcp[,5] <- cqi_cpos2str(cplem, cps)
  tabcp[,6] <- cqi_cpos2str(cpfrm, cps)

  iddt <- cqi_cpos2struc(cpattr, cps)
  tabcp[,7] <- cqi_struc2str(cpattr, iddt)
  if (num == 1) {
    tabcp[,7] <- as.numeric(tabcp[,7])
  }
  names(tabcp) <- c("cpos", "idlem", "idfrm", "idpos", "lem", "frm", "attr")

  cat("\n", "Effectif (tokens) : ", nrow(tabcp), "\n")
  cat("Nombre de lemmes : ", length(unique(tabcp[,2])), "\n")
  cat("Nombre de formes : ", length(unique(tabcp[,3])), "\n\n")

  class(tabcp) <- "VBASE"
  return(tabcp)

}
##############################################
VBZ <- function(vbase, D="", F="", val="", hz="") {

  # en entrée un objet créé par VBase()
  # troncature éventuelle par les dates
  # et/ou un autre attribut
  # possibilité de tirage au hasard
  # création des objets >> zipfR
  # le paramètre hz est utilisé pour générer un tirage au hasard
  # selon une distribution de poisson ; pour obtenir un effectif
  # donné, procéder par essais successifs
  # en gros, la proportion souhaitée (essayer entre 0.1 et 1)

  library(rcqp, quietly=TRUE, warn.conflicts=FALSE)
  library(zipfR, quietly=TRUE, warn.conflicts=FALSE)
  options(warn=-1)
  gc()

  if (!inherits(vbase, "VBASE")) stop("en entrée : un objet produit par VB0()", call.=F)
  vbase <- unclass(vbase)
  vbase <- as.data.frame(vbase)

  if (D != "") {
    D <- as.numeric(D)
    vbase <- vbase[(vbase[,1]>=D),]
  }
  if (F != "") {
    F <- as.numeric(F)
    vbase <-vbase[(vbase[,1]<=F),]
  }
  if (val != "") {
    vbase <- vbase[(vbase[,7]==val),]
  }

  cposd <- vbase[1,1]
  cposf <- vbase[(nrow(vbase)),1]

  if (hz != "") {
    hz <- as.numeric(hz)
    rg <- cposf - cposd
    hzd <- round(rpois(rg, hz))
    vbase <- vbase[(hzd>0),]
  }

  Ltabtfl <- vec2tfl(vbase[,5])
  Ltabvgc <- vec2vgc(vbase[,5])
  Ltabspc <- vec2spc(vbase[,5])
  Ftabtfl <- vec2tfl(vbase[,6])
  Ftabvgc <- vec2vgc(vbase[,6])
  Ftabspc <- vec2spc(vbase[,6])

  cat("\n", "Nombre de tokens retenus : ", nrow(vbase))

  vb1 <- list(cposd = cposd, cposf=cposf, Ltfl=Ltabtfl, Lvgc=Ltabvgc, Lspc=Ltabspc, Ftfl=Ftabtfl, Fvgc=Ftabvgc, Fspc=Ftabspc)
  gc()
  return(vb1)
}
#####################################################
flou <- function(corp, type, larg=0, mot) {

  # fonction de recherche floue
  # choix du corpus, du type (lemma ou word)
  # et de la distance tolérée
  # NB : la méthode choisie n'est pas la plus simple,
  # pour éviter les débordements de mémoire
  # AG gpl3 avril 2013

  library(rcqp, quietly=TRUE, warn.conflicts=FALSE)
  options(warn=-1)
  atr <- cqi_attributes(corp, "p")
  if (!(type %in% atr)) stop ('paramètre "type" inutilisable !', call.=FALSE)

  ln <- size(corpus(corp))
  ln <- ln-1
  Lcpos <- (0:ln)
  cq <- paste(corp,".",type, sep="")
  Lid <- cqi_cpos2id(cq, Lcpos)
  Lstr <- cqi_id2str(cq, Lid)
  List <- unique(Lstr)
  fuzzy <- agrep(mot, List, max=larg, ignore.case=TRUE, value=TRUE)
  if (length(fuzzy)==0) stop ('pas de réponse, modifiez le paramètre "larg" !', call.=FALSE)


  # affichages différents pour word et lemma : 2 routines distinctes

  if (type == "word") {
    Leff <- Llm <- Lpos <- NULL
    cql <- paste(corp,".", "lemma", sep="")
    cqp <- paste(corp,".", "pos", sep="")
    for (i in 1:length(fuzzy)){
      Nid <- cqi_str2id(cq, fuzzy[i])
      Ncpos <- cqi_id2cpos(cq, Nid)
      Leff[i] <- as.numeric(length(Ncpos))
      Lpos[i] <- cqi_cpos2str(cqp, Ncpos[1])
      Llm[i] <- cqi_cpos2str(cql, Ncpos[1])
    }
    Ltt <- cbind(fuzzy, Leff, Lpos, Llm)
    Ltt <- Ltt[order(Ltt[,1]),]
    return (Ltt)
  }

  if (type == "lemma") {
    Leff <- Lpos <- NULL
    cqp <- paste(corp,".", "pos", sep="")
    for (i in 1:length(fuzzy)){
      Nid <- cqi_str2id(cq, fuzzy[i])
      Ncpos <- cqi_id2cpos(cq, Nid)
      Leff[i] <- as.numeric(length(Ncpos))
      Lpos[i] <- cqi_cpos2str(cqp, Ncpos[1])
    }
    Ltt <- cbind(fuzzy, Leff, Lpos)
    Ltt <- Ltt[order(Ltt[,1]),]
    return (Ltt)
  }
}
lisible <- function (x, y, lab,mn,mx, cex=.2){

  # on constitue le tab(leau de )dep(art)
  # library(circular, quietly=TRUE, warn.conflicts=FALSE)
  eps <- 0.0000000001
  tabdep <- as.data.frame(cbind(x,y,lab))
  names(tabdep) <- c("x","y","lab")
  row.names(tabdep) <- seq(1,nrow(tabdep))
  tabdep$x <- as.numeric(as.character(tabdep[,1]))
  tabdep$y <- as.numeric(as.character(tabdep[,2]))
  tabdep$lab <- as.character(tabdep$lab)
  htlet <- (mx-mn)/(30/cex)
  lglet <- htlet*.5
  H <- lglet/2
  indx <- as.numeric(row.names(tabdep))
  d2 <- (tabdep$x^2)+(tabdep$y^2)
  drt <- tabdep$x + (H*nchar(tabdep$lab))
  gau <- tabdep$x - (H*nchar(tabdep$lab))
  angl <- deg(atan(tabdep$y/tabdep$x))/.9

  tabdep <- as.data.frame(cbind(tabdep,indx,d2,drt,gau,angl))
  tt <- length(x)
  # sens <- 0
  tabfin <- tabpro <- tabdep
  rn <- (runif(tt*100))*100

  for (i in 1:tt){

    # on trie et on évacue la première ligne >> tableau final
    tabpro <- tabpro[sort.list(tabpro$d2),]
    cnt <- (tabpro[1,])
    tabfin[i,] <- cnt
    tabpro <- tabpro[-1,]

    # il faut repousser tout ce qui peut recouvrir le point actif  (cnt)
    # constitution du rub(an) formé de tous les points à écarter
    if (nrow(tabpro)==0) next
    cnt[1] <- as.numeric(as.character(cnt[1]))-(eps*sign(as.numeric(as.character(cnt[1]))))
    cnt[2] <- as.numeric(as.character(cnt[2]))-(eps*sign(as.numeric(as.character(cnt[2]))))
    ruban <- tabpro[(abs(as.numeric(tabpro$y)-as.numeric(as.character(cnt[2])))< htlet),]

    if (nrow(ruban) == 0) next
    rubg <- ruban[(ruban$x < as.numeric(as.character(cnt[1])) & ruban$drt > as.numeric(as.character(cnt[7]))),]
    rubd <- ruban[(ruban$x > as.numeric(as.character(cnt[1])) & ruban$gau < as.numeric(as.character(cnt[6]))),]
    rub <- rbind(rubg,rubd)
    rub <- unique(rub)
    if (nrow(rub) == 0) next
    n <- nrow(rub)
    r <- 1

    # on écarte tous les points du rub(an) alternativement horizontalement et verticalement, vers l'extérieur du quadrant
    # en combinant la valeur de l'angle et un nombre aléatoire (!)

    for (j in 1:n){
      if (rub[j,1]>0 & rub[j,2]>0 & rub[j,8]<rn[r]) tabpro[(tabpro[,4]==rub[j,4]),1] <- cnt[6]+(H*nchar(rub[j,3]))
      if (rub[j,1]>0 & rub[j,2]>0 & rub[j,8]>=rn[r]) tabpro[(tabpro[,4]==rub[j,4]),2] <- cnt[2]+(htlet)

      if (rub[j,1]>0 & rub[j,2]<0 & abs(rub[j,8])<rn[r]) tabpro[(tabpro[,4]==rub[j,4]),1] <- cnt[6]+(H*nchar(rub[j,3]))
      if (rub[j,1]>0 & rub[j,2]<0 & abs(rub[j,8])>=rn[r]) tabpro[(tabpro[,4]==rub[j,4]),2] <- cnt[2]-(htlet)

      if (rub[j,1]<0 & rub[j,2]<0 & rub[j,8]<rn[r]) tabpro[(tabpro[,4]==rub[j,4]),1] <- cnt[7]-(H*nchar(rub[j,3]))
      if (rub[j,1]<0 & rub[j,2]<0 & rub[j,8]>=rn[r]) tabpro[(tabpro[,4]==rub[j,4]),2] <- cnt[2]-(htlet)

      if (rub[j,1]<0 & rub[j,2]>0 & abs(rub[j,8])<rn[r]) tabpro[(tabpro[,4]==rub[j,4]),1] <- cnt[7]-(H*nchar(rub[j,3]))
      if (rub[j,1]<0 & rub[j,2]>0 & abs(rub[j,8])>=rn[r]) tabpro[(tabpro[,4]==rub[j,4]),2] <- cnt[2]+(htlet)
      r <- r+1
    }

    # on recalcule la position relative de tous les points restants
    # de manière à être sûr d'attaquer le bon point au tour suivant
    tabpro$d2 <- (tabpro$x^2) + (tabpro$y^2)
    tabpro$drt <- tabpro$x + (H*nchar(tabpro$lab))
    tabpro$gau <- tabpro$x - (H*nchar(tabpro$lab))
  }

  # on remet le tableau final dans l'ordre des lignes au départ (indx)
  tabfin <- tabfin[sort.list(tabfin$indx),]
  tabfin[,3] <- lab

  return(tabfin)
}
