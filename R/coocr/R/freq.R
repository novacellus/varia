
##########################################################
# premier groupe : évolutions
##########################################################

freqlem <- function(corp, lm, trnch=20, attr="", val="", date) {

  # recherche d'un lemme et affichage
  # d'un histogramme simple (par tranches de mêmes effectifs)
  # à utiliser pour une première observation globale
  # A GUERREAU gpl3 janvier 2013

  library(MASS, quietly=TRUE, warn.conflicts=FALSE)
  library(rcqp, quietly=TRUE, warn.conflicts=FALSE)
  options(warn=-1)
  LL <- "Lsub"
  corpdate <- paste(corp, date, sep=".")

  # vérifications élémentaires
  if ((attr!="" & val=="") | (attr=="" & val!="")) stop ("Requête incomplète !", call.=FALSE)

  provl <- paste(corp, ".lemma", sep="")
  provp <- paste(corp, ".pos", sep="")
  idprov <- cqi_regex2id(provl, lm)

  if (length(idprov)==0){
    fl <- 1
    lm2 <- paste(lm, ".*", sep="")
    idprov2 <- cqi_regex2id(provl, lm2)

    if (length(idprov2) == 0) stop(paste(lm, "n'existe pas !"), call. = FALSE)

    strprov <- cqi_id2str(provl, idprov2)
    strcposprov <- NULL
    posprov <- NULL
    freq <- NULL
    cat("\n")
    for (i in 1:length(idprov2)){
      strcposprov[i] <- cqi_id2cpos(provl, idprov2[i])
      posprov[i] <- cqi_cpos2str(provp, strcposprov[i])
      freq[i] <- cqi_id2freq(provl, idprov2[i])
      cat(strprov[i], "/", posprov[i], "/", freq[i], "\n")
    }
    cat("\n")
    if (fl > 0) stop ("Choisissez parmi ces lemmes !", call.=FALSE)
  }

  # cas simple : corpus entier
  if (attr=="" & val=="") {
    quer1 <- paste ('[lemma="', lm, '" %cd]', sep ="")
    cqi_query(corp, LL, quer1)
    sc <- paste(corp, LL, sep = ":")
    A1d <- cqi_dump_subcorpus(sc)
    ef <- nrow(A1d)
    cat("\n")
    cat(lm, "effectif total : ", ef)
    flist <- cqi_fdist1(sc, "match", "word") # formes correspondant au lemme
    sc2 <- paste(corp, ".word", sep="")
    flistb <- data.frame(cqi_id2str(sc2, flist[,1]), flist[,2])
    cat("\n", "\n")
    write.matrix(flistb)

    # découpage et affichage de l'histogramme
    ln <- size(corpus(corp))
    xl <- paste(corp, " / ", trnch, "  (effectif total = ", ln, ")", sep = "")
    mn <- paste("lemme : ", lm, " (effectif du lemme = ", ef, ")", sep = "")
    br <- seq(1, ln, length.out=trnch+1)
    decal <- (ln/trnch)/2
    br2 <- as.integer((br[1:(length(br))-1]) + decal)
    listdate <- cqi_cpos2struc(corpdate, br2)
    listdateb <- cqi_struc2str(corpdate, listdate)
    opar <- par(mar=par("mar"))
    on.exit(par(opar))
    par(mar=c(4.2,4.2,2,0.5))
    H <- hist(A1d[,1], breaks = br, border = "grey", main = "", xlab = "", ylab = "Effectif par tranche")
    axis(side=1, at=br2, labels=listdateb, line=2)
    title(main = mn, line=1, cex.main=.8, font.main=1, adj = 0)
    mtext(xl, 3, line=0,cex=.8, font=1, adj=0)
  }

  # cas complexe : seulement une partie du corpus
  # nécessite de créer un sous-corpus provisoire complet correspondant
  # au critère de choix, pour pouvoir créer des tranches égales
  # dans ce sous-corpus
  if (attr!="" & val!="") {
    cat("Le calcul de l'effectif d'une fraction du corpus total prend du temps... \n")
    LL2 <- paste(LL, "2", sep="")
    quer2 <- paste('[lemma="',lm,'" %cd]', "::match.", attr, "=\"", val, "\"", sep="")
    cqi_query(corp, LL2, quer2)
    sc2 <- paste(corp, LL2, sep = ":")
    lsc <- cqi_subcorpus_size(sc2)
    if (lsc==0) stop("Pas de réponse (valeur de l'attribut inexistante ?), reformulez la requête !", call.=FALSE)
    A1d2 <- cqi_dump_subcorpus(sc2)
    ef <- nrow(A1d2)
    indexL <- as.integer(seq(1, ef, 1))
    efI <- cbind(as.integer(A1d2[,1]), indexL)
    LL3 <- paste(LL, "3", sep="")
    quertt <- paste("[lemma=\".*\"]", "::match.", attr, "=\"", val, "\"", sep="")
    cqi_query(corp, LL3, quertt)
    sc3 <- paste(corp, LL3, sep = ":")
    A2d3 <- cqi_dump_subcorpus(sc3)
    ln2 <- nrow(A2d3)

    # ici, une ruse diabolique !
    # (nécessaire pour contourner la discontinuité des cpos)
    indexT <- as.integer(seq(1, ln2, 1)) # index auxiliaire
    corpI <- cbind(as.integer(A2d3[,1]), indexT)
    choix <- corpI[,1] %in% efI[,1]
    corpI2 <- cbind(corpI, choix)
    indbon <- corpI2[(corpI2[,3]==TRUE),2]

    # découpage et affichage de l'histogramme
    mn <- paste("lemme : ", lm, " (effectif du lemme = ", ef, ")", sep = "")
    sb <- paste(attr, " = ", val, " (effectif = ", ln2, ")", sep="")
    br <- seq(1, ln2, length.out=trnch+1)
    decal <- (ln2/trnch)/2
    br2 <- as.integer((br[1:(length(br))-1]) + decal)
    br3 <- corpI2[br2,1]
    listdate <- cqi_cpos2struc(corpdate, br3)
    listdateb <- cqi_struc2str(corpdate, listdate)
    opar <- par(mar=par("mar"))
    on.exit(par(opar))
    par(mar=c(4.2,4.2,2,0.5))
    H <- hist(indbon, breaks = br, border = "grey", main = "", xlab = "", ylab = "Effectif par tranche")
    axis(side=1, at=br2, labels=listdateb, line=2)
    title(main = mn, line=1, cex.main=.8, font.main=1, adj = 0)
    mtext(sb, 3, line=0,cex=.8, font=1, adj=0)
  }

  # calcul et affichage d'une courbe de tendance (fonction lowess())
  mov <- lowess(H$counts)
  br2 <- (br[1:(length(br))-1]) + decal
  lines(br2, mov$y, col="blue", lwd=1.5)
  EFTR <<- H$counts
  LB1 <<- br2
}

#####################################################################################

freqlem2 <- function(corp, lm, trnch=20, attr="", val="", date, D=0, F=Inf) {

  # recherche d'un lemme et affichage
  # d'un histogramme simple (par tranches de mêmes effectifs)
  # à utiliser pour une première observation globale
  # affichage *chronologique* (anamorphose des tranches)
  # A GUERREAU gpl3 janvier 2013

  library(MASS, quietly=TRUE, warn.conflicts=FALSE)
  library(rcqp, quietly=TRUE, warn.conflicts=FALSE)
  options(warn=-1)
  LL <- "Lsub"

  # vérifications
  at <- cqi_attributes(corp, "s")
  if ((date %in% at) == FALSE) stop ('Indiquez *explicitement* le nom du paramètre "date" !', call.=FALSE)
  corpdate <- paste(corp, date, sep=".")
  if ((attr!="" & val=="") | (attr=="" & val!="")) stop ("Requête incomplète !", call.=FALSE)
  ln <- size(corpus(corp))
  if (F=="Inf") {
    F <- ln-1
  }
  lgcorp <- F - D # effectif du corpus traité

  provl <- paste(corp, ".lemma", sep="")
  provp <- paste(corp, ".pos", sep="")
  idprov <- cqi_regex2id(provl, lm)

  if (length(idprov)==0){
    fl <- 1
    lm2 <- paste(lm, ".*", sep="")
    idprov2 <- cqi_regex2id(provl, lm2)

    if (length(idprov2) == 0) stop(paste(lm, "n'existe pas !"), call. = FALSE)

    strprov <- cqi_id2str(provl, idprov2)
    strcposprov <- NULL
    posprov <- NULL
    freq <- NULL
    cat("\n")
    for (i in 1:length(idprov2)){
      strcposprov[i] <- cqi_id2cpos(provl, idprov2[i])
      posprov[i] <- cqi_cpos2str(provp, strcposprov[i])
      freq[i] <- cqi_id2freq(provl, idprov2[i])
      cat(strprov[i], "/", posprov[i], "/", freq[i], "\n")
    }
    cat("\n")
    if (fl > 0) stop ("Choisissez parmi ces lemmes !", call.=FALSE)
  }

  # cas simple
  if (attr=="" & val=="") {
    quer1 <- paste ('[lemma="', lm, '" %cd]', sep ="")
    cqi_query(corp, LL, quer1)
    sc <- paste(corp, LL, sep = ":")
    A1d <- cqi_dump_subcorpus(sc)

    # choix d'un sous-corpus en fonction des bornes D / F
    cposLm <- A1d[,1]
    cposLm <- cposLm[cposLm > D & cposLm < F]
    cposLm <- as.integer(cposLm)
    ef <- length(cposLm) # effectif du lemme

    cat("\n")
    cat(lm, " : effectif total dans le (sous)corpus considéré = ", ef)

    flist <- cqi_fdist1(sc, "match", "word") # formes correspondant au lemme
    sc2 <- paste(corp, ".word", sep="")
    lid <- cqi_cpos2id(sc2, cposLm)
    forme <- cqi_id2str(sc2, lid)
    flistb <- data.frame(lid,forme)
    lifreq <- as.data.frame(aggregate(1:nrow(flistb), flistb, length))
    lifreq <- lifreq[rev(order(lifreq[,3])),]
    cat("\n", "\n")
    write.matrix(lifreq[,2:3])

    xl <- paste(corp, " / ", trnch, "  (effectif total dans le (sous)corpus considéré  = ", lgcorp, ")", sep = "")
    mn <- paste("lemme : ", lm, " (effectif du lemme = ", ef, ")", sep = "")

    br <- seq(D, F, length.out=trnch+1)
    decal <- (lgcorp/trnch)/2
    br2 <- as.integer((br[1:(length(br))-1]) + decal)
    listdate2 <- cqi_cpos2struc(corpdate, br2)
    listdateb2 <- cqi_struc2str(corpdate, listdate2)

    # il faut essayer d'élininer tout ce qui n'est pas daté
    listdate <- cqi_cpos2struc(corpdate, br)
    listdateb <- cqi_struc2str(corpdate, listdate)
    listdateb <- as.numeric(listdateb)
    listdateb <- listdateb[listdateb < 9000]
    if (nchar(tail(listdateb,1))==3) {
      listdateb <- listdateb[listdateb < 999]
    }
    trnch2 <- length(listdateb)-1
    listdateb2 <- listdateb2[1:trnch2]

    opar <- par(mar=par("mar"))
    on.exit(par(opar))
    par(mar=c(4.2,4.2,2,0.5))
    H <- hist(cposLm, breaks = br, plot = FALSE)
    efftrn <- H$counts
    efftrn <- efftrn[1:trnch2]
    EFTR <<- efftrn
    LB1 <<- listdateb
    LB2 <<- listdateb2
    eff2 <- c(efftrn,1)
    plot(listdateb, eff2, type="n", xlab="chronologie", ylab="effectif par tranche")
    for (i in 1:trnch2) {
      #i <- 1:trnch2
      rect(listdateb[i], 0, listdateb[i+1], efftrn[i], border="grey")
    }

    # axis(side=1, at=br2, labels=listdateb, line=2)
    title(main = mn, line=1, cex.main=.8, font.main=1, adj = 0)
    mtext(xl, 3, line=0,cex=.8, font=1, adj=0)
  }

  # cas complexe
  if (attr!="" & val!="") {
    cat("Le calcul de l'effectif d'une fraction du corpus total prend du temps... \n")
    LL2 <- paste(LL, "2", sep="")
    quer2 <- paste('[lemma="',lm,'" %cd]', "::match.", attr, "=\"", val, "\"", sep="")
    cqi_query(corp, LL2, quer2)
    sc2 <- paste(corp, LL2, sep = ":")
    lsc <- cqi_subcorpus_size(sc2)
    if (lsc==0) stop("Pas de réponse (valeur de l'attribut inexistante ?), reformulez la requête !", call.=FALSE)
    A1d2 <- cqi_dump_subcorpus(sc2)

    # choix d'un sous-corpus en fonction des bornes D / F
    cposLm <- A1d2[,1]
    cposLm <- cposLm[cposLm > D & cposLm < F]
    cposLm <- as.integer(cposLm)
    ef <- length(cposLm) # effectif du lemme
    cat("\n")
    cat(lm, " : effectif total dans le (sous)corpus considéré = ", ef)
    cat("\n")

    indexL <- as.integer(seq(1, ef, 1))
    efI <- cbind(cposLm, indexL)
    LL3 <- paste(LL, "3", sep="")
    quertt <- paste("[lemma=\".*\"]", "::match.", attr, "=\"", val, "\"", sep="")
    cqi_query(corp, LL3, quertt)
    sc3 <- paste(corp, LL3, sep = ":")
    A2d3 <- cqi_dump_subcorpus(sc3)

    # restriction  D / F
    cposSc <- A2d3[,1]
    cposSc <- cposSc[cposSc > D & cposSc < F]
    cposSc <- as.integer(cposSc)
    ln2 <- length(cposSc)

    # ici, une ruse diabolique !
    indexT <- as.integer(seq(1, ln2, 1))
    corpI <- cbind(as.integer(cposSc), indexT)
    choix <- corpI[,1] %in% efI[,1]
    corpI2 <- cbind(corpI, choix)
    indbon <- corpI2[(corpI2[,3]==TRUE),2]

    flist <- cqi_fdist1(sc2, "match", "word")
    sc4 <- paste(corp, ".word", sep="")
    lid <- cqi_cpos2id(sc4, cposLm)
    forme <- cqi_id2str(sc4, lid)
    flistb <- data.frame(lid,forme)
    lifreq <- as.data.frame(aggregate(1:nrow(flistb), flistb, length))
    lifreq <- lifreq[rev(order(lifreq[,3])),]
    cat("\n", "\n")
    write.matrix(lifreq[,2:3])

    mn <- paste("lemme : ", lm, " (effectif du lemme = ", ef, ")", sep = "")
    sb <- paste(attr, " = ", val, " (effectif = ", ln2, ")", sep="")

    br <- seq(1, ln2, length.out=trnch+1)
    decal <- (ln2/trnch)/2
    br2 <- as.integer((br[1:(length(br))-1]) + decal)
    br3 <- corpI2[br2,1]
    listdate <- cqi_cpos2struc(corpdate, br3)
    listdateb <- cqi_struc2str(corpdate, listdate)
    listdateb2 <- listdateb[listdateb < 9000]
    if (nchar(tail(listdateb,1))==3) {
      listdateb <- listdateb[listdateb < 999]
    }
    trnch2 <- length(listdateb)-1
    listdateb2 <- listdateb2[1:trnch2]

    opar <- par(mar=par("mar"))
    on.exit(par(opar))
    par(mar=c(4.2,4.2,2,0.5))
    H <- hist(indbon, breaks = br,plot=FALSE)
    efftrn <- H$counts
    efftrn <- efftrn[1:trnch2]
    EFTR <<- efftrn
    LB1 <<- listdateb
    LB2 <<- listdateb2
    eff2 <- c(efftrn,1)
    plot(listdateb, eff2, type="n", xlab="chronologie", ylab="effectif par tranche")
    for (i in 1:trnch2) {
      rect(listdateb[i], 0, listdateb[i+1], efftrn[i], border="grey")
    }

    axis(side=1, at=br2, labels=listdateb, line=2)
    title(main = mn, line=1, cex.main=.8, font.main=1, adj = 0)
    mtext(sb, 3, line=0,cex=.8, font=1, adj=0)
  }

  # courbe de tendance
  mov <- lowess(efftrn)
  lines(listdateb2, mov$y, col="blue", lwd=1.5)
}

#####################################################################################

freqfrm <- function(corp, frm, trnch=20, attr="", val="") {

  # recherche d'une forme et affichage
  # d'un histogramme simple (par tranches de mêmes effectifs)
  # à utiliser pour une première observation globale

  library(MASS, quietly=TRUE, warn.conflicts=FALSE)
  library(rcqp, quietly=TRUE, warn.conflicts=FALSE)
  options(warn=-1)
  LL <- "Lsub"

  # vérifications
  if ((attr!="" & val=="") | (attr=="" & val!="")) stop ("Requête incomplète !", call.=FALSE)

  provl <- paste(corp, ".word", sep="")
  provp <- paste(corp, ".pos", sep="")
  provlm <- paste(corp, ".lemma", sep="")
  idprov <- cqi_regex2id(provl, frm)

  if (length(idprov)==0){
    fl <- 1
    lm2 <- paste(frm, ".*", sep="")
    idprov2 <- cqi_regex2id(provl, lm2)

    if (length(idprov2) == 0) stop(paste(frm, "n'existe pas !"), call. = FALSE)

    strprov <- cqi_id2str(provl, idprov2)
    strcposprov <- NULL
    posprov <- NULL
    freq <- NULL
    cat("\n")
    for (i in 1:length(idprov2)){
      strcposprov[i] <- cqi_id2cpos(provl, idprov2[i])
      posprov[i] <- cqi_cpos2str(provp, strcposprov[i])
      freq[i] <- cqi_id2freq(provl, idprov2[i])
      cat(strprov[i], "/", posprov[i], "/", freq[i], "\n")
    }
    cat("\n")
    if (fl > 0) stop ("Choisissez parmi ces lemmes !", call.=FALSE)
  }

  # cas simple
  if (attr=="" & val=="") {
    quer1 <- paste ("[word='", frm, "'%cd]", sep ="")
    cqi_query(corp, LL, quer1)
    sc <- paste(corp, LL, sep = ":")
    A1d <- cqi_dump_subcorpus(sc)
    ef <- nrow(A1d)
    cat("\n")
    cat(frm, "effectif total : ", ef)

    stopos <- cqi_cpos2str(provp, A1d)
    stolm <- cqi_cpos2str(provlm, A1d)
    stodf <- as.data.frame(cbind(stopos, stolm))
    stodf[,1] <- as.character(stodf[,1])
    stodf[,2] <- as.character(stodf[,2])
    stob <- as.data.frame(aggregate(1:nrow(stodf), stodf, length))
    cat("\n", unique(stob[,1]), unique(stob[,2]))
    flist <- cqi_fdist1(sc, "match", "word")
    sc2 <- paste(corp, ".word", sep="")
    flistb <- data.frame(cqi_id2str(sc2, flist[,1]), flist[,2])
    names(flistb) <- c("forme", "effectif")
    cat("\n", "\n")
    write.matrix(flistb)

    ln <- size(corpus(corp))
    xl <- paste(corp, " / ", trnch, "  (effectif total = ", ln, ")", sep = "")
    mn <- paste("forme : ", frm, " (effectif de la forme = ", ef, ")", sep = "")
    br <- seq(1, ln, length.out=trnch+1)
    decal <- (ln/trnch)/2
    H <- hist(A1d[,1], breaks = br, border = "grey", main = mn, xlab = xl, ylab = "Effectif par tranche")
  }

  # cas complexe
  if (attr!="" & val!="") {
    cat("Le calcul de l'effectif d'une fraction du corpus total prend du temps... \n")
    LL2 <- paste(LL, "2", sep="")
    quer2 <- paste("[word=\"",frm,"\"%cd]", "::match.", attr, "=\"", val, "\"", sep="")
    cqi_query(corp, LL2, quer2)
    sc2 <- paste(corp, LL2, sep = ":")
    lsc <- cqi_subcorpus_size(sc2)
    if (lsc==0) stop("Pas de réponse (valeur de l'attribut inexistante ?), reformulez la requête !", call.=FALSE)
    A1d2 <- cqi_dump_subcorpus(sc2)
    ef <- nrow(A1d2)
    cat("\n")
    cat(frm, "effectif total : ", ef)
    stopos <- cqi_cpos2str(provp, A1d2)
    stolm <- cqi_cpos2str(provlm, A1d2)
    stodf <- as.data.frame(cbind(stopos, stolm))
    stodf[,1] <- as.character(stodf[,1])
    stodf[,2] <- as.character(stodf[,2])
    stob <- as.data.frame(aggregate(1:nrow(stodf), stodf, length))
    cat("\n", unique(stob[,1]), unique(stob[,2]))


    flist <- cqi_fdist1(sc2, "match", "word")
    sc3 <- paste(corp, ".word", sep="")
    flistb <- data.frame(cqi_id2str(sc3, flist[,1]), flist[,2])
    names(flistb) <- c("forme", "effectif")
    cat("\n", "\n")
    write.matrix(flistb)

    cat("\n")
    indexL <- as.integer(seq(1, ef, 1))
    efI <- cbind(as.integer(A1d2[,1]), indexL)
    LL3 <- paste(LL, "3", sep="")
    quertt <- paste("[word=\".*\"]", "::match.", attr, "=\"", val, "\"", sep="")
    cqi_query(corp, LL3, quertt)
    sc3 <- paste(corp, LL3, sep = ":")
    A2d3 <- cqi_dump_subcorpus(sc3)
    ln2 <- nrow(A2d3)
    # ici, une ruse diabolique !
    indexT <- as.integer(seq(1, ln2, 1))
    corpI <- cbind(as.integer(A2d3[,1]), indexT)
    choix <- corpI[,1] %in% efI[,1]
    corpI2 <- cbind(corpI, choix)
    indbon <- corpI2[(corpI2[,3]==TRUE),2]
    mn <- paste("forme : ", frm, " (effectif de la forme = ", ef, ")", sep = "")
    sb <- paste("index du corpus cible (effectif total = ", ln2, ")", sep="")
    br <- seq(1, ln2, length.out=trnch+1)
    decal <- (ln2/trnch)/2
    H <- hist(indbon, breaks = br, border = "grey", main = mn, xlab = sb, ylab = "Effectif par tranche")
  }

  # courbe de tendance
  mov <- lowess(H$counts)
  br2 <- (br[1:(length(br))-1]) + decal
  lines(br2, mov$y, col="blue", lwd=1.5)
  EFTR <<- H$counts
  LB1 <<- br2
}

#####################################################################################

freqexp <- function(corp, exp, trnch=25, attr="", lim=200) {

  # recherche d'une expression (CQL) et affichage
  # d'un histogramme simple (par tranches de mêmes effectifs)
  # enregistrement de la liste des occurrences
  # affichage réduit à 'lim' occurrences max

  library(MASS, quietly=TRUE, warn.conflicts=FALSE)
  library(rcqp, quietly=TRUE, warn.conflicts=FALSE)
  options(warn=-1)
  if (exp=="") {
    stop("pas d'expression !")
  }

  # calculs et affichage de l'histogramme
  LL <- "EXP"
  provl <- paste(corp, ".word", sep="")
  provp <- paste(corp, ".pos", sep="")
  provlm <- paste(corp, ".lemma", sep="")
  obj3 <- paste(corp, attr, sep=".")
  cqi_query(corp, LL, exp)
  sc <- paste(corp, LL, sep = ":")
  A1d <- cqi_dump_subcorpus(sc)
  res <- A1d[,1]
  ef <- nrow(A1d)
  cat("\n")
  cat(exp, "effectif total : ", ef)
  cat("\n")

  ln <- size(corpus(corp))
  xl <- paste(corp, " / ", trnch, "  (effectif total = ", ln, ")", sep = "")
  mn <- paste("expression : ", exp, " (effectif de l'expression = ", ef, ")", sep = "")
  br <- seq(1, ln, length.out=trnch+1)
  decal <- (ln/trnch)/2
  H <- hist(A1d[,1], breaks = br, border = "grey", main = mn, xlab = xl, ylab = "Effectif par tranche")
  mov <- lowess(H$counts)
  br2 <- (br[1:(length(br))-1]) + decal
  lines(br2, mov$y, col="blue", lwd=1.5)
  EFTR <<- H$counts
  LB1 <<- br2

  # recherche et affichage de l'attribut choisi
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

  ench2 <- function(x, provl){
    a <- cqi_cpos2str(provl, x)
    b <- paste(a[1:20], sep="", collapse=" ")
    return(b)
  }

  tabw <- apply(mat, 1, ench2, provl)
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

##########################################################

freqcooc <- function(corp, lm1, lm2, dis=5, coeff="dice", trnch=10, D=0, F=Inf, date, selec=1) {

  # juillet 2014 Alain GUERREAU -> ||GPL3||
  # calculs des effectifs de lm1 et lm2 et de la cooc,
  # division en tranches égales ;
  # graphique lm1 lm2 coeff_cooc
  # + courbes de tendance
  # fenêtre = dis*2

  library(MASS, quietly=TRUE, warn.conflicts=FALSE)
  library(plyr, quietly=TRUE, warn.conflicts=FALSE)
  library(rcqp, quietly=TRUE, warn.conflicts=FALSE)
  options(warn=-1)
  gcinfo(verbose=FALSE)

  # AGENDA éliminer dates 9999 / tri par attr + val


  # vérifications préalables

  if (dis < 1) {
    stop("si dis=0, pas de cooccurrent !", call.=FALSE)
  }
  vcoeff <- c("dice","poiss","daille","pmi","hyperg","ms")
  if (coeff %in% vcoeff == FALSE) stop ("coefficient inconnu, choisir : dice/poiss/daille/pmi/hyperg/ms", call.=FALSE)

  provl <- paste(corp, ".lemma", sep="")
  a <- cqi_str2id(provl, lm1)
  if (a < 0) stop(paste(lm1, "n'existe pas !"), call.=FALSE)
  b <- cqi_str2id(provl, lm2)
  if (b < 0) stop(paste(lm2, "n'existe pas !"), call.=FALSE)
  # if ((attr!="" & val=="") | (attr=="" & val!="")) stop ("Requête incomplète !", call.=FALSE)
  file <- lm
  LL1 <- "Lsub1"
  LL2 <- "Lsub2"
  LL3 <- "Lsub3"
  corpdate <- paste(corp, date, sep=".")

  N <- size(corpus(corp))
  Ntr <- N / trnch
  br <- seq.int(1, N, length.out=trnch+1)
  decal <- (N/trnch)/2
  br2 <- as.integer((br[1:(length(br))-1]) + decal)

  listdate <- cqi_cpos2struc(corpdate, br2)
  listdateb <- cqi_struc2str(corpdate, listdate)

  # cas simple
  #if (attr=="" & val=="") {

  # calculs des effectifs observés
  quer1 <- paste ('[lemma="', lm1, '" %cd]', sep ="")
  quer2 <- paste ('[lemma="', lm2, '" %cd]', sep ="")
  quer3 <- paste ('([lemma="', lm1, '" %cd][]{0,',dis,'}[lemma="', lm2, '" %cd])|([lemma="', lm2, '" %cd][]{0,',dis,'}[lemma="', lm1, '" %cd])', sep ="")
  sc1 <- paste(corp, LL1, sep = ":")
  sc2 <- paste(corp, LL2, sep = ":")
  sc3 <- paste(corp, LL3, sep = ":")

  cqi_query(corp, LL1, quer1)
  dump1 <- cqi_dump_subcorpus(sc1)[,1]
  cqi_query(corp, LL2, quer2)
  dump2 <- cqi_dump_subcorpus(sc2)[,1]
  cqi_query(corp, LL3, quer3)
  subsz <- cqi_subcorpus_size(sc3)
  if (subsz == 0) stop ("Aucune cooccurrence !", call.=FALSE)
  dump3 <- cqi_dump_subcorpus(sc3)[,1]

  ef1 <- length(dump1)
  ef2 <- length(dump2)
  ef3 <- length(dump3)
  ef311 <- as.numeric(ef1) * as.numeric(ef2) / N

  cat("\n")
  cat(lm1, ", effectif total : ", ef1, "\n")
  cat(lm2, ", effectif total : ", ef2, "\n")
  cat("cooc ",lm1, "/",lm2,  ", effectif total : ", ef3, " (indépendance =", ef311,")", "\n\n")
  mn <- paste("corpus : ", corp, " (effectif :", N, ") /", trnch,"    ",lm1, " = ", ef1,"  ", lm2, " = ", ef2, "  cooc = ",ef3,sep = "")
  # sb <- paste(lm1, " = ", ef1,"    ", lm2, " = ", ef2, "    cooc = ",ef3, sep="")

  R1 <- hist(dump1, breaks=br, plot=FALSE)$counts
  C1 <- hist(dump2, breaks=br, plot=FALSE)$counts
  O11 <- hist(dump3, breaks=br, plot=FALSE)$counts
  RCO <- cbind(R1,C1,O11)
  E11 <- R1*C1/Ntr
  O11mE11 <- O11-E11
  mindiff <- min(O11mE11)
  if (mindiff < 0) cat ("tranche(s) sans cooccurrence ! \n\n")

  # calculs des coefficients
  c.dice <- (O11*1e+10)/(C1+R1)^selec
  c.poiss <- O11*(log(O11)-log(E11)-1)
  c.daille <- log((O11^3)/E11)
  c.pmi <- log(O11/E11)
  c.hyperg <- -phyper(O11, C1, (N-C1), (2*dis*R1), lower.tail=FALSE, log.p=TRUE)
  Or <- O11/R1
  Oc <- O11/C1
  OCR <- cbind(Or,Oc)
  c.ms <- apply(OCR,1,min)

  HDD <- cbind(RCO,c.dice,c.poiss,c.daille,c.pmi,c.hyperg,c.ms,O11, E11)
  cat(N,Ntr,ef1,ef2,ef3,ef311,sum(E11),"\n")

  # transformation en indices
  R1 <- R1/mean(R1)*100
  C1 <- C1/mean(C1)*100
  O11 <- O11/mean(O11)*100
  c.dice <- c.dice/mean(c.dice)*100
  c.poiss <- c.poiss/mean(abs(c.poiss))*100
  c.daille <- c.daille/mean(c.daille)*100
  c.pmi <- c.pmi/mean(c.pmi)*100
  c.hyperg <- c.hyperg/mean(c.hyperg)*100
  c.ms <- c.ms/mean(c.ms)*100

  # choix du coefficient
  if (coeff=="dice") c.coeff <- c.dice
  if (coeff=="poiss") c.coeff <- c.poiss
  if (coeff=="daille") c.coeff <- c.daille
  if (coeff=="pmi") c.coeff <- c.pmi
  if (coeff=="hyperg") c.coeff <- c.hyperg
  if (coeff=="ms") c.coeff <- c.ms

  # préparation du graphique
  mima <- range(c(R1,C1,O11))
  opar <- par(mar=par("mar"))
  on.exit(par(opar))
  par(mar=c(4.2,4.2,2,0.5))
  plot(mima, main = "", xlab = "", ylab = "Indices (base 100)",xlim=c(1,N),ylim=mima, type="n")
  axis(side=1, at=br2, labels=listdateb, line=2)
  title(main = mn, line=1, cex.main=.8, font.main=1, adj = 0)
  # mtext(sb, 3, line=0,cex=.8, font=1, adj=0)

  # affichage des valeurs observées
  lines(br2,R1,type="o", col="blue", lty=3, pch=21)
  lines(br2,C1, type="o", col="red", lty=3, pch=21)
  lines(br2,O11, type="o", col= "green", lty=2, pch=19)

  # affichage des tendances et de la légende
  mov <- lowess(c.coeff)
  lines(br2, mov$y, col="seagreen", lwd=4)
  mov1 <- lowess(R1)
  lines(br2, mov1$y, col="blue", lwd=1)
  mov2 <- lowess(C1)
  lines(br2, mov2$y, col="red", lwd=1)
  legend("bottomright", legend=c(lm1,lm2,coeff), cex=.7, lty=c(1,1,1),lwd=c(1,1,4), col=c("blue","red","seagreen"))

  return(HDD)
  #}
}

#####################################################################################

stablem <- function(corp="PL", lm, trnch="", mod="R") {

  # étude de la stabilité d'un lemme dans un corpus
  # moyenne et variances cumulées
  # NB les courbes aléatoires sont différentes à chaque essai :
  # les résultats de la fonction sample() sont différents à chaque éxécution
  # >> bootstrap rapide
  # AG 9/2014   GPL3

  library(rcqp, quietly=TRUE, warn.conflicts=FALSE)
  options(warn=-1)
  t1 <- Sys.time()

  LL <- "Lsub"

  # vérification
  provl <- paste(corp, ".lemma", sep="")
  idprov <- cqi_regex2id(provl, lm)
  if (length(idprov)==0) stop ("lemme inexistant dans ce corpus !", call.=FALSE)

  # récupération des occurrences
  quer1 <- paste ('[lemma="', lm, '"]', sep ="")
  cqi_query(corp, LL, quer1)
  sscorp <- paste(corp, LL, sep = ":")
  sscorp.dump <- cqi_dump_subcorpus(sscorp)
  sscorp.dump <- as.vector(sscorp.dump[,1])
  sscorp.ef <- length(sscorp.dump)
  corp.size <- size(corpus(corp))
  sscorp.var <- NULL
  sscorp.moy <- NULL
  if (trnch=="") {
    trnch <- round(corp.size/12000) # par défaut, tranches de 12000 tokens
  }

  # découpage
  br <- round(seq(1, corp.size, length.out=trnch+1))
  H <- hist(sscorp.dump, breaks=br, plot=FALSE)
  val.trnch <- H$counts
  rm(H)
  gc()

  # calcul de la variance cumulée
  if (mod=="A") { # mise en ordre aléatoire (avec bootstrap rapide)
    cat("Les tirages aléatoires successifs demandent du temps ! \n")
    sscorp.vartab <- matrix(0, nrow=trnch, ncol=100)
    for (j in 1:100) {
      val.trnch<- sample(val.trnch)
      for (i in 1:trnch) {
        sscorp.vartab[i,j] <- var(val.trnch[1:i])
      }
    }
    sscorp.var <- apply(sscorp.vartab, 1, mean)

    sscorp.moytab <- matrix(0, nrow=trnch, ncol=100)
    for (j in 1:100) {
      val.trnch<- sample(val.trnch)
      for (i in 1:trnch) {
        sscorp.moytab[i,j] <- mean(val.trnch[1:i])
      }
    }
    sscorp.moy <- apply(sscorp.moytab, 1, mean)
  }

  if (mod=="R") {
    for (i in 1:trnch) {
      sscorp.var[i] <- var(val.trnch[1:i])
    }

    for (i in 1:trnch) {
      sscorp.moy[i] <- mean(val.trnch[1:i])
    }
  }

  sstitre <- "ordre réel (dans le corpus)"
  if (mod=="A") {
    sstitre <- "ordre aléatoire (bootstrap rapide)"
  }

  # affichages
  layout(matrix(1:4, 2,2, byrow=TRUE))
  titre2 <- paste("corpus : ",corp," ; lemme : ", lm, " ; effectif total = ", sscorp.ef, sep="")
  H <- hist(sscorp.dump, breaks = br, freq=TRUE, border = "grey", main = "", xlab = "ordre du corpus", ylab = "Effectif par tranche")
  title(main = titre2, line=1.5, cex.main=1, font.main=1, adj = 0)
  mov <- lowess(H$counts)
  L <<- mov$y
  decal <- (corp.size/trnch)/2
  br2 <- (br[1:(length(br)-1)]) + decal
  B <<- br2
  lines(br2, mov$y, col="green", type="l", lwd=1.5)

  #dev.new()
  titranal <- "évolution de la moyenne cumulée"
  plot(x=(1:trnch),y=sscorp.moy, type="b", xlab = "nombre de tranches cumulées considérées", ylab="moyenne", col="blue", sub=sstitre)
  title(main = titre2, line=1.5, cex.main=1, font.main=1, adj = 0)
  mtext(titranal, 3, line=0.3,cex=1, font=1, adj=0)
  #dev.new()
  titranal <- "effectifs par tranche"
  hist(val.trnch, xlab="histogramme des fréquences par tranche", ylab="effectifs", main="", sub=sstitre)
  title(main = titre2, line=1.5, cex.main=1, font.main=1, adj = 0)
  mtext(titranal, 3, line=0.3,cex=1, font=1, adj=0)
  #
  titranal <- "évolution de la variance cumulée"
  plot(x=(1:trnch),y=sscorp.var, type="b", xlab = "nombre de tranches cumulées considérées", ylab="variance", col="red", sub=sstitre)
  title(main = titre2, line=1.5, cex.main=1, font.main=1, adj = 0)
  mtext(titranal, 3, line=0.3,cex=1, font=1, adj=0)

  stabvar <- as.list(NULL)
  stabvar$val <- val.trnch
  stabvar$cum <- sscorp.var

  t2 <- Sys.time()
  td <- difftime(t1,t2, units="secs")
  cat("\n","temps écoulé :",round(td), "secondes \n")

  return (stabvar)
}

