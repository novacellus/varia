
coocE <- function(tc, trnch=3, cex=.8) {

  # septembre 2012 AG -> ||GPL3||
  # table des lemmes cooccurrents divisée en tranches chronologiques
  # de même taille, classées par centres de gravité,
  # graphe de l'évolution des fréquences
  # + graphe Bertin-Cibois

  library(ade4, quietly=TRUE, warn.conflicts=FALSE)
  if (!inherits(tc, "COOCA"))
    stop("en entrée : un objet produit par coocA()", call.=F)
  attach(tc, warn.conflicts=FALSE)
  # tabdisLm <- distmatc
  idutils <- cooc1[,1]
  ef <- length(cposLm)


  #### moyennes réciproques (sur colonnes)

  mrecip1 <- function(x) {
    li.m <- rep(0, nrow(x))
    for (j in 1:ncol(x)){
      for (i in 1:nrow(x)){
        li.m[i] <- li.m[i]+(x[i,j]*j)
      }
    }
    li.m <- li.m / rowSums(x)
    x <- x[sort.list(li.m),]
    return(x)
  }

  ### affichage méthode Bertin-Cibois
  ### (assocplot() modifié)

  repfac2 <- function (x, col = c("red", "white"), border= c("red", "blue"), main = NULL,
                       xlab = NULL, ylab = NULL)  {
    opar <- par
    par(mar = c(2.5,4.8,2.0,0.1))
    n <- sum(x)
    f <- x[, NCOL(x):1]
    e <- outer(rowSums(f), colSums(f), "*")/n
    d <- as.matrix(f - e)/sqrt(e)
    e <- sqrt(e)
    x.w <- apply(e, 1, max)
    y.h <- apply(d, 2, max) - apply(d, 2, min)
    ylim <- c(0, sum(x.w))
    xlim <- c(0, sum(y.h))
    plot.new()
    plot.window(xlim, ylim, log = "")
    x.r <- cumsum(x.w)
    x.m <- (c(0, x.r[-NROW(f)]) + x.r)/2
    y.u <- cumsum(y.h)
    y.m <- y.u - apply(pmax(d, 0), 2, max)
    x.m <- rev(x.m)
    y.m <- rev(y.m)
    z <- expand.grid(x.m, y.m)
    rect(z[, 2], z[, 1] - e/2, z[, 2] + d, z[, 1] + e/2,col = col[1 + (d < 0)], border= border[1+(d<0)], lwd=2)
    txtitr <- paste("corpus : ", corpus, "   bornes : ", debfin, "   lemme : ", lemme, "   fréqu. : ", effectif, "   fenêtre = +/-", dist, "   nb cooccurrents =  ", nb, sep="")
    titranal <- "Graphe de Bertin-Cibois de la répartion des cooccurrents par tranches chronologiques (écarts à l'indépendance)"
    axis(2, at = x.m, labels = rownames(f), tick = FALSE, cex.axis=0.8, las=2)
    axis(1, at = y.m, labels = colnames(f), tick = FALSE, cex.axis=1.0, las=1)
    abline(v = y.m, lty = 2)
    ndn <- names(dimnames(f))
    if (length(ndn) == 2) {
      if (is.null(xlab))
        xlab <- ndn[1]
      if (is.null(ylab))
        ylab <- ndn[2]
    }
    title(main = txtitr, line=1, xlab = xlab, ylab = ylab, cex.main=.8, font.main=1, adj=0)
    mtext(text=titranal, 3, line=0, cex=.8, font=1, adj=0)
    par <- opar
  }
  ########################################################################

  # création d'un tableau lexical sélectionné et trié ======

  tablexT <- matrix(1, nrow=length(cposLm), ncol=nb)
  for (i in 1:nb) {
    tabpro <- (tlex %in% idutils[i])
    tabpro <- matrix(tabpro, nrow=length(cposLm), ncol=dist*2)
    stab <- apply(tabpro, 1, sum)
    tablexT[,i] <- stab
  }
  colnames(tablexT) <- cooc1[,4]
  rownames(tablexT) <- cposLm
  T <- tablexT

  if (ef < 100) {
    cat("\n" , "ATTENTION ! moins de 100 occurrences du lemme choisi", "\n")
    cat("   Résultats à interpréter avec une extrême prudence !")
  }

  # découpage et traitement du tableau
  # (on découpe le tableau en n groupes comportant le même nombre de lignes)
  T <- tablexT
  n <- trnch
  eff <- nrow(T)
  rt <- seq(1, eff, 1)
  brn <- quantile(rt, probs=seq(0,1,1/n), type=1)
  L <- list(NULL)
  for (i in 1:n) {
    L[[i]] <- T[(brn[i]:brn[i+1]),]
  }

  rp <- (rep(0, (n*ncol(tablexT)) ))
  TabC <- matrix(rp, ncol=n, nrow=ncol(tablexT))
  for (i in 1:n) {
    S <- apply(L[[i]], 2, sum)
    TabC[,i] <- S
  }
  rownames(TabC) <- cooc1[,4]
  cnam <- NULL
  for (i in 1:trnch) {
    rn <- paste("P", i, sep="")
    cnam <- c(cnam, rn)
  }
  colnames(TabC) <- cnam
  txtitr <- paste("corpus : ", corpus, "   bornes : ", debfin, "   lemme : ", lemme, "   fréqu. : ", effectif, sep="")

  # premier graphique : chronologie générale
  if (attr=="" & val=="") {
    opar <- par(mar=par("mar"))
    on.exit(par(opar))
    par(mar=c(4.0,3.4,1.7,0.5))
    trnch2 <- 20
    cposbrn <- cposLm[brn]
    br <- seq(1, NN, length.out=trnch2+1)
    decal <- (NN/trnch2)/2
    #mn <- paste("lemme : ", lemmme, " (effectif du lemme = ", effectif, ")", sep = "")
    sb <- paste("index du corpus cible (effectif total = ", NN, ")", sep="")
    H <- hist(cposLm, breaks = br, border = "grey", main = "", xlab = sb, ylab = "Effectif par tranche", cex.lab=.8, cex.axis=.7, mgp=c(2,1,0))
  }

  if (attr!="" & val!="") {
    opar <- par(mar=par("mar"))
    on.exit(par(opar))
    par(mar=c(4.0,3.4,1.7,0.5))
    trnch2 <- 20
    cposbrn <- cposLm[brn]
    ef <- length(cposLm)
    indexL <- as.integer(seq(1, ef, 1))
    efI <- cbind(as.integer(cposLm, indexL))
    # ici, une ruse diabolique !
    ln2 <- nrow(A2d3)
    indexT <- as.integer(seq(1, ln2, 1))
    corpI <- cbind(as.integer(A2d3[,1]), indexT)
    choix <- corpI[,1] %in% efI[,1]
    corpI2 <- cbind(corpI, choix)
    indbon <- corpI2[(corpI2[,3]==TRUE),2]
    cposbrn <- indbon[brn]
    #mn <- paste("lemme : ", lm, " (effectif du lemme = ", ef, ")", sep = "")
    sb <- paste("index du corpus cible (effectif total = ", NN, ")", sep="")
    br <- seq(1, ln2, length.out=trnch2+1)
    decal <- (ln2/trnch2)/2
    H <- hist(indbon, breaks = br, border = "grey", main = "", xlab = sb, ylab = "Effectif par tranche", cex.lab=.8, cex.axis=.7, mgp=c(2,1,0))
  }

  mov <- lowess(H$counts)
  br2 <- (br[1:(length(br))-1]) + decal
  lines(br2, mov$y, col="blue", lwd=1.5)
  for (i in 2:(trnch)){
    points(cposbrn[i] ,0, pch="|",cex=3, col="turquoise")
  }
  title(main = txtitr, line=0, cex.main=.8, font.main=1, adj = 0)

  # mise en ordre du tableau et calcul des écarts à l'indépendance
  # affichage des écarts ordonnés
  TAB <- mrecip1(TabC)
  dev.new()
  repfac2(TAB)
  chi2 <- chisq.test(TAB)
  dif <- round(chi2$observed-chi2$expected)
  scT <- apply(TAB, 2, sum)
  TAB <- rbind(TAB, scT)
  srT <- apply(TAB, 1, sum)
  TAB <- cbind(TAB, srT)

  # affichages et sauvegardes
  sink(file=lemme, append=TRUE)
  cat("\n", date(), "\n\n")
  print(TAB)
  cat("\n\n")
  print(dif)
  sink(file=NULL)
  cat("\n\n")
  result <- cbind(row.names(TAB), TAB)
  write.matrix(result, file="")
  res <- list(observe=TAB, difference=dif)
  class(res) <- "COOCC"

  return(res)
}
