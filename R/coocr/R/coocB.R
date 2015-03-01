coocB <- function(tc, cex=.7, marg=4, prox="dice", selec=1){

# septembre 2012 AG -> ||domaine public||
# cooccurrents de deuxième ordre
# création de 2 lexicogrammes :
# > cooc dans le contexte uniquement
# > cooc dans le corpus entier

if (!inherits(tc, "COOCA"))
   stop("en entrée : un objet produit par coocA()", call.=F)
library(ade4, quietly=TRUE, warn.conflicts=FALSE)
library(circular, quietly=TRUE, warn.conflicts=FALSE)
library(MASS, quietly=TRUE, warn.conflicts=FALSE)
library(rcqp, quietly=TRUE, warn.conflicts=FALSE)
library(plyr, quietly=TRUE, warn.conflicts=FALSE)
options(warn = -1)
t1 <- Sys.time()
graphics.off()

# mises en place générales
#attach(t, warn.conflicts=FALSE)
corp <- tc$corpus
GDbon <- tc$cooc0
lnbs <- (tc$nb^2 - tc$nb)/2
ef <- tc$effectif
lico <- tc$cooc1[,4]
fr <- tc$cooc1[,5]
w1 <- as.numeric(tc$cooc1[,7])
attr <- tc$attr
val <- tc$val
poscooc <- tc$poscooc
pospiv <- tc$pospiv
nb <- tc$nb
dist <- tc$dist
effectif <- tc$effectif
cposLm <- tc$cposLm
lgcorp <- tc$lgcorp
tlex <- tc$tlex
txtitr <- paste("corpus : ", tc$corpus, "   bornes : ", tc$debfin, "   stock cible :  ", tc$NN,  "   lemme : ", tc$lemme, "   fréqu. : ", tc$effectif, "   fenêtre = +/-", tc$dist, "   nb cooccurrents =  ", tc$nb, sep="")
clr <- c("green3", "blue2", "red2", "magenta", "red", "darkcyan", "darkorange", "yellow4", "darkorchid", "darksalmon", "steelblue4", "sienna2", "ivory4")
df <- list(corpus=tc$corpus, debfin=tc$debfin, lgcorp=tc$lgcorp, lemme=tc$lemme, effectif=tc$effectif, nb=tc$nb, dis=tc$dist, cposLm=tc$cposLm, tlex=tc$tlex, cooc1=tc$cooc1, NN=tc$NN, A2d3 =tc$A2d3)

####################################################
# COOCCURRENCES DE DEUXIÈME ORDRE **T**
# = dénombrement des paires cooccurrentes
# dans le corpus entier

# -> effectifs dans un tableau simple
# calcul de la distance de Dice pour toutes les paires

c <- corpus(corp)
Clm <- paste(corp, ".lemma", sep="")

larg <- length(lico)
vdis1 <- vdis2 <- vdis3 <- vdis4 <- vdis5 <- vdis6 <- 0
n <- 0

###### !

# cas simple
if (attr=="" & val=="") {
for (i in 1:(larg-1)) {
  prem <- i + 1
  cat("*")
  for (j in prem:larg) {
    n <- n + 1
    a <- lico[i]
    b <- lico[j]
    dem <- paste("MU(meet[lemma = '" , a , "'][lemma = '" , b , "']"," -", marg+1, " " , marg+1,")", sep="")
    N <- subcorpus(c,dem)
    eff <- size(N)
    vdis3[n] <- eff
    vdis1[n] <- lico[i]
    vdis2[n] <- lico[j]
    vdis4[n] <- fr[i]
    vdis5[n] <- fr[j]
    }
  }
}

# ici une possible vectorisation
# effet pratique quasi nul...
#tbin <- cbind(vdis1, vdis2)
#calcooc <- function(x, marg=marg) {
#  a <- x[1]
#  b <- x[2]
#  dem1 <- paste("[lemma = '" , a , "'][]{0," , marg , "}[lemma = '" , b , "']",  sep="")
#  dem2 <- paste("[lemma = '" , b , "'][]{0," , marg , "}[lemma = '" , a , "']",  sep="")
#  N1 <- subcorpus(c, dem1)
#  N2 <- subcorpus(c, dem2)
#  eff <- size(N1) + size(N2)
#  return(eff)
#}
#vdis3 <- apply(tbin, 1, calcooc, marg)

# cas complexe
if (attr!="" & val!="") {
for (i in 1:(larg-1)) {
  prem <- i + 1
  cat("*")
  for (j in prem:larg) {
    n <- n + 1
#    cpt <- cpt + 1
    a <- lico[i]
    b <- lico[j]
    dem1 <- paste("[lemma = '" , a , "'][]{0," , marg , "}[lemma = '" , b , "']", "::match.", attr, "=\"", val, "\"", sep="")
    dem2 <- paste("[lemma = '" , b , "'][]{0," , marg , "}[lemma = '" , a , "']", "::match.", attr, "=\"", val, "\"", sep="")
    N1 <- subcorpus(c, dem1)
    N2 <- subcorpus(c, dem2)
    eff <- size(N1) + size(N2)
    vdis3[n] <- eff
    vdis1[n] <- lico[i]
    vdis2[n] <- lico[j]
    vdis4[n] <- fr[i]
    vdis5[n] <- fr[j]
    }
  }
}

vdis6 <- vdis3 * 1000 / (vdis4 + vdis5)^selec
vdis6 <- sprintf("%06.3f", vdis6)
vdis <- cbind(vdis1, vdis2, vdis3, vdis4, vdis5, vdis6)
vdis <- as.data.frame(vdis, stringsAsFactors=FALSE)
names(vdis) <- c("lm1", "lm2", "cooc", "ef_lm1", "ef_lm2", "Dice")

# création de la matrice des distances
cat("\n")
dice <- vdis[,6]
bid <- rep(0,larg^2)
tabdisLm <- matrix(bid, nrow=larg, ncol=larg)
n <- 0
for (i in 1:(larg-1)) {
  prem <- i + 1
  for (j in prem:larg) {
    n <- n + 1
    tabdisLm[i,j] <- dice[n]
    tabdisLm[j,i] <- dice[n]
    }
  }
tabdisLm <- as.numeric(tabdisLm)
tabdisLm <- matrix(tabdisLm, nrow=larg, ncol=larg)
rownames(tabdisLm) <- lico
colnames(tabdisLm) <- lico

vdis <- vdis[rev(order(vdis[,6])),]
df$coocT <- vdis
df$distmatT <- tabdisLm
tabdisLm <- log1p(tabdisLm ^ .2)

# création d'un tableau de contingence des contextes
# une ligne par contexte, nbs colonnes (=coocccurrents retenus)
# le 1 signifie que le lemme est présent dans le contexte considéré
idlem <- tc$cooc1[,c(1,4)]
idlem <- idlem[(order(idlem[,2])),]
idutils <- idlem[,1]
tlexb <- tlex %in% idutils
tlexb <- matrix(tlexb, ncol=ncol(tlex), nrow=nrow(tlex))
tlexS <- apply(tlexb, 1, sum)
tlexu <- tlex[(tlexS>1),]
tlexcl <- alply(tlexu, 1, match, idutils)
tlexcl <- t(as.data.frame(tlexcl))

cretab <- function(x) {
  lig0 <- rep(0,50)
  lig1 <- x[!is.na(x)]
  lig0[lig1] <- 1
  return (lig0)
}

tabsup <- alply(tlexcl, 1, cretab)
tabsup <- t(as.data.frame(tabsup))
tabsup <- tabsup[,(1:nb)]
colnames(tabsup) <- idlem[,2]

# analyse en composantes principales
#  fréquences totales (dans le corpus)
tab.pca <- dudi.pca(tabdisLm, scannf=FALSE, nf=4)
contrib <- inertia.dudi(tab.pca, row.inertia=TRUE, col.inertia=TRUE)
df$ACPt <- tab.pca
df$CONTRt <- contrib
contrabs <- apply(contrib$col.cum, 1, sum)
df$CUMTOTt <- sort(contrabs)
#(hisden(contrabs))
T0 <- tab.pca$co
ymax <- max(T0[,2])
ymin <- min(T0[,2])
Tbon <- lisible(T0[,1],T0[,2],lab=row.names(T0),mn=ymin, mx=ymax,cex=(cex+.4))
T0[,1] <- Tbon[,1]
T0[,2] <- Tbon[,2]
T0[,5] <- tc$cooc1[,2]
poscoo <- unlist(strsplit(poscooc, " "))
nbpos <- length(poscoo)
for (i in 1:nbpos) {
  T0[(T0[,5]==poscoo[i]),6] <- clr[i]
}

# calcul des segments
if (prox=="spearman") {
  cortabdist <- cor(tabdisLm, method="spearman")}
if (prox=="dice") {
  cortabdist <- tabdisLm }
cortabdist <- 1 - cortabdist
rempl <- as.vector(rep(1, lnbs*3))
discor <- matrix(rempl, nrow=lnbs, ncol=3)
n <- 0
for (i in 1:(nb-1)) {
  prem <- i + 1
  for (j in prem:nb) {
    n <- n + 1
    discor[n,3] <- cortabdist[i,j]
    discor[n,1] <- rownames(cortabdist)[i]
    discor[n,2] <- colnames(cortabdist)[j]
    }
  }

discor <- discor[(order(discor[,3])),]
discor <- as.data.frame(discor)
R1r <- match(discor[,1], rownames(T0))
R2r <- match(discor[,2], rownames(T0))
discor[,4] <- T0[R1r, 1]
discor[,5] <- T0[R1r, 2]
discor[,6] <- T0[R2r, 1]
discor[,7] <- T0[R2r, 2]
discor1 <- discor[(1:12),]
discor2 <- discor[13:nb,]

# affichage des points-colonnes (lemmes)
# graphics.off()
dev.new()
titranal <- "ACP sur tableau des distances de Dice entre co-cooccurrents (corpus cible entier) + distances de Dice"
opar <- par(mar=par("mar"))
on.exit(par(opar))
par(mar=c(0.5,0.5,1.7,0.5))
(plot(T0[,1], T0[,2], type="n", xlab="", ylab="", axes=FALSE, frame.plot=TRUE))
axis(side=1, labels=FALSE, tick=0)
axis(side=2, labels=FALSE, tick=0)
title(main = txtitr, line=1, cex.main=.8, font.main=1, adj = 0)
mtext(titranal, 3, line=0,cex=.8, font=1, adj=0)
segments(discor1[,4], discor1[,5], discor1[,6], discor1[,7], lwd=1.5, col="grey")
segments(discor2[,4], discor2[,5], discor2[,6], discor2[,7], lwd=.4, col="grey")

text(T0[,1], T0[,2], labels=row.names(T0), col=T0[,6], cex=cex)

# calcul des points-lignes
ttabsup <- t(tabsup)
suplignes <- supcol(tab.pca, ttabsup)


########################################################
# COOCCURRENCES DE DEUXIÈME ORDRE **C**
# = dénombrement des paires cooccurrentes
# calculs dans le contexte uniquement

if (ef < 200) {
  df$lgcorp <- lgcorp
  class(df) <- "COOCB"
  return(df)
}

# valeurs fixes
larg <- dist * 2
npair <- larg * (larg-1) / 2
long1 <- effectif
long2 <- effectif * npair
n <- 0
BINOM1 <- NULL
BINOM2 <- NULL
TBINOM <- matrix(seq(0, length(cposLm)*larg), nrow=length(cposLm), ncol=larg)

# création du tableau des ids des cooccurrents dans la fenêtre
for (i in dist:1) {
  cposG <- cposLm - i
  Gid <- cqi_cpos2id(Clm, cposG)
  n <- n+1
  TBINOM[,n] <- Gid
  }
  for (i in 1:dist) {
  cposD <- cposLm + i
  Did <- cqi_cpos2id(Clm, cposD)
  n <- n+1
  TBINOM[,n] <- Did
  }

# création de la suite des paires de cooccurrents
for (i in 1:(larg-1)) {
  cat (".", sep = "")
  prem <- i + 1
  for(j in prem:larg) {
    BINOM1 <- c(BINOM1,TBINOM[,i])
    BINOM2 <- c(BINOM2,TBINOM[,j])
  }
}
BINOM<- cbind(BINOM1, BINOM2)

# premier nettoyage
BINOM <- BINOM[(BINOM[,1]!=BINOM[,2]),]
BINOM <- BINOM[(BINOM[,1]%in%idutils),]
BINOM <- BINOM[(BINOM[,2]%in%idutils),]

if (nrow(BINOM) == 0) {
  cat("Pas de co-cooccurrents pour cette requête...")
  return (df)
  }

# permutations B-A -> A-B (si B < A)
# de manière à compter ensemble A-B et B-A
BINOMa <- BINOM[(BINOM[,1]<BINOM[,2]),]
BINOMb <- BINOM[(BINOM[,1]>BINOM[,2]),]
BINOMa <- as.data.frame(BINOMa)
BINOMb <- as.data.frame(BINOMb)
BINOMb[,3] <- BINOMb[,2]
BINOMb[,4] <- BINOMb[,1]
BINOMb <- BINOMb[,3:4]
names(BINOMb) <- c("BINOM1", "BINOM2")
BINOM  <- rbind(BINOMa, BINOMb)

# comptage des fréquences des paires
BINOM <- as.data.frame(BINOM)
liste <- as.data.frame(aggregate(1:nrow(BINOM), BINOM, length))
lnmax <- (nb^2 - nb)/2
df$densite <- nrow(liste)/lnmax

# compléments au tableau général et mise en forme
# liste <- liste[liste[,3]>4,]
liste[,4] <- cqi_id2str(Clm, liste[,1])
liste[,5] <- cqi_id2str(Clm, liste[,2])
i1 <- match(liste[,1], GDbon[,1])
liste[,6] <- GDbon[i1,5]
i2 <- match(liste[,2], GDbon[,1])
liste[,7] <- GDbon[i2,5]
Dice <- (liste[,3] * 1e+10) / (liste[,6] + liste[,7])^2
N <- sqrt(sum(Dice^2))
Dice <- Dice / N
liste[,8] <- Dice
liste <- liste[rev(order(liste[,8])),]

liste <- liste[(order(liste[,4])),]
nam <- c("id1", "id2", "F.cooc", "lem1", "lem2", "F.lm1", "F.lm2", "Dice")
names(liste) <- nam
or <- seq(1, nrow(liste))
row.names(liste) <- or

df$cooc2 <- liste
cooc2b <- liste[,c(3,4,5)]
cooc2b <- cooc2b[rev(order(cooc2b[,1])),]
cat("\n")
write.matrix(cooc2b[1:30,])

# CRÉATION DE LA MATRICE DE DISTANCES ENTRE LES COOCCURRENTS =====
#  = distances entre les "co-cooccurrents
li1 <- liste[,c(1,4)]
li2 <- liste[,c(2,5)]
nli <- c("ids", "lem")
names(li1) <- nli
names(li2) <- nli
listep1 <- rbind(li1, li2)
listep1 <- unique(listep1)
listep1 <- listep1[order(listep1[,2]),]
row.names(listep1) <- seq(1, nrow(listep1))
xn <- rep(0, (nrow(listep1)^2))
tabdisLm <- matrix(xn, nrow=nrow(listep1), ncol=nrow(listep1))
colnames(tabdisLm) <- listep1[,2]
rownames(tabdisLm) <- listep1[,2]

for (i in 1:nrow(liste)) {
  tabdisLm[(liste[i,4]),(liste[i,5])] <- liste[i,8]
  tabdisLm[(liste[i,5]),(liste[i,4])] <- liste[i,8]
  }
df$distmatc <- tabdisLm
tabdisLm <- log1p(tabdisLm ^ .2)

w <- match(colnames(tabdisLm), tc$cooc1[,4])
w1 <- tc$cooc1[w,7]
w1 <- as.numeric(w1)
df$rowweights <- w1

# analyse en composantes principales
# fréquences dans le contexte uniquement
tab.pca <- dudi.pca(tabdisLm, scannf=FALSE, nf=4)
contrib <- inertia.dudi(tab.pca, row.inertia=TRUE, col.inertia=TRUE)
df$ACPc <- tab.pca
df$CONTRc <- contrib
contrabs <- apply(contrib$col.cum, 1, sum)
df$CUMTOTc <- sort(contrabs)
T0 <- tab.pca$co
ymax <- max(T0[,2])
ymin <- min(T0[,2])
Tbon <- lisible(T0[,1],T0[,2],lab=row.names(T0),mn=ymin, mx=ymax,cex=(cex+.4))
T0[,1] <- Tbon[,1]
T0[,2] <- Tbon[,2]
Tr <- match(row.names(tabdisLm), tc$cooc1[,4])
T0[,5] <- tc$cooc1[Tr,2]
poscoo <- unlist(strsplit(poscooc, " "))
nbpos <- length(poscoo)
for (i in 1:nbpos) {
  T0[(T0[,5]==poscoo[i]),6] <- clr[i]
}

# calcul des segments
if (prox=="spearman") {
  cortabdist <- cor(tabdisLm, method="spearman") }
if (prox=="dice") {
  cortabdist <- tabdisLm }

cortabdist <- 1 - cortabdist
rempl <- as.vector(rep(1, lnbs*3))
discor <- matrix(rempl, nrow=lnbs, ncol=3)
n <- 0
nb <- nrow(cortabdist)
for (i in 1:(nb-1)) {
  prem <- i + 1
  for (j in prem:nb) {
    n <- n + 1
    discor[n,3] <- cortabdist[i,j]
    discor[n,1] <- rownames(cortabdist)[i]
    discor[n,2] <- colnames(cortabdist)[j]
    }
  }

discor <- discor[(order(discor[,3])),]
discor <- as.data.frame(discor)
R1r <- match(discor[,1], rownames(T0))
R2r <- match(discor[,2], rownames(T0))
discor[,4] <- T0[R1r, 1]
discor[,5] <- T0[R1r, 2]
discor[,6] <- T0[R2r, 1]
discor[,7] <- T0[R2r, 2]
discor1 <- discor[(1:12),]
discor2 <- discor[13:nb,]

# affichage
dev.new()
titranal <- "ACP sur tableau des distances de Dice entre co-cooccurrents (contexte) + distances de Dice"
opar <- par(mar=par("mar"))
on.exit(par(opar))
par(mar=c(0.5,0.5,1.7,0.5))
(plot(T0[,1], T0[,2], type="n", xlab="", ylab="", axes=FALSE, frame.plot=TRUE))
axis(side=1, labels=FALSE, tick=0)
axis(side=2, labels=FALSE, tick=0)
title(main = txtitr, line=1, cex.main=.8, font.main=1, adj = 0)
mtext(titranal, 3, line=0,cex=.8, font=1, adj=0)
segments(discor1[,4], discor1[,5], discor1[,6], discor1[,7], lwd=1.5, col="grey")
segments(discor2[,4], discor2[,5], discor2[,6], discor2[,7], lwd=.4, col="grey")

text(T0[,1], T0[,2], labels=row.names(T0), col=T0[,6], cex=cex)

# calcul des points-lignes
tabsup <- tabsup[,(listep1[,2])]
ttabsup <- t(tabsup)
df$tabsup2 <- tabsup
suplignes <- supcol(tab.pca, ttabsup)
df$tabsuplignes2 <- suplignes


t2 <- Sys.time()
td <- difftime(t1,t2, units="secs")
cat("\n","temps écoulé :",round(td), "secondes \n")

return(df)
}
