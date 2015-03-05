#' Wrapper function for computing coocurrences (mainly for interactive use)
#' The function parses user input, checks for existing corpora, formulates query phrase
#' @param corp Corpus name
#' @param piv Searched word
#' @param poscooc POS of the coocurences
#' @param dis Coocurence window
#' @param coeff Measure
#' @param selec ???
#' @param nbs Number of coocurrences
#' @param D Start (index of corpus position)
#' @param F (index of corpus position)
#' @param attr Structural attribute to filter the results
#' @param val Value of the structural attribute
#' @param corp_lemma Name of the positional attribute for lemma
#' @param corp_pos Name of the positional attribute for POS
#' @return 'cooca' object
#' @examples
#' cooc <- coocA(corp="PATROLOGIA",piv="hortus",dis=3)
#' @export
coocA <- function(corp, piv, poscooc="QLF SUB VBE",
                  dis=5, coeff="dice", selec=1, nbs=30,
                  D=0, F=Inf, attr="", val="",
                  corp_lemma="lemma", corp_pos="pos"
                  ) {
  t1 <- Sys.time()
# vérifications de base
  # du corpus
  if(! .check_corpus_name(corp))
     stop(.error_message("no_such_corpus",corpus=corp))
  # des attributes
  p_attr_list <- .list_attrs(corp)$p
  if(!corp_lemma %in% p_attr_list)
    stop(.error_message("no_p_attr",corp_lemma=corp_lemma,p_attr_list=p_attr_list))
  corp.lm <- paste(corp, ".", corp_lemma , sep="")

  if(!corp_pos %in% p_attr_list)
    stop(.error_message("no_p_attr",corp_pos=corp_pos,p_attr_list=p_attr_list))
  corp.ps <- paste(corp, ".", corp_pos , sep="")

  #du pivot --> id
  piv.id <- find_word(corp=corp,corp_lemma=corp_lemma,piv=piv)

  # du coefficient
  lst.coeff <- c("dice","poiss","daille","pmi","hyperg","mins")
  if (! coeff %in% lst.coeff) stop ("coefficient inconnu, choisir : daille/dice/hyperg/mins/pmi/poiss", call.=FALSE)

  # de la borne supérieure
  CORP <- rcqp::corpus(corp) # corpus complet
  CORP.ln <- rcqp::size(CORP)
  lgcorp <- CORP.ln
  if (F == Inf) F <- CORP.ln
  if (F > CORP.ln) F <- CORP.ln
  lgcorp <- F-D

  # de la paire attribut/valeur (validité de la paire : AGENDUM !!!)

  #
  # établissement de la liste des cpos du pivot
  #
  LL <- "Lsub" # nom générique de sous-corpus
  A2 <- ""
  dump.pivot <- find_words(piv.id=piv.id,corp.lm=corp.lm,piv=piv,corp=corp)
  #if ((attr!="" & val=="") | (attr=="" & val!="")) {
  #  stop ("Requête incomplète !", call.=FALSE)
  #  cat ("\n Vérifications simples : OK ! \n")
  #} else if (attr!="" & val!="") {
  # cat("La recherche du pivot dans une fraction du corpus total prend du temps (beaucoup)... \n")
  # LL2 <- "Lsub2"
  # quer2 <- paste("[lemma='", piv, "' %cd]", "::match.", attr, "=\"", val, "\"", sep="")
  # rcqp::cqi_query(corp, LL2, quer2)
  # sscorp2 <- paste(corp, LL2, sep = ":")
  # if (rcqp::cqi_subcorpus_size(sscorp2)==0) stop ("Pas de réponse, vérifiez la requête !", call.=F)
  # dump.pivot <- rcqp::cqi_dump_subcorpus(sscorp2)[,1]
  # rcqp::cqi_drop_subcorpus(sscorp2)
  #} else if (attr=="" & val=="") {

   # dump.pivot <-rcqp::cqi_id2cpos(corp.lm,piv.id)
  #}
  gc()
  cpos.pivot <- dump.pivot[dump.pivot>D & dump.pivot <F] # tranche simple
  eff.pivot <- length(cpos.pivot)
  if (eff.pivot == 0) stop("Aucune occurrence du pivot dans la sélection !")
  cat (" cpos du pivot : OK ! \n")
  pos.piv <- rcqp::cqi_cpos2str(corp.ps ,cpos.pivot[1])

  #
  # établissement de la liste des POS retenus (lst.pos.r)
  # acrobatie nécessaire pour le cas des pos subdivisés ("radical"+subdivision)
  # e.g. VER = VER.cond + VER.futu + VER.impe + VER.impf .....
  # > il faut passer des "c.hoisis" aux "r.éels"
  #

  lst.pos.tt <- rcqp::cqp_flist(CORP, "pos") # liste complète des POS du corpus
  lst.pos.tt <- sort(names(lst.pos.tt))
  pos.coo <- unlist(strsplit(poscooc, " ")) # liste des POS choisis
  nbpos.c <- length(pos.coo) # nombre des POS choisis
  lst.pos.ri <- NULL
  lst.pos.r <- NULL
  lst.pos.c <- NULL

  for (i in 1:nbpos.c) {  # on traite "radical" par "radical"
    lpp <- grep(pos.coo[i], lst.pos.tt) # vecteur des POS correspondant à un "radical" (indices)
    lpp2 <- rep(pos.coo[i], length(lpp))# liste parallèle des POS choisis
    lst.pos.ri <- c(lst.pos.ri, lpp) # constitution par étape du vecteur (indices) des POS r.éels
    lst.pos.c <- c(lst.pos.c, lpp2) # idem pour les POS c.hoisis (vecteur strings)
  }

  lst.pos.r <- lst.pos.tt[lst.pos.ri] # passage à un vecteur de strings
  corr.pos.rc <- cbind(as.vector(lst.pos.r),as.vector(lst.pos.c)) # tableau de correspondance réels/choisis (strings)
  nb.pos.r <- 1
  if (length(lst.pos.r) > 1) {     # mises en ordre
    corr.pos.rc <- corr.pos.rc[order(corr.pos.rc[,1]),]
    lst.pos.r <- corr.pos.rc[,1]
    nb.pos.r <- length(lst.pos.r)
  }
  cat (" Liste des POS réels : OK ! \n")

  #
  # construction de la table des cooccurrents (ids)
  # 2 sorties : * tableau lexical (une ligne par occ du pivot) * vecteur cumulatif
  #
  dis.matrix <- c(-dis:-1,1:dis)
  ids.tablex <- vapply(dis.matrix, '+', cpos.pivot, FUN.VALUE = numeric(length(cpos.pivot)) )
  ids.GD <- rcqp::cqi_cpos2id(corp.lm,cpos = as.vector(ids.tablex) )
  ps.GD <- rcqp::cqi_cpos2str(corp.ps,as.vector(ids.tablex) )
  ids.tablex <- matrix(as.numeric(ids.GD), nrow=eff.pivot, ncol=(dis*2))
  cooc.GD <- cbind(ids.GD,ps.GD) # ids et pos correspondants, 2 colonnes ([1,] "1"    "PON"
  #cat("***-")

  #
  # tri des ids en fonction des POS sélectionnés
  #

  cooc.GD.pos <- matrix(FALSE, nrow=eff.pivot, ncol=(dis*2)) # on joue sur le fait qu'une matrice est un vecteur...
  for (i in 1:nb.pos.r) {
    logx <- ps.GD == lst.pos.r[i] # vecteur ajoute TRUE (=1) pour chaque POS sélectionné
    cooc.GD.pos <- cooc.GD.pos + logx #
  }
  cooc.GD.pos <- as.logical(cooc.GD.pos) # vecteur de sélection
  cooc.GD.tr <- cooc.GD[cooc.GD.pos,] # on ne garde que ce qui est TRUE
  for (i in 1:nrow(corr.pos.rc)) {
    asub <- corr.pos.rc[i,1]
    cooc.GD.tr[(cooc.GD.tr[,2]==asub),2] <- corr.pos.rc[i,2] # table à 2 colonnes : ids/POS (noms choisis) (triée)
    cat("***-")
  }

  #
  # calcul des fréquences par lemme et troncature
  #

  cooc.GD.tr <- as.data.frame(cooc.GD.tr)
  cooc.final <- plyr::count(cooc.GD.tr,c("ids.GD","ps.GD")) # on fait la somme par lemme > 3e colonne
  cooc.final <- cooc.final[rev(order(cooc.final[,3])),]
  cooc.final <- cooc.final[cooc.final[,3] > 2,] # troncature : au moins 3 cooccurrents
  cat("***-")


  # compléments et mise en forme du tableau
  id.n <- as.character(cooc.final[,1])
  id.n <- as.numeric(id.n)
  cooc.final[,1] <- id.n
  cooc.final[,4] <- rcqp::cqi_id2str(corp.lm, id.n) # on ajoute le lemme (string) > 4e colonne
  cooc.final <- cooc.final[(cooc.final[,4]!=piv),]
  cooc.final <- cooc.final[(cooc.final[,4]!="--"),]
  cooc.final <- cooc.final[(cooc.final[,4]!="-"),]
  cooc.final <- cooc.final[(cooc.final[,4]!="unknown"),]

  sup <- rep(1, nrow(cooc.final)) # nettoyage des lemmes parasites restants (?)
  cooc.final[,5] <- sup
  cooc.final <- cooc.final[order(cooc.final[,4]),]
  cat("***01-")

  if (nrow(cooc.final)<2) stop(paste(piv, "/", eff.pivot,"/", "pas de cooccurrent significatif !"), call.=FALSE)

  for (i in 1:(nrow(cooc.final)-1)) {
    if (cooc.final[i,4] == cooc.final[i+1,4]){
      cooc.final[i,5] <- 0
      cooc.final[i+1,5] <- 0
    }
  }
  cat("**02-")
  cooc.final <- cooc.final[(cooc.final[,5]==1),]

  #
  # calcul des fréquences *totales* des lemmes retenus dans le corpus considéré
  # difficulté : les sélections de sous-corpus (par cpos ou/et par attribut)
  #

  if (attr=="" & val=="" & D==0 & F==CORP.ln) { #    cas élémentaire : corpus complet
    cat("**1-", "\n")
    cooc.final[,5] <- rcqp::cqi_id2freq(corp.lm, cooc.final[,1])
    eff.tt.final <- CORP.ln
    cat(eff.tt.final, "\n")
  }

  if (attr=="" & val=="" & (D!=0 | F!=CORP.ln)) { #    simple troncature par cpos
    cat("**2-")
    options(scipen=999)  # supprimer la notation scientifique (pb avec cqp)
    quer5 <- paste("abc:[lemma=\".*\" %cd]::abc >=",D," & abc <=",F, sep="")
    LL3 <- "Lsub3"
    rcqp::cqi_query(corp, LL3, quer5) # création du sous-corpus ad hoc
    sscorp3 <- paste(corp, LL3, sep=":")
    eff.tt.final <- rcqp::cqi_subcorpus_size(sscorp3)
    cat(eff.tt.final, "\n")
    lst.freq <- rcqp::cqi_fdist1(sscorp3, "match", "lemma") # liste de fréquences du ss-corpus sélectionné
    cat("***-")
    lst.str <- rcqp::cqi_id2str(corp.lm, lst.freq[,1]) # on passe aux strings
    cat("***-")
    lst.freq <- as.data.frame(cbind(lst.freq, lst.str))
    lst.freq <- lst.freq[(lst.freq[,1] %in% cooc.final[,1]),] # on ne retient que les lemmes qui cooccurrent
    cat("***-", "\n")
    lst.freq <- lst.freq[(order(lst.freq[,3])),] # ordre alpha des lemmes
    cooc.final[,5] <- as.numeric(as.character(lst.freq[,2])) # on ajoute les fréquences globales au tableau
    cat("***-", "\n")
  }

  if (attr!="" & val!="") { # sélections
    cat("**3-", "\n")
    options(scipen=999)  # supprimer la notation scientifique (pb avec cqp)
    quer5 <- paste("abc:[lemma=\".*\" %cd & _.", attr, "=\"", val, "\"]::abc >=",D," & abc <=",F, sep="")
    LL3 <- "Lsub3"
    rcqp::cqi_query(corp, LL3, quer5) # création du sous-corpus ad hoc
    lgcorp <- rcqp::cqi_subcorpus_size(paste(corp,LL3, sep=":"))
    cat("***-")
    sscorp3 <- paste(corp, LL3, sep=":")
    eff.tt.final <- rcqp::cqi_subcorpus_size(sscorp3)
    cat(eff.tt.final, "\n")
    lst.freq <- rcqp::cqi_fdist1(sscorp3, "match", "lemma") # liste de fréquences du ss-corpus sélectionné
    cat("***-")
    lst.str <- rcqp::cqi_id2str(corp.lm, as.numeric(lst.freq[,1])) # on passe aux strings
    cat("***-")
    lst.freq <- as.data.frame(cbind(lst.freq, lst.str))
    lst.freq <- lst.freq[(lst.freq[,1] %in% cooc.final[,1]),] # on ne retient que les lemmes qui cooccurrent
    cat("***-", "\n")
    lst.freq <- lst.freq[(order(lst.freq[,3])),] # ordre alpha des lemmes
    cooc.final[,5] <- as.numeric(as.character(lst.freq[,2])) # on ajoute les fréquences globales au tableau
    A2 <- rcqp::cqi_dump_subcorpus(sscorp3)

  }

  #
  # filtrage des cooccurrents "significatifs" par divers coefficients-filtres
  #

  eff.piv <- rep(eff.pivot, nrow(cooc.final))
  cooc.final[,6] <- as.numeric(eff.piv)
  E11 <- as.numeric(cooc.final[,5])*as.numeric(cooc.final[,6])/as.numeric(eff.tt.final)
  cooc.final[,7] <- E11
  nams.tab <- c("ids.GD","ps.GD", "O11", "strings.cooc", "C1", "R1", "E11")
  names(cooc.final) <- nams.tab

  # coefficients de Dice-Soerensen, de Poisson-Stirling, de Daille, hypergéométrique, PMI, smin
  # (Stefan EVERT, _The Statistics of Word Cooccurrences_, 2005, pp. 75-94 // www.collocations.de)

  attach(as.data.frame(cooc.final),warn.conflicts=FALSE) # manoeuvre risquée...
  c.dice <- (O11*100000)/(C1+R1)^selec
  c.poiss <- O11*(log(O11)-log(E11)-1)
  c.daille <- log((O11^3)/E11)
  c.pmi <- log(O11/E11)
  c.hyperg <- as.numeric(as.character(round(-phyper(O11, C1, (eff.tt.final-C1), (2*dis*R1), lower.tail=FALSE, log.p=TRUE), digits=10)))
  Or <- O11/R1
  Oc <- O11/C1
  OCR <- cbind(Or,Oc)
  c.ms <- apply(OCR,1,min)
  coeff.tt <- cbind(c.dice,c.poiss,c.daille,c.pmi,c.hyperg,c.ms)

  detach()
  cooc.final <- cbind(cooc.final,coeff.tt)

  cat(paste(piv, pos.piv, eff.pivot, sep=" / "))
  cat("\n")

  #
  # mise en forme de l'affichage
  #

  if (coeff=="dice") numcol <- 8
  if (coeff=="poiss") numcol <- 9
  if (coeff=="daille") numcol <- 10
  if (coeff=="pmi") numcol <- 11
  if (coeff=="hyperg") numcol <- 12
  if (coeff=="mins") numcol <- 13
  resultat <- cooc.final[rev(order(cooc.final[,numcol])),]
  nb.lignes <- min(nbs, nrow(resultat))
  resultat <- resultat[1:nb.lignes,]
  ord.res <- order(resultat[,2], -resultat[,numcol])
  resultat <- resultat[ord.res,]

  cooc.affich.p <- as.data.frame(cbind(resultat[,2],resultat[,4],resultat[,3],resultat[,5],resultat[,numcol]))
  cooc.affich <- cooc.affich.p[1,]
  nm.affich <- c("POS", "lemme", "eff.cooc", "eff.lemme", coeff)
  cooc.affich <- rbind(nm.affich,rep(" ", 5), cooc.affich)

  for (i in 2:nrow(cooc.affich.p)) {
    if (cooc.affich.p[i,1]==cooc.affich.p[i-1,1]) cooc.affich <- rbind(cooc.affich, cooc.affich.p[i,])
    if (cooc.affich.p[i,1]!=cooc.affich.p[i-1,1]) cooc.affich <- rbind(cooc.affich, rep(" ",5))
  }
  MASS::write.matrix(cooc.affich,file="")

  # création de la liste finale / résultat
  # bizarreries : noms de variables d'une version antérieure, conservés pour correspondance avec les programmes suivants

  cooc0 <- data.frame(id=cooc.final[,1],pos=cooc.final[,2],cooc=cooc.final[,3], lem=cooc.final[,4],eff.lem=cooc.final[,5],eff.pivot=cooc.final[,6],Dice=cooc.final[,8],PoissonS=cooc.final[,9],Daille=cooc.final[,10],HyperG=cooc.final[,12],PMI=cooc.final[,11],cum="")
  cooc1 <- data.frame(id=resultat[,1],pos=resultat[,2],cooc=resultat[,3], lem=resultat[,4],eff.lem=resultat[,5],eff.pivot=resultat[,6],Dice=resultat[,8],PoissonS=resultat[,9],Daille=resultat[,10],HyperG=resultat[,12],PMI=resultat[,11],cum="")
  DF <- paste(D,F, sep="//")
  result <- list(call=match.call(),cooc0=cooc0,cooc1=cooc1,A2d3=A2,lgcorp=lgcorp,lemme=piv,cposLm=cpos.pivot,corpus=corp,effectif=eff.pivot,debfin=DF,pospiv=pos.piv,poscooc=poscooc,dist=dis,attr=attr,val=val,nb=nbs,NN=CORP.ln,tlex=ids.tablex)

  t2 <- Sys.time()
  td <- difftime(t1,t2, units="secs")
  cat("\n","temps écoulé :",round(td), "secondes \n")

  class(result) <- "COOCA"
  return (result)
}
