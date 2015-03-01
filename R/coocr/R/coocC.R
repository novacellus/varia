################################################################################

coocC <- function(x, y) {

  # mise en ordre de deux listes de cooccurrents
  # pour comparaison simplifiée

  if ( (!inherits(x, "COOCA")) | (!inherits(y, "COOCA")) )
    stop("en entrée : deux objets produits par coocA()", call.=F)

  gau <- x$cooc1[,c(3,4,7)]
  dro <- y$cooc1[,c(3,4,7)]
  gd <- merge(gau, dro, by.x=2, by.y=2, all.x=TRUE, all.y=TRUE, incomparables= 0)
  gd[is.na(gd)] <- 0
  gd <- gd[order(gd[,5]),]
  gdA <- gd[(gd[,4]==0),]
  gdA <- gdA[(order(gdA[,3])),]
  gd <- gd[(gd[,4]!=0),]
  gd <- gd[rev(order(gd[,3])),]
  gdB <- gd[(gd[,3]==0),]
  gdB <- gdB[rev(order(gdB[,5])),]
  gd <- gd[(gd[,3]!=0),]
  rA <- rank(gd[,3])
  rB <- rank(gd[,5])
  rdiff <- rB - rA
  gd <- gd[order(rdiff),]
  gdf <- rbind(gdA,gd,gdB)
  return(gdf)
}
