# Different checks

#' Wrapper function listing corpora available on the system
.list_corpora <- function() {
  as.vector( rcqp::cqi_list_corpora() )
}
#' Wrapper function listing attributes
#' @return List of positional and structural attributes defined on the corpus
.list_attrs <- function(corpus) {
  list(p=rcqp::cqi_attributes(corpus,"p"),s=rcqp::cqi_attributes(corpus,"s"))
}
#' Checks corpus name submitted by user.
#' @param corpus_name Corpus name
#' @return If corpus exists, returns TRUE; if not, returns FALSE.
.check_corpus_name <- function(corpus_name) {
  list_corpus <- .list_corpora()
  if (!corpus_name %in% list_corpus)
    return(FALSE)
  TRUE
}

#' Checks if word exists in the corpus
#' @return If word exists in the corpus, returns its id; if not, returns useful errors
find_word <- function(corp,corp_lemma,piv) {
  corp.lm <-paste(corp, ".", corp_lemma , sep="")
  id.piv <- rcqp::cqi_str2id(corp.lm, piv)
  if (id.piv < 0) {  # si lemme inconnu, recherche simple par joker affixé
    id.prov <- rcqp::cqi_regex2id(corp.lm, piv) # id provisoire
    if (length(id.prov)==0){
      piv2 <- paste(piv, ".*", sep="")    # on ajoute un joker
      id.prov <- rcqp::cqi_regex2id(corp.lm, piv2) # vecteur des id correspondant à piv+joker
    }
    if (length(id.prov)==0) stop(paste(piv, ": lemme inconnu ! utilisez la fonction flou() !"),call.=FALSE)
    str.prov <- rcqp::cqi_id2str(corp.lm, id.prov)
    cpos.prov <- NULL
    pos.prov <- NULL
    freq.prov <- NULL
    for (i in 1:length(id.prov)){
      cpos.prov[i] <- rcqp::cqi_id2cpos(corp.lm, id.prov[i])
      pos.prov[i] <- rcqp::cqi_cpos2str(corp.ps, cpos.prov[i])
      freq.prov[i] <- rcqp::cqi_id2freq(corp.lm, id.prov[i])
      cat(str.prov[i], "/", pos.prov[i], "/", freq.prov[i], "\n")
    }
    if (length(id.prov) > 0) stop ("Choisissez parmi ces lemmes !", call.=FALSE)
  }
  return(id.piv)
}
#' Find all occurrences of word by its id
#' @return Vector of word corpus positions
find_words <-function(piv.id,attr="",val="",...){
  dots <-list(...)
  if ( attr == "" && val =="") {
    dump.pivot <-rcqp::cqi_id2cpos(dots$corp.lm,piv.id)
  } else if( !(attr == "" || val == "" ) ) {
      cat("La recherche du pivot dans une fraction du corpus total prend du temps (beaucoup)... \n")
      LL2 <- "Lsub2"
      quer2 <- paste("[lemma='", piv, "' %cd]", "::match.", attr, "=\"", val, "\"", sep="")
      rcqp::cqi_query(corp, LL2, quer2)
      sscorp2 <- paste(corp, LL2, sep = ":")
      if (rcqp::cqi_subcorpus_size(sscorp2)==0) stop ("Pas de réponse, vérifiez la requête !", call.=F)
      dump.pivot <- rcqp::cqi_dump_subcorpus(sscorp2)[,1]
      rcqp::cqi_drop_subcorpus(sscorp2)
  } else if ( attr  == "" || val == "" ) {
      stop ("Requête incomplète !", call.=FALSE)
      cat ("\n Vérifications simples : OK ! \n")
  }
  return(dump.pivot)
}
