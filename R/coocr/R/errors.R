#' Error handling
#' @type Error type
#'
.error_message <- function (type, ...) {
  dots <- list(...)

  if(type == "no_such_corpus") {
    corpora <- .list_corpora()
    msg <-paste (
        dots$corpus, ": aucun corpus de ce nom !", "\n",
        "Vous pouvez choisir entre:", "\n",
        paste(corpora,collapse = "\n"),
      sep="")
  }
  if(type == "no_p_attr") {
    msg <-paste(
      dots$corp_lemma,dots$corp_pos,
      ": no such positional attribute","\n",
      "Definiez l'attribute",
      if( length(dots$corp_lemma) > 0 ) {"corp_lemma"}  else {dots$corp_pos},
      "comme un de suivants:","\n",
      paste("-",dots$p_attr_list,collapse="\n"),
      sep="" )
  }
  msg
}
