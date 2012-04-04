tv.eco <- function () paste('Deprecated, use tv.traits instead.')
  
tv.traits <- function (db, eco = 'ecodbase.dbf', refl)
{
    tv_home <- tv.home()
    if(missing(refl))  refl <- if(missing(db)) tv.refl() else tv.refl(db)
    ecodb <- read.dbf(file.path(tv_home, 'Species', refl, eco))
    empty <- function(x) all(is.na(x) | x == 0)
    na <- apply(ecodb, 2, empty)
    if(any(na)) {cat("\n The following columns contain no data and are omitted: \n")
                cat(names(ecodb)[na])
                ecodb <- ecodb[, !na]
                }
    cat("\n\n Changing character fields into logical, integer or numericals if appropriate: \n")
# ecoDB <- apply(ecodb, 2, function(x) type.convert(as.character(x)))
# doesnt work 
    ecoDB <- ecodb
    for(i in 1:ncol(ecodb)) if(is.factor(ecodb[,i])) {
      ecoDB[,i] <- as.character(ecodb[,i])
      ecoDB[,i] <- type.convert(ecoDB[,i]) }

    for(i in 1:ncol(ecoDB))  if(class(ecodb[,i]) != class(ecoDB[,i])) cat('Class of', names(ecoDB)[i], 'changed to ', class(ecoDB[,i]), '\n')
 
#     if(rm.dupl) {
#       taxa <- tax('all', refl=refl)
#       ecoDB$LETTERCODE <- taxa$LETTERCODE[match(ecoDB$TAXNR, taxa$SPECIES_NR)]
#       ecoDB$SYNONYM <- taxa$SYNONYM[match(ecoDB$TAXNR, taxa$SPECIES_NR)]
#       tab <- table(ecoDB$LETTERCODE)
#       ind <- ecoDB$LETTERCODE %in% names(tab)[tab>1] & ecoDB$SYNONYM
#       ecoDB <- ecoDB[!ind,]
#       rownames(ecoDB) <- ecoDB$LETTERCODE
#     }
    ecoDB
}
