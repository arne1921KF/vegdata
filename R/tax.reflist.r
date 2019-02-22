store <- local({
  refList <- list()
  function(name, value) {
    if(missing(name)) return(refList)
    if(missing(value)) return(refList[[name]])
    refList[name] <<- list(value)
  }
})


# Load taxonomic reference list
load.taxlist <- function(refl, reflist.type = c('Turboveg', 'EDIT'), detailed = FALSE, check = FALSE, hybrid, ...) {
  reflist.type <- match.arg(reflist.type, c('Turboveg', 'EDIT'))
  if(missing(hybrid)) hybrid = 'substitute'
  if(reflist.type == 'Turboveg') {
    args <- list(check = check)
    tv_home <- do.call("tv.home", args)
    if(detailed) dbf <- 'tax.dbf' else dbf <- 'species.dbf'
    reflist.path <- file.path(tv_home, 'Species', refl, dbf)
#    print(reflist.path)
    reflist <- paste(refl, ifelse(detailed,'.detailed',''), sep='')
    if(is.null(store(reflist))) { # load.species(refl=refl, detailed = detailed)
      supportedReflists <- c('GermanSL 1.0', 'GermanSL 1.1', 'GermanSL 1.2', 'GermanSL 1.3', 'GermanSL 1.4')
      supportedReflists <- c(supportedReflists, sub(' ', '', supportedReflists))
      supportedReflists <- c(supportedReflists, tolower(supportedReflists))
      if(!file.exists(reflist.path)) {
        message(paste('Taxonomic reference list file', reflist.path, 'does not exist.'))
        if(refl %in% supportedReflists) {
          message('\nTaxonomic list (',dbf, ') of reflist', refl, ' not available.\n\n')
          tfile <- tempfile()
           if(grepl('GermanSL', refl)) {
            version <- substr(refl, 10, nchar(refl))
            m <- try(download.file(paste('https://germansl.infinitenature.org/GermanSL', version, 'GermanSL.zip',sep='/'), tfile), silent=TRUE)
            if(m == 0) unzip(tfile, exdir= file.path(tv_home, 'Species')) else 
              unzip(file.path(path.package('vegdata'), 'tvdata','Species','TaxrefExample.zip'), exdir = file.path(tv_home, 'Species'))
           }
          } else message('\nTaxonomic list ', refl, ' not supported for automatic download.\n')
      }
      
      if(file.exists(reflist.path)) species <- read.dbf(reflist.path, as.is=TRUE) else stop(paste('Reference list file', reflist.path, 'does not exist.'))
      names(species) <- TCS.replace(names(species)) # replaces several column names, depending on list source. source: tax.names.R
      species$TaxonName <- taxname.abbr(..., x = species$TaxonName, hybrid = hybrid) # replaces several characters and abbreviations. source: tax.names.R 
      if(detailed) { # i.e. if tax.dbf is loaded
        species$TaxonConcept <- taxname.abbr( ..., x=species$TaxonConcept, hybrid = hybrid)
        if('VernacularName' %in% names(species) & Sys.info()['sysname'] != 'SunOS') species$VernacularName <- iconv(species$VernacularName, from='CP437', to='UTF-8') # possible also: to = '' on windows locale 1252, ISO-8859-1
        if('NameAuthor' %in% names(species) & Sys.info()['sysname'] != 'SunOS') species$NameAuthor <- iconv(species$NameAuthor, from='CP437', to='UTF-8')
        if('AccordingTo' %in% names(species) & Sys.info()['sysname'] != 'SunOS') species$AccordingTo <- iconv(species$AccordingTo, from='CP437', to='UTF-8') # NB: col is broken in v.1.3 tax.dbf data source!
        if('TaxonName' %in% names(species) & Sys.info()['sysname'] != 'SunOS') species$TaxonName <- iconv(species$TaxonName, from='UTF-8', to='UTF-8') # of course this is col is present!
        if('NACHWEIS' %in% names(species) & Sys.info()['sysname'] != 'SunOS') species$NACHWEIS <- iconv(species$NACHWEIS, from='UTF-8', to='UTF-8') # 
        if('BEGRUEND' %in% names(species) & Sys.info()['sysname'] != 'SunOS') species$BEGRUEND <- iconv(species$BEGRUEND, from='UTF-8', to='UTF-8') #
        if('TaxonConcept' %in% names(species) & Sys.info()['sysname'] != 'SunOS') species$TaxonConcept <- iconv(species$TaxonConcept, from='UTF-8', to='UTF-8') # cf. Persicaria × condensata
          if(refl == "GermanSL 1.3") {
          # load germanSL 1.4, replace AccordingTo with AccordingTo (or SECUNDUM, respectively)
          # TO BE RUN OUTSIDE OF FUNCTION IN PRODUCTIVE CODE   
            # secundum <- load.taxlist("GermanSL 1.4", detailed=TRUE) %>% #using the most recent German SL version here!
            # select(TaxonUsageID, TaxonConceptID, AccordingTo) 
            # tax_dbf <- tax_dbf %>%
            #   select(-AccordingTo) %>% 
            #   left_join(secundum, by=c("TaxonUsageID", "TaxonConceptID"))
            # rm(secundum)
          } 
         }  else { # i.e. if species.dbf is loaded - legacy code?
        if('VernacularName' %in% names(species)& Sys.info()['sysname'] != 'SunOS') species$VernacularName <- iconv(species$VernacularName, from='CP437', to='UTF-8') # alternatively from = ''?
        if('TaxonName' %in% names(species) & Sys.info()['sysname'] != 'SunOS') species$TaxonName <- iconv(species$TaxonName, from='UTF-8', to='UTF-8')
        }
      if(refl %in% supportedReflists && detailed==FALSE) species <- species[,c('TaxonUsageID','LETTERCODE','TaxonName', 'VernacularName','SYNONYM', 'TaxonConceptID')] else {
        include <- !names(species) %in% c('SHORTNAME')
        species <- species[, include]
      }
       store(reflist, species)
    } else  species <- store(reflist)
  } else stop(c('Until now only reference list type "Turboveg" is supported. If you want to use taxval with other reflists please contact me.'))
  	
   return(species)
}
# Note: Since AccordingTo it broken in GermanSL 1.3, it should be replaced with data from GermanSL 1.4. This means loading the v1.4 as well!

tv.refl <- function(refl, db, tv_home) {
#   capwords <- function(s, strict = FALSE) {
#       cap <- function(s) paste(toupper(substring(s,1,1)), {s <- substring(s,2); if(strict) toupper(s) else s}, sep = "", collapse = " " )
#       sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
#   }
  if(missing(tv_home)) tv_home <- tv.home()
  if(!missing(db)) {
    if(file.access(file.path(tv_home, 'Data', db[1], 'tvwin.dbf')) == 0) {
      refl <- read.dbf(file.path(tv_home, 'Data', db[1], 'tvwin.dbf'), as.is = TRUE)$FLORA
    } else {
      dbattr <- file.path(tv_home, 'Data', db[1], 'tvwin.set')
      if(file.access(dbattr) == 0) refl <-  sub('A\002', '', readBin(dbattr,what='character', n=3)[3]) else 
    stop('Database attribute file tvwin.set from database "', db[1], '" not available. Please specify name of taxonomic reference list!') 
  } } else  
    if(!missing(refl)) {
      rli <- list.dirs(path = file.path(tv_home, "Species"), full.names = TRUE, recursive = FALSE)
      rli <- sapply(rli, function(x) substring(x, nchar(tv_home) + 10), USE.NAMES = FALSE)
      if(length(rli) > 0) refl <- match.arg(refl, rli)
    } else refl <- 'GermanSL 1.3'
#  if(!exists(refl)) refl <- fun(gsub(' ','',refl))
  if(tolower(substr(refl, 1,8)) == 'germansl') refl <- paste('GermanSL', substring(refl, 9, nchar(refl)), sep='')
  return(refl)
 }
