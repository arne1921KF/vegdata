# Load taxonomic reference list
load.taxlist <- function(refl, reflist.type= 'Turboveg', verbose, ...) {
  reflist.type <- match.arg(reflist.type, c('Turboveg'))
  if(reflist.type == 'Turboveg') {
    tv_home <- tv.home()
    reflist <- paste(refl, ifelse(verbose,'.verbose',''), sep='')
    if(is.null(store(reflist))) { # load.species(refl=refl, verbose = verbose)
      dbf <- if(verbose) 'tax.dbf' else 'species.dbf'
      supportedReflists <- c('GermanSL 1.0', 'GermanSL 1.1', 'GermanSL 1.2', 'Czsk 0.1')
      supportedReflists <- c(supportedReflists, sub(' ', '', supportedReflists))
      supportedReflists <- c(supportedReflists, tolower(supportedReflists))
      reflist.path <- file.path(tv_home, 'Species', refl, dbf)

      if(file.access(reflist.path)) {
        if(refl %in% supportedReflists) {
          cat('\nTaxonomic list (',dbf, ') of reflist version', refl, 'not available.\n\n')
            storage <- file.path(tv_home, 'Species')
          tfile <- tempfile()
           if(grepl('GermanSL', refl)) {
            version <- paste("version", substr(refl, 10, nchar(refl)), sep = "")
            download.file(paste('http://geobot.botanik.uni-greifswald.de/download/GermanSL',version,'GermanSL.zip',sep='/'), tfile)
            unzip(tfile, exdir= storage)
           }
           if(grepl('Czsk', refl)) {
            version <- paste("version", substr(refl, 6, nchar(refl)), sep = "")
            download.file(paste('http://geobot.botanik.uni-greifswald.de/download/CZSK',version, 'Czsk.zip',sep='/'), tfile)
            unzip(tfile, exdir= storage)
           }
          }  else stop('\nTaxonomic list ', refl, ' not supported.\n')
      } else storage <- file.path(tv_home, 'Species')

      species <- read.dbf(file.path(storage, refl, dbf), as.is=TRUE)
      names(species) <- TCS.replace(names(species))
      species$TaxonName <- taxname.abbr(species$TaxonName)
      if(verbose) {
        species$TaxonConcept <- taxname.abbr(species$TaxonConcept)
        species$VernacularName <- iconv(species$VernacularName, from='UTF8', to='')
        species$AUTHOR <- iconv(species$AUTHOR, from='UTF8', to='')
      }  else {
        species$VernacularName <- iconv(species$VernacularName, from='CP850', to='')
        species$AUTHOR <- iconv(species$AUTHOR, from='CP850', to='')
      }
      if(refl %in% supportedReflists && verbose==FALSE) species <- species[,c('TaxonUsageID','LETTERCODE','TaxonName', 'VernacularName','SYNONYM', 'TaxonConceptID')] else {
        include <- !names(species) %in% c('SHORTNAME')
        species <- species[, include]
      }
       store(reflist, species)
    } else  species <- store(reflist)      
  } else stop('Only reflisttype Turboveg implemented until now')
  
  return(species)
}

store <- local({
  refList <- list()
  function(name, value) {
    if(missing(name)) return(refList)
    if(missing(value)) return(refList[[name]])
    refList[name] <<- list(value)
  }
})

# As dBase is an old DOS format, Umlaute  are  stored  using  a  different  code  table
#    (namely ASCII) than most modern unices (namely ANSI).
taxname.abbr <- function(x, hybrid=FALSE) {
#  loc <- Sys.getlocale(category='LC_CTYPE')
#  Sys.setlocale("LC_ALL","C")
#  print('Executing taxname.abbr ...')
    x <- sub('\ ag[.]', ' agg.', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ ssp[.]', ' subsp.', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ v[.]\ ', ' var. ', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ sv[.]\ ', ' subvar. ', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ Sec[.]\ ', ' sect. ', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ Ser[.]\ ', ' ser. ', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ Subs[.]\ ', ' subsect. ', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ spec[.]', ' species', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ s[.]l[.]', ' s. l.', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ s[.]str[.]', ' s. str.', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ s[.]\ str[.]', ' sensustricto', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ s[.]\ l[.]', ' sensulato', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ s[.]lat[.]', ' sensulato', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ s[.] lat[.]', ' sensulato', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ s[.]\ ', ' subsp. ', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ sensustricto', ' s. str.', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ sensulato', ' s. l.', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ f[.]\ ', ' fo. ', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ sp[.]', ' spec.', x, perl=TRUE, useBytes=TRUE)
    x <- sub('\ nothosubsp[.]' , '\ nothossp.', x, perl=TRUE, useBytes=TRUE)
    if(hybrid)
      x <- sub('\ x.' , '\ ', x, perl=TRUE, useBytes=TRUE)
    #  Sys.setlocale(category='LC_CTYPE', locale=loc)
   return(x)  
}


taxname.simplify <- function(x, genus=TRUE, epithet=TRUE) {
#    x <- 'Sëlixae calcarae subsp. holdae'
    x <- gsub('\U00EB', 'e', x, perl=TRUE, useBytes=TRUE)
    x <- gsub('\U00CF', 'i', x, perl=TRUE, useBytes=TRUE)
    x <- gsub('ii', 'i', x, perl=TRUE, useBytes=TRUE)
    x <- gsub('nn', 'n', x, perl=TRUE, useBytes=TRUE)
    x <- gsub('ph', 'p', x, perl=TRUE, useBytes=TRUE)
    x <- gsub('rh', 'h', x, perl=TRUE, useBytes=TRUE)
    x <- gsub('th', 't', x, perl=TRUE, useBytes=TRUE)
    x <- gsub('tt', 't', x, perl=TRUE, useBytes=TRUE)
    x <- gsub('y', 'i', x, perl=TRUE, useBytes=TRUE)
if(epithet) {
    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('ae\\b', '', substr(x, regexpr('\ ', x), nchar(x))), sep='')
    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('arum\\b', '', substr(x, regexpr('\ ', x), nchar(x))), sep='')
    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('ea\\b', '', substr(x, regexpr('\ ', x), nchar(x))), sep='')
    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('ei\\b', '', substr(x, regexpr('\ ', x), nchar(x))), sep='')
    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('eos\\b', '', substr(x, regexpr('\ ', x), nchar(x))), sep='')
    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('ia\\b', '', substr(x, regexpr('\ ', x), nchar(x))), sep='')
    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('ium\\b', '', substr(x, regexpr('\ ', x), nchar(x))), sep='')
    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('ius\\b', '', substr(x, regexpr('\ ', x), nchar(x))), sep='')
    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('orum\\b', '', substr(x, regexpr('\ ', x), nchar(x))), sep='')

    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('a\\b', '', substr(x, regexpr('\ ', x), nchar(x))), sep='')
    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('e\\b', '', substr(x, regexpr('\ ', x), nchar(x))), sep='')
    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('ens\\b', '', substr(x, regexpr('\ ', x), nchar(x))), sep='')
    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('es\\b', '', substr(x, regexpr('\ ', x), nchar(x))), sep='')
    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('i\\b', '', substr(x, regexpr('\ ', x), nchar(x))), sep='')
    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('is\\b', '', substr(x, regexpr('\ ', x), nchar(x))), sep='')
    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('on\\b', '', substr(x, regexpr('\ ', x), nchar(x))), sep='')
    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('um\\b', '', substr(x, regexpr('\ ', x), nchar(x))), sep='')
    x <- paste(substr(x, 1, regexpr('\ ', x)-1), gsub('us\\b', '', substr(x, regexpr('\ ', x), nchar(x))), sep='')
}
if(genus) {
    x <-	paste(sub('a$', '', substr(x, 1, regexpr('\ ', x)-1)), substr(x, regexpr('\ ', x), nchar(x)), sep='')
    x <-	paste(sub('as$', '', substr(x, 1, regexpr('\ ', x)-1)), substr(x, regexpr('\ ', x), nchar(x)), sep='')
    x <-	paste(sub('e$', '', substr(x, 1, regexpr('\ ', x)-1)), substr(x, regexpr('\ ', x), nchar(x)), sep='')
    x <-	paste(sub('es$', '', substr(x, 1, regexpr('\ ', x)-1)), substr(x, regexpr('\ ', x), nchar(x)), sep='')
    x <-	paste(sub('eus$', '', substr(x, 1, regexpr('\ ', x)-1)), substr(x, regexpr('\ ', x), nchar(x)), sep='')
    x <-	paste(sub('is$', '', substr(x, 1, regexpr('\ ', x)-1)), substr(x, regexpr('\ ', x), nchar(x)), sep='')
    x <-	paste(sub('on$', '', substr(x, 1, regexpr('\ ', x)-1)), substr(x, regexpr('\ ', x), nchar(x)), sep='')
    x <-	paste(sub('u$', '', substr(x, 1, regexpr('\ ', x)-1)), substr(x, regexpr('\ ', x), nchar(x)), sep='')
    x <-	paste(sub('um$', '', substr(x, 1, regexpr('\ ', x)-1)), substr(x, regexpr('\ ', x), nchar(x)), sep='')
    x <-	paste(sub('us$', '', substr(x, 1, regexpr('\ ', x)-1)), substr(x, regexpr('\ ', x), nchar(x)), sep='')
}
return(x)
}
# taxname.simplify(x, Gattungsendung=TRUE, Artendung=TRUE)


TCS.replace <- function(x) {
## Turboveg  
  x <- replace(x, x=='ABBREVIAT', 'TaxonName')
  x <- replace(x, x=='taxonName', 'TaxonName')
  x <- replace(x, x=='Taxon', 'TaxonName')
  x <- replace(x, x=='SPECIES_NR', 'TaxonUsageID')
  x <- replace(x, x=='VALID_NAME', 'TaxonConcept')
  x <- replace(x, x=='VALID_NR', 'TaxonConceptID')
  x <- replace(x, x=='AGG_NAME', 'IsChildTaxonOf')
  x <- replace(x, x=='AGG', 'IsChildTaxonOfID')
  x <- replace(x, x=='SECUNDUM', 'AccordingTo')
  x <- replace(x, x=='NATIVENAME', 'VernacularName')
  x <- replace(x, x=='RANG', 'TaxonRank')
  x <- replace(x, x=='CLASSIFICA', 'Classification')
## Florkart Germany (BfN lists)
  x <- replace(x, x=='TAXNAME', 'TaxonName')
  x <- replace(x, x=='sipnr', 'TaxonConceptID')
  x <- replace(x, x=='SIPNR', 'TaxonConceptID')
  x <- replace(x, x=='namnr', 'TaxonUsageID')
  x <- replace(x, x=='NAMNR', 'TaxonUsageID')
  x <- replace(x, x=='TAXNR', 'TaxonUsageID')
  x <- replace(x, x=='AGG_NAME', 'IsChildTaxonOf')
  x <- replace(x, x=='AGGNR', 'IsChildTaxonOfID')
  x <- replace(x, x=='RANK', 'TaxonRank')  
  x <- replace(x, x=='Rank', 'TaxonRank')  
## ESveg
  x <- replace(x, x=="taxonCode", 'TaxonUsageID')
  x <- replace(x, x=="observationCode", "RELEVE_NR")
  x <- replace(x, x=="ObservationID", "RELEVE_NR")
  x <- replace(x, x=="stratumCode", "LAYER")
  x <- replace(x, x=="Stratum", "LAYER")
  x <- replace(x, x=="Percentage_mean", "COVER_PERC")
  x <- replace(x, x=="coverPercent", "COVER_PERC")
## CDM
  x <- replace(x, x=='Taxon', 'TaxonName')
  x <- replace(x, x=='Taxon.ID', 'TaxonUsageID')  
  return(x)
}



#########################################################
## function tax
#########################################################

"tax" <- function(...) UseMethod("tax")

tax.default <- function(x, refl, verbose = FALSE, syn = TRUE, concept = NULL, strict = FALSE, vernacular = FALSE, simplify=FALSE, quiet = FALSE, ...) {
	tv_home <- tv.home()
#########################################################
## internal functions
#########################################################

concept.FUN <- function(species, concept, dbf, ...) {
  cat('\n Will use taxon concept', concept, '.\n\n')
  verbose=TRUE
  species$TaxonName <- as.character(species$TaxonName)
  species$TaxonConcept <- as.character(species$TaxonConcept)
  species$IsChildTaxonOf <- as.character(species$IsChildTaxonOf)
  species$AccordingTo <- as.character(species$AccordingTo)
  conc <- read.dbf(file.path(tv_home, 'Species', refl, paste(concept,'dbf',sep='.')), as.is=TRUE)
  co <- conc[match(species$TaxonUsageID, conc$TaxonUsageID, nomatch = 0),]
  species[match(conc$TaxonUsageID,species$TaxonUsageID),c('SYNONYM','TaxonConceptID','IsChildTaxonOfID')] <- co[match(conc$TaxonUsageID,co$TaxonUsageID),c('SYNONYM','TaxonConceptID','IsChildTaxonOfID')]
  species$TaxonName[match(conc$TaxonUsageID,species$TaxonUsageID,nomatch = 0)] <- co$TaxonName[match(conc$TaxonUsageID,co$TaxonUsageID,nomatch = 0)]
  species$TaxonConcept[match(conc$TaxonUsageID,species$TaxonUsageID,nomatch = 0)] <- co$TaxonConcept[match(conc$TaxonUsageID,co$TaxonUsageID,nomatch = 0)]
  species$TaxonRank[match(conc$TaxonUsageID,species$TaxonUsageID,nomatch = 0)] <- co$TaxonRank[match(conc$TaxonUsageID,co$TaxonUsageID,nomatch = 0)]
  species$IsChildTaxonOf[match(conc$TaxonUsageID,species$TaxonUsageID,nomatch = 0)] <- co$IsChildTaxonOf[match(conc$TaxonUsageID,co$TaxonUsageID,nomatch = 0)]
  species$AccordingTo[match(conc$TaxonUsageID,species$TaxonUsageID,nomatch = 0)] <- co$AccordingTo[match(conc$TaxonUsageID,co$TaxonUsageID,nomatch = 0)]
}

# Subsetting
select.taxa <- function(x, species, strict, vernacular = FALSE, simplify = FALSE, genus = TRUE, epithet = TRUE, ...) {
  if(is.factor(x)) x <- as.character(x)
  if(simplify) strict <- TRUE
  if(is.numeric(x) | is.integer(x))
    l <- species[match(x, species$TaxonUsageID),]  ## Tax numbers
  if(is.character(x)) {
    if(nchar(x[1]) == 7 & x[1] == toupper(x[1]))  {  ## Lettercode
      l <- species[species$LETTERCODE %in% x,] 
    } else {		## Taxnames
    	l <- species[species$TaxonUsageID==(-1),]
      for(i in 1:length(x))
        if(vernacular) {
          l <- if(!strict) rbind(l, species[grep(x[i], species$VernacularName, useBytes=TRUE), ]) else 
            rbind(l, species[match(species$VernacularName, x[i], nomatch = 0) > 0, ])
        }
      else {
        l <- if(!strict) 
        	# if(simplify) rbind(l, species[grep(x[i], species$TaxonName, useBytes=TRUE), ])
        		rbind(l, species[grep(x[i], species$TaxonName, useBytes=TRUE), ]) else 
        	if(simplify) {
        		# print(taxname.simplify(species$TaxonName, genus=TRUE))
        		rbind(l, species[match(taxname.simplify(species$TaxonName, genus, epithet, ...), taxname.simplify(x[i], genus, epithet, ...), nomatch = 0) > 0, ]) 
        		} else
          rbind(l, species[match(species$TaxonName, x[i], nomatch = 0) > 0, ])
      }}
  }
  if(length(l) == 0) stop('No species found!')
  return(l)
}

##### end of tax internal functions #####

if(missing(refl)) refl <- tv.refl(refl, tv_home=tv_home)
if(!quiet) cat('Reference list used:', refl, '\n')	
species <- load.taxlist(refl, reflist.type='Turboveg', verbose=verbose)

### Filter
# if(!is.null(concept)) species <- concept.FUN(species, concept)
if(x[1] != 'all') species <- select.taxa(x, species, strict, vernacular, simplify, ...)
  if(!syn) species <- species[species$SYNONYM == FALSE,]
if(simplify) species$simpleSpeciesName <- taxname.simplify(species$TaxonName, ...)
return(species)
}
#  ls(pos='package:vegdata')

# getS3method('tax', 'default')
tax.veg <- function(veg, ...) {
  if(is.null(attr(veg, 'taxreflist'))) stop('Object must have attribute taxreflist.')
  taxa <- tax.default(names(veg), syn=FALSE, ...)
  cat('Taxa in vegetation matrix:\n\n')
  return(taxa$TaxonName[match(names(veg), taxa$LETTERCODE)])
}

