listIsos <- scan('isotopeList.txt', character())

for (i in listIsos){
  iso <- toupper(i)
  existing <- dir('decayLists/')
  
  if (length(which(existing == toupper(iso))) > 0){
    next
  } else {
    message(paste('Getting data for', iso, sep = ' '))
    source('makeIsoRDS.R')
  }
}

