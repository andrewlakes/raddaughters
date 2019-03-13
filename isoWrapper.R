listIsos <- scan('isotopeList.txt', character())

for (i in listIsos){
  iso <- toupper(i)
  existing <- dir('decayLists/')
  
  if (length(which(existing == toupper(iso))) > 0){
    next
  } else {
    print(iso)
    source('makeIsoRDS.R')
  }
}