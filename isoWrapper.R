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


# Isotope error list
# I131 'Error in seq.default(ndk) : 'from' must be a finite number
#     In addition: Warning message:
#     In eval(ei, envir) : NAs introduced by coercion

# no 'you have reached a stable isotope!' output for:
# 89SR, 244CM
