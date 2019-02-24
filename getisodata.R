library(stringr)
library(PeriodicTable)
library(plyr)
data("periodicTable")

Isotopes <- list()

isofile <- 'isotopes/JEFF33-rdd_all.asc'

#which(regexpr('213BI DECAY', readLines(con = 'isotopes/JEFF33-rdd_all.asc')) == 1)

linesplit <- function(x) unlist(strsplit(x, split = " "))[
  which(unlist(strsplit(x, split = " ")) != "")]


iso <- '225AC' # input the parent isotope! 
branchThreshold <- 0.0001

addIso <- list()
addIso$isotope <- iso
addIso$decayLevel <- 0
addIso$A <- as.numeric(str_extract(iso, "[0-9]+"))
addIso$symb <- (str_extract(iso, "[aA-zZ]+"))
addIso$symb <- periodicTable$symb[
  which(match(tolower(periodicTable$symb), tolower(addIso$symb)) == 1)]
addIso$Z = periodicTable$numb[
  which(match(periodicTable$symb, addIso$symb) == 1)]
addIso$masterYield <- 1


hfl <- linesplit(grep(iso, readLines(con = 'isotopes/JEFF33-rdd_all.asc'), value = TRUE)[1])

if (hfl[4] == 'Y'){
  t12 <- as.numeric(hfl[3])*365
} else if (hfl[4] == 'D'){
  t12 <- as.numeric(hfl[3])
} else if (hfl[4] == 'H') {
  t12 <- as.numeric(hfl[3])/24
} else if (hfl[4] == 'M') {
  t12 <- as.numeric(hfl[3])/24/60
} else if (hfl[4] == 'S') {
  t12 <- as.numeric(hfl[3])/24/60/60
} else if (hfl[4] == 'MS') {
  t12 <- as.numeric(hfl[3])/24/60/60/1000
}

addIso$t12 <- t12

linematches <- grep(iso, readLines(con = 'isotopes/JEFF33-rdd_all.asc'), value = FALSE)
# ndk_line <- grep(iso, readLines(con = 'isotopes/JEFF33-rdd_all.asc'), value = FALSE)[1]+1
ndk <- as.integer(linesplit(readLines(con = isofile)[linematches[1]+1])[5])

Decays <- list()
Alpha <- list()
Beta <- list()
Positron <- list()
EC <- list()
IT <- list()
for (i in seq(ndk)) {
  outp <- linesplit(readLines(con = 'isotopes/JEFF33-rdd_all.asc')[linematches[i+1]])
  outq <- linesplit(readLines(con = 'isotopes/JEFF33-rdd_all.asc')[linematches[i+1]+1])
  
  if (outp[2] == 'A') {
    Alpha$Q <- as.numeric(outp[6])
    Alpha$BranchYield <- as.numeric(outq[2])/100
    Alpha$BrQ <- Alpha$Q * Alpha$BranchYield
    Alpha$Daughter <- paste(addIso$A-4, periodicTable$symb[which(periodicTable$numb == addIso$Z-2)], sep = '')
  } 
  if (outp[2] == 'B-') {
    Beta$Q <- as.numeric(outp[6])
    Beta$BranchYield <- as.numeric(outq[2])/100
    Beta$BrQ <- Beta$Q * Beta$BranchYield
    Beta$Daughter <- paste(addIso$A, periodicTable$symb[which(periodicTable$numb == addIso$Z+1)], sep = '')
  } 
  if (outp[2] == 'B+') {
    Positron$Q <- as.numeric(outp[6])
    Positron$BranchYield <- as.numeric(outq[2])/100
    Positron$BrQ <- Positron$Q * Positron$BranchYield
    Positron$Daughter <- paste(addIso$A, periodicTable$symb[which(periodicTable$numb == addIso$Z-1)], sep = '')
  } 
  if (outp[2] == 'EC') {
    EC$Q <- as.numeric(outp[6])
    EC$BranchYield <- as.numeric(outq[2])/100
    EC$BrQ <- EC$Q * EC$BranchYield
    EC$Daughter <- paste(addIso$A, periodicTable$symb[which(periodicTable$numb == addIso$Z-1)], sep = '')
  } 
  if (outp[2] == 'IT') {
    IT$Q <- as.numeric(outp[6])
    IT$BranchYield <- as.numeric(outq[2])/100
    IT$BrQ <- EC$Q * EC$BranchYield
    IT$Daughter <- paste(addIso$A, addiso$symb, sep = '')
  } 
}

Decays$Alpha <- Alpha
Decays$Beta <- Beta
Decays$Positron <- Positron
Decays$EC <- EC
Decays$IT <- IT
Decays <- Decays[-which(lapply(Decays, length) == 0)]
addIso$Decays <- Decays

Isotopes[[length(Isotopes)+1]] <- addIso
names(Isotopes)[length(Isotopes)] <- iso



#---- Starting to add daughter isotopes ----

dk_levs <- 9

declevs <- array()
for (i in seq(dk_levs)) {
  for (n in seq(length(Isotopes))) {
    declevs[n] <- Isotopes[[n]]$decayLevel
  }
  
  isos <- which(declevs == i-1)
  
  for (j in isos) {
    for (dk in seq(length(Isotopes[[j]]$Decays))) {
      iso <- toupper(Isotopes[[j]]$Decays[[dk]]$Daughter)
      
      addIso <- list()
      addIso$isotope <- iso
      addIso$decayLevel <- i
      addIso$A <- as.numeric(str_extract(iso, "[0-9]+"))
      addIso$symb <- (str_extract(iso, "[aA-zZ]+"))
      addIso$symb <- periodicTable$symb[
        which(match(tolower(periodicTable$symb), tolower(addIso$symb)) == 1)]
      addIso$Z = periodicTable$numb[
        which(match(periodicTable$symb, addIso$symb) == 1)]
      addIso$masterYield <- Isotopes[[j]]$Decays[[dk]]$BranchYield * Isotopes[[j]]$masterYield
      
      if (addIso$masterYield < branchThreshold) {
        paste('The ', addIso$A, addIso$symb, ' branch has been terminated due to insufficient yield',
              sep = '')
        next
      }
      
      hfl <- linesplit(grep(iso, readLines(con = 'isotopes/JEFF33-rdd_all.asc'), value = TRUE)[1])
      if (length(hfl) == 0) {
        print('You have reached a stable isotope!')
      }
      
      if (hfl[4] == 'Y'){
        t12 <- as.numeric(hfl[3])*365
      } else if (hfl[4] == 'D'){
        t12 <- as.numeric(hfl[3])
      } else if (hfl[4] == 'H') {
        t12 <- as.numeric(hfl[3])/24
      } else if (hfl[4] == 'M') {
        t12 <- as.numeric(hfl[3])/24/60
      } else if (hfl[4] == 'S') {
        t12 <- as.numeric(hfl[3])/24/60/60
      } else if (hfl[4] == 'MS') {
        t12 <- as.numeric(hfl[3])/24/60/60/1000
      }
      
      addIso$t12 <- t12
      
      linematches <- grep(iso, readLines(con = 'isotopes/JEFF33-rdd_all.asc'), value = FALSE)
      # ndk_line <- grep(iso, readLines(con = 'isotopes/JEFF33-rdd_all.asc'), value = FALSE)[1]+1
      ndk <- as.integer(linesplit(readLines(con = isofile)[linematches[1]+1])[5])
      
      Decays <- list()
      Alpha <- list()
      Beta <- list()
      Positron <- list()
      EC <- list()
      IT <- list()
      
      for (i in seq(ndk)) {
        outp <- linesplit(readLines(con = 'isotopes/JEFF33-rdd_all.asc')[linematches[i+1]])
        outq <- linesplit(readLines(con = 'isotopes/JEFF33-rdd_all.asc')[linematches[i+1]+1])
        
        if (outp[2] == 'A') {
          Alpha$Q <- as.numeric(outp[6])
          Alpha$BranchYield <- as.numeric(outq[2])/100
          Alpha$BrQ <- Alpha$Q * Alpha$BranchYield
          Alpha$Daughter <- paste(addIso$A-4, periodicTable$symb[which(periodicTable$numb == addIso$Z-2)], sep = '')
        } 
        if (outp[2] == 'B-') {
          Beta$Q <- as.numeric(outp[6])
          Beta$BranchYield <- as.numeric(outq[2])/100
          Beta$BrQ <- Beta$Q * Beta$BranchYield
          Beta$Daughter <- paste(addIso$A, periodicTable$symb[which(periodicTable$numb == addIso$Z+1)], sep = '')
        } 
        if (outp[2] == 'B+') {
          Positron$Q <- as.numeric(outp[6])
          Positron$BranchYield <- as.numeric(outq[2])/100
          Positron$BrQ <- Positron$Q * Positron$BranchYield
          Positron$Daughter <- paste(addIso$A, periodicTable$symb[which(periodicTable$numb == addIso$Z-1)], sep = '')
        } 
        if (outp[2] == 'EC') {
          EC$Q <- as.numeric(outp[6])
          EC$BranchYield <- as.numeric(outq[2])/100
          EC$BrQ <- EC$Q * EC$BranchYield
          EC$Daughter <- paste(addIso$A, periodicTable$symb[which(periodicTable$numb == addIso$Z-1)], sep = '')
        } 
        if (outp[2] == 'IT') {
          IT$Q <- as.numeric(outp[6])
          IT$BranchYield <- as.numeric(outq[2])/100
          IT$BrQ <- EC$Q * EC$BranchYield
          IT$Daughter <- paste(addIso$A, addiso$symb, sep = '')
        } 
      }
      
      Decays$Alpha <- Alpha
      Decays$Beta <- Beta
      Decays$Positron <- Positron
      Decays$EC <- EC
      Decays$IT <- IT
      Decays <- Decays[-which(lapply(Decays, length) == 0)]
      addIso$Decays <- Decays
      
      if (length(Isotopes) > 1) {
        if (any(match(names(Isotopes),addIso$isotope), na.rm = TRUE)) {
          Isotopes[[which(match(names(Isotopes),addIso$isotope) == 1)]]$masterYield <- 
            Isotopes[[which(match(names(Isotopes),addIso$isotope) == 1)]]$masterYield + addIso$masterYield
        } else {
          Isotopes[[length(Isotopes)+1]] <- addIso
          names(Isotopes)[length(Isotopes)] <- iso
        }
      } else {
        Isotopes[[length(Isotopes)+1]] <- addIso
        names(Isotopes)[length(Isotopes)] <- iso
      }
    }
  }  
}


Isotopes$`213BI`$Decays$Beta$Daughter
