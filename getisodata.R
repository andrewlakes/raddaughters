library(stringr)
library(PeriodicTable)
library(plyr)
data("periodicTable")


Isotopes <- list()

isofile <- 'isotopes/JEFF33-rdd_all.asc' # this is the master isotope data file

# custom function for splitting lines at whitespaces and unlisting
linesplit <- function(x) unlist(strsplit(x, split = " "))[
  which(unlist(strsplit(x, split = " ")) != "")]

iso <- '227TH' # input the parent isotope! 
branchThreshold <- 0.0001 # input the branch percentage threshold for abandoning a branch

# start making the new isotope in the addIso list
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

# get the half-life data, check for units and convert to days 
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

# check for which lines match the decay information for the isotope of interest
linematches <- grep(iso, readLines(con = 'isotopes/JEFF33-rdd_all.asc'), value = FALSE)

# number of decays as read from the data file
ndk <- as.integer(linesplit(readLines(con = isofile)[linematches[1]+1])[5])

Decays <- list()
Alpha <- list()
Beta <- list()
Positron <- list()
EC <- list()
IT <- list()

# search through the data file and get the lines that match the isotope, find the decay type(s)
for (i in seq(ndk)) {
  outp <- linesplit(readLines(con = 'isotopes/JEFF33-rdd_all.asc')[linematches[i+1]])
  outq <- linesplit(readLines(con = 'isotopes/JEFF33-rdd_all.asc')[linematches[i+1]+1])
  
  if (outp[2] == 'A') {
    Alpha$Q <- as.numeric(outp[6])
    Alpha$branchYield <- as.numeric(outq[2])/100
    Alpha$BrQ <- Alpha$Q * Alpha$branchYield
    Alpha$daughter <- paste(addIso$A-4, periodicTable$symb[which(periodicTable$numb == addIso$Z-2)], sep = '')
  } 
  if (outp[2] == 'B-') {
    Beta$Q <- as.numeric(outp[6])
    Beta$branchYield <- as.numeric(outq[2])/100
    Beta$BrQ <- Beta$Q * Beta$branchYield
    Beta$daughter <- paste(addIso$A, periodicTable$symb[which(periodicTable$numb == addIso$Z+1)], sep = '')
  } 
  if (outp[2] == 'B+') {
    Positron$Q <- as.numeric(outp[6])
    Positron$branchYield <- as.numeric(outq[2])/100
    Positron$BrQ <- Positron$Q * Positron$branchYield
    Positron$daughter <- paste(addIso$A, periodicTable$symb[which(periodicTable$numb == addIso$Z-1)], sep = '')
  } 
  if (outp[2] == 'EC') {
    EC$Q <- as.numeric(outp[6])
    EC$branchYield <- as.numeric(outq[2])/100
    EC$BrQ <- EC$Q * EC$branchYield
    EC$daughter <- paste(addIso$A, periodicTable$symb[which(periodicTable$numb == addIso$Z-1)], sep = '')
  } 
  if (outp[2] == 'IT') {
    IT$Q <- as.numeric(outp[6])
    IT$branchYield <- as.numeric(outq[2])/100
    IT$BrQ <- EC$Q * EC$branchYield
    IT$daughter <- paste(addIso$A, addiso$symb, sep = '')
  } 
}

# populate the list Decays, and remove the empty ones
Decays$Alpha <- Alpha
Decays$Beta <- Beta
Decays$Positron <- Positron
Decays$EC <- EC
Decays$IT <- IT
Decays <- Decays[-which(lapply(Decays, length) == 0)]
addIso$Decays <- Decays

# add the new isotope to the Isotopes list
Isotopes[[length(Isotopes)+1]] <- addIso
names(Isotopes)[length(Isotopes)] <- iso


#---- Starting to add daughter isotopes ----

dk_levs <- 9 # choose the number of levels of the decay chain to follow

declevs <- array()

# make an array of isotopes at a given decay level (i.e. no repeating)
for (i in seq(dk_levs)) {
  for (n in seq(length(Isotopes))) {
    declevs[n] <- Isotopes[[n]]$decayLevel
  }
  
  isos <- which(declevs == i-1)
  
  # repeat loop above
  for (j in isos) {
    for (dk in seq(length(Isotopes[[j]]$Decays))) {
      iso <- toupper(Isotopes[[j]]$Decays[[dk]]$daughter)
      
      addIso <- list()
      addIso$isotope <- iso
      addIso$decayLevel <- i
      addIso$A <- as.numeric(str_extract(iso, "[0-9]+"))
      addIso$symb <- (str_extract(iso, "[aA-zZ]+"))
      addIso$symb <- periodicTable$symb[
        which(match(tolower(periodicTable$symb), tolower(addIso$symb)) == 1)]
      addIso$Z = periodicTable$numb[
        which(match(periodicTable$symb, addIso$symb) == 1)]
      addIso$masterYield <- Isotopes[[j]]$Decays[[dk]]$branchYield * Isotopes[[j]]$masterYield
      
      # check to see if the branch is insignificant
      if (addIso$masterYield < branchThreshold) {
        print(paste('The ', addIso$A, addIso$symb, ' branch has been terminated due to insufficient yield',
              sep = ''))
        next
      }
      
      # check to see if there is a half-life!  If not, the daughter isotope is stable by JEFF definitions
      hfl <- linesplit(grep(iso, readLines(con = 'isotopes/JEFF33-rdd_all.asc'), value = TRUE)[1])
      if (length(hfl) == 0) {
        stop('You have reached a stable isotope!')
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
          Alpha$branchYield <- as.numeric(outq[2])/100
          Alpha$BrQ <- Alpha$Q * Alpha$branchYield
          Alpha$daughter <- paste(addIso$A-4, periodicTable$symb[which(periodicTable$numb == addIso$Z-2)], sep = '')
        } 
        if (outp[2] == 'B-') {
          Beta$Q <- as.numeric(outp[6])
          Beta$branchYield <- as.numeric(outq[2])/100
          Beta$BrQ <- Beta$Q * Beta$branchYield
          Beta$daughter <- paste(addIso$A, periodicTable$symb[which(periodicTable$numb == addIso$Z+1)], sep = '')
        } 
        if (outp[2] == 'B+') {
          Positron$Q <- as.numeric(outp[6])
          Positron$branchYield <- as.numeric(outq[2])/100
          Positron$BrQ <- Positron$Q * Positron$branchYield
          Positron$daughter <- paste(addIso$A, periodicTable$symb[which(periodicTable$numb == addIso$Z-1)], sep = '')
        } 
        if (outp[2] == 'EC') {
          EC$Q <- as.numeric(outp[6])
          EC$branchYield <- as.numeric(outq[2])/100
          EC$BrQ <- EC$Q * EC$branchYield
          EC$daughter <- paste(addIso$A, periodicTable$symb[which(periodicTable$numb == addIso$Z-1)], sep = '')
        } 
        if (outp[2] == 'IT') {
          IT$Q <- as.numeric(outp[6])
          IT$branchYield <- as.numeric(outq[2])/100
          IT$BrQ <- EC$Q * EC$branchYield
          IT$daughter <- paste(addIso$A, addiso$symb, sep = '')
        } 
      }
      
      Decays$Alpha <- Alpha
      Decays$Beta <- Beta
      Decays$Positron <- Positron
      Decays$EC <- EC
      Decays$IT <- IT
      Decays <- Decays[-which(lapply(Decays, length) == 0)]
      addIso$Decays <- Decays
      
      # check to see if two branches are converging, and if so add the masterYield parameters
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


Isotopes$`213BI`$Decays$Beta$daughter
