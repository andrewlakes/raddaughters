library(stringr)
library(PeriodicTable)
library(plyr)
data("periodicTable")

Isotopes <- list()

isofile <- 'isotopes/JEFF33-rdd_all.asc'

#which(regexpr('213BI DECAY', readLines(con = 'isotopes/JEFF33-rdd_all.asc')) == 1)

linesplit <- function(x) unlist(strsplit(x, split = " "))[
  which(unlist(strsplit(x, split = " ")) != "")]


iso <- '213BI'

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

dk_levs <- 1

isos <- as.numeric(laply(Isotopes, function(x) which(x$decayLevel == dk_levs-1)))


for (dk in seq(length(Isotopes[[isos]]$Decays))) {
  iso <- toupper(Isotopes[[isos]]$Decays[[dk]]$Daughter)
  
  addIso <- list()
  addIso$isotope <- iso
  addIso$decayLevel <- dk_levs
  addIso$A <- as.numeric(str_extract(iso, "[0-9]+"))
  addIso$symb <- (str_extract(iso, "[aA-zZ]+"))
  addIso$symb <- periodicTable$symb[
    which(match(tolower(periodicTable$symb), tolower(addIso$symb)) == 1)]
  addIso$Z = periodicTable$numb[
    which(match(periodicTable$symb, addIso$symb) == 1)]
  addIso$masterYield <- Isotopes[[isos]]$Decays[[dk]]$BranchYield * Isotopes[[isos]]$masterYield
  
  
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
  
}   










addIso$Decays$Beta$Daughter

Isotopes$`213BI`$Decays$Beta$Daughter


laply(Isotopes, function(x) which(x$decayLevel == 0))



Isotopes


periodicTable

paste(addIso$A-4, periodicTable$symb[which(periodicTable$numb == addIso$Z-2)], sep = '')
paste(addIso$A, periodicTable$symb[which(periodicTable$numb == addIso$Z+1)], sep = '')







outp <- linesplit(readLines(con = 'isotopes/JEFF33-rdd_all.asc')[linematches[i+1]])
outq <- linesplit(readLines(con = 'isotopes/JEFF33-rdd_all.asc')[linematches[i+1]+1])

if (outp[2] == 'A') {
  Alpha$E <- as.numeric(outp[6])
  Alpha$BranchYield <- as.numeric(outq[2])/100
} 
if (outp[2] == 'B-') {
  Beta$E <- as.numeric(outp[6])
  Beta$BranchYield <- as.numeric(outq[2])/100
} 
if (outp[2] == 'EC') {
  EC$E <- as.numeric(outp[6])
  EC$BranchYield <- as.numeric(outq[2])/100
} 













else if (dk[i] == 'B-') { 
  linenum <- which(regexpr('Mean B-', readLines(con = paste('isotopes/',isotope, sep = ''))) == 1)
  line <- unlist(strsplit(readLines(con = paste('isotopes/',isotope, sep = ''))[linenum], split = ' '))
  Beta$Mean <- as.numeric(line[which(line != "")][4])
  Beta$SD <- as.numeric(line[which(line != "")][6])
} 
outp[2]












grep(iso, readLines(con = 'isotopes/JEFF33-rdd_all.asc'), value = FALSE)






grep(, readLines(con = 'isotopes/JEFF33-rdd_all.asc'))

Decays <- list('Alpha' = Alpha, 'Beta' = Beta, 'Gamma' = Gamma, 'Xray' = Xray, 'CE' = CE)


Decays$Gamma$Mean










which(regexpr('Mean Gamma', readLines(con = paste('isotopes/',isotope, sep = ''))) == 1)

startE <- which(regexpr('Mean Gamma', readLines(con = paste('isotopes/',isotope, sep = ''))) == 1)
endE <- which(regexpr('Mean Recoil', readLines(con = paste('isotopes/',isotope, sep = ''))) == 1)
energies <- readLines(con = paste('isotopes/',isotope, sep = ''))[startE:endE]

spacesplit <- function(x) unlist(strsplit(x, split = " "))
es <- sapply(X = energies, FUN = spacesplit)
es <- es[which(es != "")]

Alpha <- list()
Alpha$Mean <- 5



Isotopes$Ac225$Decays$Alpha$MeanE




a <- unlist(strsplit(energies[which(regexpr('Alpha', energies) != -1)], split = ' '))
alpha$mean <- as.numeric(a[which(a != '')][4])
alpha$sd <- as.numeric(a)

Emean <- array(as.numeric(a[which(a != '')][4]))
Esd <- array(as.numeric(a[which(a != '')][6]))



  which(unlist(strsplit(energies[which(regexpr('Alpha', energies) != -1)], split = ' ')) != "")



addiso <- list(name = paste(unlist(strsplit(isotope, split = '_'))[2],unlist(strsplit(isotope, split = '_'))[3],
                            sep = ''), 
               Z = as.numeric(unlist(strsplit(isotope, split = '_'))[1]),
               A = as.numeric(unlist(strsplit(isotope, split = '_'))[3]),
               sym = unlist(strsplit(isotope, split = '_'))[2],
               t12d = t12, 
               halflife = paste(hfl[3],hfl[4])) 





which(regexpr('Decay Mode', readLines(con = paste('isotopes/',isotope, sep = ''))) == 1)








iso_params <- list(Z = "numeric", 
                   A = "numeric",
                   numT12 = "numeric",
                   unitT12 = "character",
                   id = 'character')




paste(unlist(strsplit(isotope, split = '_'))[2],unlist(strsplit(isotope, split = '_'))[3],
      sep = '')



Isotopes

isotope <- '064_Gd_153'
isotope <- '071_Lu_177'
isotope <- '084_Po_220'

download.file(paste('https://t2.lanl.gov/nis/data/endf/decayVII.1/',isotope, sep = ''), 
              destfile = paste('isotopes/',isotope, sep = ''))




which(regexpr('Parent half-life', readLines(con = paste('isotopes/',isotope, sep = ''))) == 1)
which(pos == 1)



t12 <- strsplit(readLines(con = paste('isotopes/',isotope, sep = '')), split = ' ')
t12n <- as.numeric(t12[3])



Isotopes <- list()
#Isotopes[[1]] <-list(id = 'Ac225', Z = 89)
#names(Isotopes)[1] <- 'Ac225'

isotope <- '064_Gd_153'
isotope <- '089_Ac_225'
isotope <- '083_Bi_213'
isotope <- '097_Bk_249'

download.file('http://www.data-explorer.com/content/data/periodic-table-of-elements-csv.zip', 
              destfile = 'isotopes/periodicTable.zip')


http://www.oecd-nea.org/dbdata/jeff/jeff33/
  
  download.file(paste('https://t2.lanl.gov/nis/data/endf/decayVII.1/',isotope, sep = ''), 
                destfile = paste('isotopes/',isotope, sep = ''))

hfl <- strsplit(readLines(con = paste('isotopes/',isotope, sep = '')), split = ' ')[[
  which(regexpr('Parent half-life', readLines(con = paste('isotopes/',isotope, sep = ''))) == 1)
  ]]

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


# make the decay info

Alpha <- list()
Beta <- list()
Gamma <- list()
Xray <- list()
CE <- list()


dk <- strsplit(readLines(con = paste('isotopes/',isotope, sep = '')), split = ' ')[[
  which(regexpr('Decay Mode', readLines(con = paste('isotopes/',isotope, sep = ''))) == 1)
  ]]
dk <- dk[3:(min(which(dk == ""))-1)]
dk <- unlist(strsplit(dk, split = ','))

ndk <- length(dk)
i = 1
if (dk[i] == 'A') {
  linenum <- which(regexpr('Mean Alpha', readLines(con = paste('isotopes/',isotope, sep = ''))) == 1)
  line <- unlist(strsplit(readLines(con = paste('isotopes/',isotope, sep = ''))[linenum], split = ' '))
  Alpha$Mean <- as.numeric(line[which(line != "")][4])
  Alpha$SD <- as.numeric(line[which(line != "")][6])
} else if (dk[i] == 'B-') { 
  linenum <- which(regexpr('Mean B-', readLines(con = paste('isotopes/',isotope, sep = ''))) == 1)
  line <- unlist(strsplit(readLines(con = paste('isotopes/',isotope, sep = ''))[linenum], split = ' '))
  Beta$Mean <- as.numeric(line[which(line != "")][4])
  Beta$SD <- as.numeric(line[which(line != "")][6])
} 
linenum <- which(regexpr('Mean Gamma', readLines(con = paste('isotopes/',isotope, sep = ''))) == 1)
line <- unlist(strsplit(readLines(con = paste('isotopes/',isotope, sep = ''))[linenum], split = ' '))
Gamma$Mean <- as.numeric(line[which(line != "")][4])
Gamma$SD <- as.numeric(line[which(line != "")][6])

linenum <- which(regexpr('Mean X-Ray', readLines(con = paste('isotopes/',isotope, sep = ''))) == 1)
line <- unlist(strsplit(readLines(con = paste('isotopes/',isotope, sep = ''))[linenum], split = ' '))
Xray$Mean <- as.numeric(line[which(line != "")][4])
Xray$SD <- as.numeric(line[which(line != "")][6])

linenum <- which(regexpr('Mean CE', readLines(con = paste('isotopes/',isotope, sep = ''))) == 1)
line <- unlist(strsplit(readLines(con = paste('isotopes/',isotope, sep = ''))[linenum], split = ' '))
CE$Mean <- as.numeric(line[which(line != "")][4])
CE$SD <- as.numeric(line[which(line != "")][6])

