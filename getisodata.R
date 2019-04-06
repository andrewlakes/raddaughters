library(stringr)
library(PeriodicTable)
library(plyr)
library(ggplot2)
#library(plotly)
library(Scale)
library(DiagrammeR)
library(visNetwork)


data("periodicTable")


# download.file('http://www.data-explorer.com/content/data/periodic-table-of-elements-csv.zip', 
#               destfile = 'isotopes/periodicTable.zip')

#### URL for source http://www.oecd-nea.org/dbdata/jeff/jeff33/

setwd('./')

Isotopes <- list()


#Not sure why there are two projects, raddaughters, and Rad_daughters?
#isofile <- '~/Projects/Rad_daughters/isotopes/JEFF33-rdd_all.asc' # this is the master isotope data file
isofile <- '~/raddaughters/JEFF33-rdd_all.asc'

# custom function for splitting lines at whitespaces and unlisting
linesplit <- function(x) unlist(strsplit(x, split = " "))[
  which(unlist(strsplit(x, split = " ")) != "")]

iso <- '225AC' # input the parent isotope! 
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
hfl <- linesplit(grep(iso, readLines(con = isofile), value = TRUE)[1])

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

# add the specific activity in mCi / ug
addIso$SA <- (log(2)/(addIso$t12*24*60)) * (6.02214076E23 / addIso$A / 1E6 / 2.22E9)

# check for which lines match the decay information for the isotope of interest
linematches <- grep(iso, readLines(con = isofile), value = FALSE)

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
  outp <- linesplit(readLines(con = isofile)[linematches[i+1]])
  outq <- linesplit(readLines(con = isofile)[linematches[i+1]+1])
  
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
    IT$daughter <- paste(addIso$A, addIso$symb, sep = '')
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

dk_levs <- 100 # choose the number of levels of the decay chain to follow

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
      hfl <- linesplit(grep(iso, readLines(con = isofile), value = TRUE)[1])
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
      
      # add the specific activity in mCi / ug
      addIso$SA <- (log(2)/(addIso$t12*24*60)) * (6.02214076E23 / addIso$A / 1E6 / 2.22E9)
      
      linematches <- grep(iso, readLines(con = isofile), value = FALSE)
      # ndk_line <- grep(iso, readLines(con = 'isotopes/JEFF33-rdd_all.asc'), value = FALSE)[1]+1
      ndk <- as.integer(linesplit(readLines(con = isofile)[linematches[1]+1])[5])
      
      Decays <- list()
      Alpha <- list()
      Beta <- list()
      Positron <- list()
      EC <- list()
      IT <- list()
      
      for (i in seq(ndk)) {
        outp <- linesplit(readLines(con = isofile)[linematches[i+1]])
        outq <- linesplit(readLines(con = isofile)[linematches[i+1]+1])
        
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
          IT$daughter <- paste(addIso$A, addIso$symb, sep = '')
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
          #Isotopes[[which(match(names(Isotopes),addIso$isotope) == 1)]]$masterYield <- 
          #  Isotopes[[which(match(names(Isotopes),addIso$isotope) == 1)]]$masterYield + addIso$masterYield
          Isotopes[[length(Isotopes)+1]] <- addIso
          names(Isotopes)[length(Isotopes)] <- paste(iso, dk, sep = '_')
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

# check for duplicate isotope entries, merge the masterYield value, and delete the duplicate 
for (i in seq(length(Isotopes))){
  if (length(unlist(strsplit(names(Isotopes[i]), split = '_'))) == 2) {
    match <- which(names(Isotopes) == unlist(strsplit(names(Isotopes[i]), split = '_'))[1])
    Isotopes[[match]]$masterYield <- Isotopes[[match]]$masterYield + Isotopes[[i]]$masterYield
    Isotopes[i] <- list(NULL)
  }
}

#if there are no dupes, ignore -> had error original way
if (length(which(lapply(Isotopes, length) == 0)) != 0){Isotopes <- Isotopes[-c(which(lapply(Isotopes, length) == 0))]}
#original:
#Isotopes <- Isotopes[-c(which(lapply(Isotopes, length) == 0))]




#saveRDS(Isotopes, paste('decayLists/', Isotopes[[1]]$isotope, sep = ''))




#make a list practice:
#practice = list(c('color', "things", "n stuff"), list(1,2,3), matrix(1:6, nrow=2, ncol=3))
#names(practice) = c("1st thing", "second thing", "third thing")
#names(practice[[2]]) = c("1", "2", "3")




####Instead of DiagrammeR -> using visNetwork####


#Ignore insignificant nodes....
#first to determine ALL nodes


#IsotopesTerm = NA
#for (n in 1:length(Isotopes)){
#1) if decay daughter does not match any Isotopes[#], this is a terminal isotope, and create a new NODE
#    if ((length(which((str_detect(str_detect(names(Isotopes), coll(Isotopes[[n]]$Decays$Alpha$daughter, ignore_case=TRUE)), coll('TRUE'))))) == 0) 
#      & 
#      (length(which((str_detect(str_detect(names(Isotopes), coll(Isotopes[[n]]$Decays$Beta$daughter, ignore_case=TRUE)), coll('TRUE'))))) == 0) 
#      &
#      (length(which((str_detect(str_detect(names(Isotopes), coll(Isotopes[[n]]$Decays$Positron$daughter, ignore_case=TRUE)), coll('TRUE'))))) == 0) 
#      &
#      (length(which((str_detect(str_detect(names(Isotopes), coll(Isotopes[[n]]$Decays$EC$daughter, ignore_case=TRUE)), coll('TRUE'))))) == 0) 
#      &
#      (length(which((str_detect(str_detect(names(Isotopes), coll(Isotopes[[n]]$Decays$IT$daughter, ignore_case=TRUE)), coll('TRUE'))))) == 0)) 
#            {IsotopesTerm[n] = names(Isotopes[n])}
#}

#Find what the final isotope is from those found:
#IsotopesTerminal = NA
#for (n in which(!is.na(IsotopesTerm) )){
#  IsotopesTerminal[n] = Isotopes[[n]][[9]][[1]][[4]]
#}

#clean the final isotope list
#IsotopesTerminal = na.omit(unique(IsotopesTerminal))


#set up names - pull from Isotopes master list
#nodesnames = names(Isotopes)

nodesnames = NA

for (n in 1:(length(Isotopes))){
  nodesnames[n] = rbind(paste(Isotopes[[n]]$A,Isotopes[[n]]$symb))
}

#add in terminal isotopes
#for (n in 1:(length(IsotopesTerminal))){
#  nodesnames[length(Isotopes)+n] = rbind(IsotopesTerminal[n])
#}

#remove spaces
nodesnames = str_replace_all(string=nodesnames, pattern=" ", repl="")

#Node attributes
nodesid = 1:length(nodesnames)
nodeslabel = nodesnames
nodesvalue = 1
nodesshape = 'circle'
nodestitle = paste0("<p><b>", 1:length(nodesnames),"</b><br>Node title goes here</p>")
nodeslevel = edgesfrom
#nodescolor = 'green'

#create master nodes data.frame
nodes = data.frame(id = nodesid, 
                   label = nodeslabel, 
                   value = nodesvalue, 
                   shape = nodesshape, 
                   title = nodestitle, 
                   shadow = TRUE ,
                   level = nodeslevel)







#find edges
#numbers in IsotopesEdgeto correspond to isotopes in nodes list order.
IsotopesEdgeto = matrix(NA, nrow = 5, ncol = length(nodes[,2]))
colnames(IsotopesEdgeto) = nodes[,2]
rownames(IsotopesEdgeto) = c('Alpha', 'Beta', 'Positron', 'EC', 'IT')




for (i in 1:(length(nodes))){
  if ((!is.null(Isotopes[[i]]$Decays$Alpha$daughter))&(isTRUE(Isotopes[[i]]$Decays$Alpha$daughter %in% nodes[,2]))) {IsotopesEdgeto[1,i] = which((str_detect(str_detect(names(Isotopes), coll(Isotopes[[i]]$Decays$Alpha$daughter, ignore_case=TRUE)), coll('TRUE'))))}
  if ((!is.null(Isotopes[[i]]$Decays$Beta$daughter))&(isTRUE(Isotopes[[i]]$Decays$Beta$daughter %in% nodes[,2]))) {IsotopesEdgeto[2,i] = which((str_detect(str_detect(names(Isotopes), coll(Isotopes[[i]]$Decays$Beta$daughter, ignore_case=TRUE)), coll('TRUE'))))}
  if ((!is.null(Isotopes[[i]]$Decays$Positron$daughter))&(isTRUE(Isotopes[[i]]$Decays$Positron$daughter %in% nodes[,2]))) {IsotopesEdgeto[3,i] = which((str_detect(str_detect(names(Isotopes), coll(Isotopes[[i]]$Decays$Beta$daughter, ignore_case=TRUE)), coll('TRUE'))))}
  if ((!is.null(Isotopes[[i]]$Decays$EC$daughter))&(isTRUE(Isotopes[[i]]$Decays$EC$daughter %in% nodes[,2]))) {IsotopesEdgeto[4,i] = which((str_detect(str_detect(names(Isotopes), coll(Isotopes[[i]]$Decays$Beta$daughter, ignore_case=TRUE)), coll('TRUE'))))}
  if ((!is.null(Isotopes[[i]]$Decays$IT$daughter))&(isTRUE(Isotopes[[i]]$Decays$IT$daughter %in% nodes[,2]))) {IsotopesEdgeto[5,i] = which((str_detect(str_detect(names(Isotopes), coll(Isotopes[[i]]$Decays$Beta$daughter, ignore_case=TRUE)), coll('TRUE'))))}
  }
  


#now split into groups based on decay mode
#can create third dimension if you want on probability later#

#create the vectors that will be inserted into the final visNet
#numbers still correlate with isotope order of 'nodes'

#edgesto = matrix(NA, nrow = length(which(!is.na(IsotopesEdgeto))), ncol = 1)
#edgesfrom = matrix(NA, nrow = length(which(!is.na(IsotopesEdgeto))), ncol = 1)
edgesto = NA
edgesfrom = NA




for (j in 1:length(IsotopesEdgeto[1,])) {
  for (i in 1:length(IsotopesEdgeto[,1])) {
    if (!is.na(IsotopesEdgeto[i,j])) {edgesfrom = c(edgesfrom, j)}
    if (!is.na(IsotopesEdgeto[i,j])) {edgesto = c(edgesto, IsotopesEdgeto[i,j])}
    
  }
}

#remove first NA from rows

edgesto = edgesto[-1]
edgesfrom = edgesfrom[-1]






#test values -> move after find edges later
#Edges variables
edgeslabel = 'hi'
edgeslength = 1
edgeswidth = 1
edgesarrows = "to"
edgesdashes = TRUE
edgestitle = paste("Edge Name HERERE", 1:length(edgesfrom))
edgessmooth = FALSE
edgesshadow = TRUE

edges = data.frame(from = edgesfrom, 
                   to = edgesto, 
                   label = edgeslabel, 
                   length = edgeslength, 
                   width = edgeswidth,
                   arrows = edgesarrows, 
                   dashes = edgesdashes, 
                   title = edgestitle, 
                   smooth = edgessmooth,
                   shadow = edgesshadow)



#test

visNetwork(nodes, edges, width = "100%")%>%
  visHierarchicalLayout()


