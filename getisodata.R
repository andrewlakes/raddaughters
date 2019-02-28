library(stringr)
library(PeriodicTable)
library(plyr)
library(ggplot2)
library(plotly)
library(Scale)
library(DiagrammeR)


data("periodicTable")


# download.file('http://www.data-explorer.com/content/data/periodic-table-of-elements-csv.zip', 
#               destfile = 'isotopes/periodicTable.zip')

#### URL for source http://www.oecd-nea.org/dbdata/jeff/jeff33/


Isotopes <- list()

isofile <- 'isotopes/JEFF33-rdd_all.asc' # this is the master isotope data file

# custom function for splitting lines at whitespaces and unlisting
linesplit <- function(x) unlist(strsplit(x, split = " "))[
  which(unlist(strsplit(x, split = " ")) != "")]

iso <- '227AC' # input the parent isotope! 
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

# add the specific activity in mCi / ug
addIso$SA <- (log(2)/(addIso$t12*24*60)) * (6.02214076E23 / addIso$A / 1E6 / 2.22E9)

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
      
      # add the specific activity in mCi / ug
      addIso$SA <- (log(2)/(addIso$t12*24*60)) * (6.02214076E23 / addIso$A / 1E6 / 2.22E9)
      
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


Isotopes$`213BI`$Decays$Beta$daughter



#Create graphic diagram of decay chain

#Show in box: isotope, decay mode, half-life, 
#Show near decay line: energy, and percentage



#make a list practice:

practice = list(c('color', "things", "n stuff"), list(1,2,3), matrix(1:6, nrow=2, ncol=3))
names(practice) = c("1st thing", "second thing", "third thing")
names(practice[[2]]) = c("1", "2", "3")



#####DiagrammeR practice#####


B = "
digraph Isotopes {


graph [overlap = true, fontsize = 10]


node [shape = box, style = filled, penwidth = 2.0, color = 'black', alpha = 50,
fillcolor = '#DDFFEB', fontname = Helvetica]
A; B; C; D; E; F

edge[color=black]
A->B B->C C->D D->E E->F

}

"
grViz(B)

A = "
digraph Isotopes {
      
      
      graph [overlap = true, fontsize = 10]
      
      
      node [shape = box, style = filled, penwidth = 2.0, color = 'black', alpha = 50,
            fillcolor = '#DDFFEB', fontname = Helvetica]
      A [label = '@@1']; B; C; D; E; F
      
      node [shape = circle, fixedsize = true, width = 0.9, penwidth = 2.0,]
      1; 2; 3; 4; 5; 6; 7; 8
      
      edge[color=black]
      A->1 B->2 B->3 B->4 C->A
      1->D E->A 2->4 1->5 1->F
      E->6 4->6 5->7 6->7 3->8
      C->B
      }
      
      [1]: 'top'
      "
grViz(A)


#Using Ac227 Isotopes list:

#These can be combined later
#sequence is
#   preamble
#   nodes
#   amble
#   edges



#                                                   Start

#annoying things:
apostrophe = "'"
openbracket = "["
closebracket = "]"
colon = ":"

#preamble:

preamble = "digraph { graph [overlap = true, fontsize = 10]
                      node [shape = box, style = filled, penwidth = 2.0, color = 'black', alpha = 50,
                            fillcolor = '#DDFFEB', fontname = Helvetica]
                            "

#insert nodes -> e.g. isotopes
#pull from Isotopes

nodes = NA

for (n in 1:length(Isotopes)){
  nodes[n] = rbind(Isotopes[[n]]$isotope)
}

#add apostrophe's to nodes
nodes = paste(apostrophe,nodes,apostrophe)

#add in bracket numbers for calling [#]
brackets = NA
for (n in 1:length(Isotopes)){
  brackets[n] = paste(openbracket,rbind(n),closebracket,colon)
  
}


#add in ";" after each isotope
semicolon = ";"
nodes = paste(nodes,semicolon)
nodes = paste(nodes,collapse="")

#amble

amble = "edge[color=black]
          "
edges = "227AC->223FR"

curlyclose = "}"
  
IsotopeDiagram = paste(preamble,nodes,amble,edges,curlyclose)
grViz(IsotopeDiagram)



"digraph { graph [overlap = true, fontsize = 10]\n      \n                      node [shape = box, style = filled, penwidth = 2.0, color = 'black', alpha = 50,\n                            fillcolor = '#DDFFEB', fontname = Helvetica]\n                             227AC ;223FR ;227TH ;223RA ;223RA ;219RN ;219RN ;215PO ;211PB ;211BI ;207TL ;211PO ; edge[color=black]\n           227AC->223FR }"
