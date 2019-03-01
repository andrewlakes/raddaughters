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

setwd('./')

Isotopes <- list()

isofile <- '~/Projects/Rad_daughters/isotopes/JEFF33-rdd_all.asc' # this is the master isotope data file

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

#saveRDS(Isotopes, './decayLists/227AC')

# check for duplicate isotope entries, merge the masterYield value, and delete the duplicate 
for (i in seq(length(Isotopes))){
  if (length(unlist(strsplit(names(Isotopes[i]), split = '_'))) == 2) {
    match <- which(names(Isotopes) == unlist(strsplit(names(Isotopes[i]), split = '_'))[1])
    Isotopes[[match]]$masterYield <- Isotopes[[match]]$masterYield + Isotopes[[i]]$masterYield
    Isotopes[i] <- list(NULL)
  }
}
Isotopes <- Isotopes[-c(which(lapply(Isotopes, length) == 0))]



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
1[label='@@1'] 2[label='@@2'] 3[label='@@3'] 4[label='@@4'] 5[label='@@5'] 6[label='@@6'] 7[label='@@7'] 8[label='@@8'] 9[label='@@9'] 10[label='@@10'] 11[label='@@11'] 12[label='@@12']

edge[color=black]
1->2 2->3 3->4 4->5 5->6

}
[1]:'ac'
[2]:'th'
[3]:'th'
[4]:'th'
[5]:'th'
[6]:'6'
[7]:'7'
[8]:'8'
[9]:'9'
[10]:'10'
[11]:'11'
[12]:'12'

"
grViz(B)

B = "digraph { graph [overlap = true, fontsize = 10]                      
node [shape = box, style = filled, penwidth = 2.0, color = 'black', alpha = 50, fillcolor = '#DDFFEB', fontname = Helvetica]
1[label='@@1'] 2[label='@@2'] 3[label='@@3'] 4[label='@@4'] 5[label='@@5'] 6[label='@@6'] 7[label='@@7'] 8[label='@@8'] 9[label='@@9']
edge[color=black] 1->2 
} 
[1]: '227AC' \n [2]: '223FR' \n [3]: '227TH' \n [4]: '223RA' \n [5]: '223RA' \n [6]: '219RN' \n [7]: '219RN' \n [8]: '215PO' \n [9]: '211PB' "










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





#sequence is
#   preamble
#   TerminalIsotopes
#   nodes
#   amble
#   edges
#   nodesnames

#  @@# calls at the bottom after curly bracket corresponding number
#       for subsets of same number, use e.g. @@1-1, @@1-2, @@2, @@3
#       corresponds to [1] \n [2] \n [3]




#                                                   Start

#word soup:
apostrophe = "'"
openbracket = "["
closebracket = "]"
colon = ":"
semicolon = ";"
atat = "@@"
curlyclose = "}"
slashn = "\n"
arrow = "->"

#preamble:
preamble = "digraph { graph [overlap = true, fontsize = 10]
                      node [shape = box, style = filled, penwidth = 2.0, color = 'black', alpha = 50,
                            fillcolor = '#DDFFEB', fontname = Helvetica]
                            "


#edges first to determine ALL nodes
IsotopesTerm = NA
for (n in 1:length(Isotopes)){
#1) if decay daughter does not match any Isotopes[#], this is a terminal isotope, and create a new NODE
    if ((length(which((str_detect(str_detect(names(Isotopes), coll(Isotopes[[n]]$Decays$Alpha$daughter, ignore_case=TRUE)), coll('TRUE'))))) == 0) 
      & 
      (length(which((str_detect(str_detect(names(Isotopes), coll(Isotopes[[n]]$Decays$Beta$daughter, ignore_case=TRUE)), coll('TRUE'))))) == 0) 
      &
      (length(which((str_detect(str_detect(names(Isotopes), coll(Isotopes[[n]]$Decays$Positron$daughter, ignore_case=TRUE)), coll('TRUE'))))) == 0) 
      &
      (length(which((str_detect(str_detect(names(Isotopes), coll(Isotopes[[n]]$Decays$EC$daughter, ignore_case=TRUE)), coll('TRUE'))))) == 0) 
      &
      (length(which((str_detect(str_detect(names(Isotopes), coll(Isotopes[[n]]$Decays$IT$daughter, ignore_case=TRUE)), coll('TRUE'))))) == 0)) 
            {IsotopesTerm[n] = names(Isotopes[n])}
}

#Find what the final isotope is from those found:
IsotopesTerminal = NA
for (n in which(!is.na(IsotopesTerm) )){
  IsotopesTerminal[n] = Isotopes[[n]][[9]][[1]][[4]]
}

#clean the final isotope list
IsotopesTerminal = na.omit(unique(IsotopesTerminal))




#insert nodes -> e.g. isotopes
nodes = NA

for (n in 1:(length(Isotopes)+length(IsotopesTerminal))){
  nodes[n] = rbind(paste(atat,n))
}


#add letters for nodes and make @@ as a label
#syntax
nodewords = "[label = '"

for (n in 1:(length(Isotopes)+length(IsotopesTerminal))){
  nodes[n] = rbind(paste(n,nodewords,nodes[n],apostrophe,closebracket))
}

#remove spaces
nodes = str_replace_all(string=nodes, pattern=" ", repl="")

#collapse vectors to single strings
nodes = paste(nodes,collapse=" ")



#set up names - pull from Isotopes master list
#nodesnames = names(Isotopes)

nodesnames = NA

for (n in 1:(length(Isotopes))){
  nodesnames[n] = rbind(paste(Isotopes[[n]]$A,Isotopes[[n]]$symb))
}

#add in terminal isotopes
for (n in 1:(length(IsotopesTerminal))){
  nodesnames[length(Isotopes)+n] = rbind(IsotopesTerminal[n])
}

#add apostrophe's and \n to nodes
nodesnames = paste(apostrophe,nodesnames,apostrophe,slashn)

#remove spaces
nodesnames = str_replace_all(string=nodesnames, pattern=" ", repl="")


#add in bracket numbers for calling [#]
brackets = NA
for (n in 1:(length(Isotopes)+length(IsotopesTerminal))){
  brackets[n] = paste(openbracket,rbind(n),closebracket,colon)
  
}

#remove spaces
brackets = str_replace_all(string=brackets, pattern=" ", repl="")

#combine brackets and nodes names
nodesnames = paste(brackets,nodesnames)

#collapse vectors to single strings
nodesnames = paste(nodesnames,collapse=" ")

#amble
amble = "edge[color=black]"

#edges
edges = "1->2 2->8 4->12"







#2) if alpha is present continue
#3) find $daughter string, and find which Isotopes[#] it is 
#4) output #-># into a vector
#5) if beta is present, continue
#6) find $daughter string, and find which Isotopes[#] it is 
#7) output #-># into a vector
#8) win



IsotopeDiagram = paste(preamble,nodes,amble,edges,slashn,curlyclose,slashn,nodesnames)
grViz(IsotopeDiagram)
