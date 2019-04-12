library(stringr)
library(PeriodicTable)
library(plyr)
library(ggplot2)
#library(plotly)
library(Scale)
library(DiagrammeR)
library(visNetwork)
library(tidyr)


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

iso <- '227AC' # input the parent isotope! 

Isotopes <- readRDS(paste('decayLists/', iso, sep = ''))







#set up names - pull from Isotopes master list
nodesnames = names(Isotopes)

nodesnames = NA

for (n in 1:(length(Isotopes))){
  nodesnames[n] = rbind(paste(Isotopes[[n]]$A,Isotopes[[n]]$symb))
}

# #add in terminal isotopes
# for (n in 1:(length(IsotopesTerminal))){
#  nodesnames[length(Isotopes)+n] = rbind(IsotopesTerminal[n])
# }

#remove spaces
nodesnames = str_replace_all(string=nodesnames, pattern=" ", repl="")

nodesdata = NA

for (i in 1:length(Isotopes)) {
  nodesdata[i] = paste("<p><b>",Isotopes[[i]]$isotope,"</b><p>", "Z = ",Isotopes[[i]]$Z, "</b><p>", "Half-life (Days) = ", Isotopes[[i]]$t12, "</b><p>", "SA (Ci/mg) = ", Isotopes[[i]]$SA)
}



#Node attributes
nodesid = 1:length(nodesnames)
nodeslabel = nodesnames
nodesvalue = 1
nodesshape = 'circle'
nodestitle = nodesdata
#temp insert, 'edgesfrom' is reevaluated later
nodeslevel = 1:length(nodesnames)#edgesfrom
#nodescolor = 'green'

#create master nodes data.frame
nodes = data.frame(id = nodesid, 
                   label = nodeslabel, 
                   value = nodesvalue, 
                   shape = nodesshape, 
                   title = nodestitle, 
                   shadow = TRUE,
                   level = nodeslevel)







#find edges
#numbers in IsotopesEdgeto correspond to isotopes in nodes list order.
IsotopesEdgeto = matrix(NA, nrow = 5, ncol = length(nodes[,2]))
colnames(IsotopesEdgeto) = nodes[,2]
rownames(IsotopesEdgeto) = c('Alpha', 'Beta', 'Positron', 'EC', 'IT')


#do a name check, and/or for terminal isotope name in 'IsotopesTerminal'


for (i in 1:(length(nodes[,2]))){
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



#edgesfrom is now complete, redo 'nodes level' based on isotope Z in 'nodes$title

nodeslevels = NA

for (i in 1:length(nodes[,1])){
  nodes$level[i] = -Isotopes[[i]]$Z
  
}





#Vairable width based on probability

#for each to-from pair, uses order found in 'nodes' which matches isotope order. Pull directly from 'Isotopes'
#ncol is number of types of destructions, alpha beta positron EC IT.
edgesthickness = matrix(NA, nrow = length(nodes[,2]), ncol = 5)

for (i in 1:length(edgesthickness[,1])) {
  if ((!is.null(Isotopes[[i]]$Decays$Alpha$branchYiel))&(isTRUE(Isotopes[[i]]$Decays$Alpha$daughter %in% nodes[,2]))) {edgesthickness[i,1] = Isotopes[[i]]$Decays$Alpha$branchYield}
  if ((!is.null(Isotopes[[i]]$Decays$Beta$branchYiel))&(isTRUE(Isotopes[[i]]$Decays$Beta$daughter %in% nodes[,2]))) {edgesthickness[i,2] = Isotopes[[i]]$Decays$Beta$branchYield}
  if ((!is.null(Isotopes[[i]]$Decays$Positron$branchYiel))&(isTRUE(Isotopes[[i]]$Decays$Positron$daughter %in% nodes[,2]))) {edgesthickness[i,3] = Isotopes[[i]]$Decays$Positron$branchYield}
  if ((!is.null(Isotopes[[i]]$Decays$EC$branchYiel))&(isTRUE(Isotopes[[i]]$Decays$EC$daughter %in% nodes[,2]))) {edgesthickness[i,4] = Isotopes[[i]]$Decays$EC$branchYield}
  if ((!is.null(Isotopes[[i]]$Decays$IT$branchYiel))&(isTRUE(Isotopes[[i]]$Decays$IT$daughter %in% nodes[,2]))) {edgesthickness[i,5] = Isotopes[[i]]$Decays$IT$branchYield}
  
}

rownames(edgesthickness) = colnames(IsotopesEdgeto)
colnames(edgesthickness) = c('Alpha', 'Beta', 'Positron', 'EC', 'IT')
#Transpose
edgesthickness = t(edgesthickness)
#now edgesthickness should match 'IsotopesEdgeto'

#collapse list into order of 'nodes[,2]'

edgesthicknessorder = matrix(NA, nrow = length(nodes[,2]), ncol = 1)

for (j in 1:length(edgesthickness[1,])) {
  for (i in 1:length(edgesthickness[,1])) {
    if (!is.na(edgesthickness[i,j])) {edgesthicknessorder = c(edgesthicknessorder, edgesthickness[i,j])}
  }
}

#collapse NA's
edgesthicknessorder = na.omit(edgesthicknessorder)
#make thicker for visual effect
edgesthicknessorder = edgesthicknessorder*5

#add colors to each type of destruction
#color master
colorkey = c("blue", "green", "orange", "purple", "red")
colorkey = t(colorkey)
colnames(colorkey) = rownames(edgesthickness)


#describe each destruction type
edgesthicknesscolor = matrix(NA, nrow = length(nodes[,2]), ncol = 1)

for (j in 1:length(edgesthickness[1,])) {
  for (i in 1:length(edgesthickness[,1])) {
    if (!is.na(edgesthickness[i,j])) {edgesthicknesscolor = c(edgesthicknesscolor, rownames(edgesthickness)[i])}
  }
}

#collapse NA's
edgesthicknesscolor = na.omit(edgesthicknesscolor)

#correlate name with 'colorkey'

edgesthicknesscolorname = NA
for (i in 1:length(edgesthicknesscolor)){
  edgesthicknesscolorname[i] = colorkey[match(edgesthicknesscolor[i], colnames(colorkey))]
}


#combine all into one
edgesthicknesscombo = cbind(edgesthicknesscolor, edgesthicknesscolorname, edgesthicknessorder)






#test values -> move after find edges later
#Edges variables
edgeslabel = 'hi'
edgeslength = 1
edgeswidth = edgesthicknesscombo[,3]
edgesarrows = "to"
edgesdashes = FALSE
edgestitle = paste("Edge Name HERERE", 1:length(edgesfrom))
edgessmooth = FALSE
edgesshadow = TRUE
edgescolor = edgesthicknesscombo[,2]

edges = data.frame(from = edgesfrom, 
                   to = edgesto, 
                   label = edgeslabel, 
                   length = edgeslength, 
                   width = edgeswidth,
                   arrows = edgesarrows, 
                   dashes = edgesdashes, 
                   title = edgestitle, 
                   smooth = edgessmooth,
                   shadow = edgesshadow,
                   color = edgescolor)



#test

visNetwork(nodes, edges, width = "100%", height = "100%")%>%
  visHierarchicalLayout()%>%
  visNodes(color = list(highlight = "white"))%>%
  visHierarchicalLayout(levelSeparation = 75)%>%
  visInteraction(dragNodes = TRUE,
                 hideEdgesOnDrag = TRUE,
                 dragView = TRUE, 
                 zoomView = TRUE,
                 navigationButtons = TRUE,
                 tooltipDelay = 50,
                 keyboard = TRUE)


