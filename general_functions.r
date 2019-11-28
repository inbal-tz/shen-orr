find_many2many_groups = function(data, min_datapoints){   #fits panther
  require(readr)
  require(dplyr)
  require(plyr)
  result_df = data.frame()
  group_num=1
  many2many <- data %>% filter(ortholog_type == "n:m")
  mouse_geneList = unique(many2many$mouse_entrez_ID)
  while (length(mouse_geneList)>0){
    human_orthologs_group = c()
    mouse_orthologs_group = c(mouse_geneList[1])
    to_add_mouse = c()
    to_add_human = setdiff(unique((data %>% filter(mouse_entrez_ID %in% mouse_orthologs_group))$human_entrez_ID),human_orthologs_group)
    while (length(to_add_mouse)>0 | length(to_add_human)>0){
      human_orthologs_group = unique(c(human_orthologs_group,to_add_human))
      mouse_orthologs_group = unique(c(mouse_orthologs_group,to_add_mouse))
      to_add_mouse = setdiff(unique((data %>% filter(human_entrez_ID %in% human_orthologs_group))$mouse_entrez_ID),mouse_orthologs_group)
      to_add_human = setdiff(unique((data %>% filter(mouse_entrez_ID %in% mouse_orthologs_group))$human_entrez_ID),human_orthologs_group)
    }
    mouse_geneList=mouse_geneList[! mouse_geneList %in% mouse_orthologs_group]
    orth_with_sufficient_data = setNames(data.frame(matrix(ncol = length(human_orthologs_group), nrow = length(mouse_orthologs_group)),row.names = mouse_orthologs_group), human_orthologs_group)
    for (i in mouse_orthologs_group){
      for (j in human_orthologs_group){
        orth_with_sufficient_data[toString(i),toString(j)]= length(unique((data %>% filter(mouse_entrez_ID == i & human_entrez_ID == j))$comp_file))
      }
    }
    orth_with_sufficient_data = orth_with_sufficient_data[!apply(orth_with_sufficient_data < min_datapoints, 1, all), !apply(orth_with_sufficient_data < min_datapoints, 2, all), drop = FALSE]
    
    temp <- data %>% filter(mouse_entrez_ID %in% rownames(orth_with_sufficient_data) & human_entrez_ID %in% colnames(orth_with_sufficient_data)) %>% select(c(human_entrez_ID,mouse_entrez_ID)) %>% distinct()
    if(nrow(temp)>0) {
      temp['group']=group_num
      result_df = rbind(result_df,temp)
      group_num=group_num+1
    }
  }
  return(result_df) 
}

disease_to_shape <- function(diseaseList){
  diseases = c("Adenocarcinoma","Cancer","EcoliInfection","EpilepsyAutism","Huntington","LeptospiralInfection","Lupus","MALTlymphoma","OvarianCancer","StaphInfection","Alpha1Anti","BrainTumor","Burns","CysticFibrosis","DMD","DownSyndrome","Gaucher","IBD","Injury","MyotonicDystrophy","RA","AtopicDermatitis","ProstateCancer","ASPS")
  shapes = c(15,16,8,15,16,8,15,16,8,15,16,8,15,16,8,15,16,8,15,16,8,15,16,8)
  return(shapes[match(diseaseList,diseases)])
}

disease_to_color <- function(diseaseList){
  diseases=c("Adenocarcinoma","Cancer","EcoliInfection","EpilepsyAutism","Huntington","LeptospiralInfection","Lupus","MALTlymphoma","OvarianCancer","StaphInfection","Alpha1Anti","BrainTumor","Burns","CysticFibrosis","DMD","DownSyndrome","Gaucher","IBD","Injury","MyotonicDystrophy","RA","AtopicDermatitis","ProstateCancer","ASPS")
  colors=c("black","black","black","red","red","red","blue","blue","blue","pink","pink","pink","green","green","green","yellow","yellow","yellow","purple","purple","purple","orange","orange","orange")
  return(colors[match(diseaseList,diseases)])
}

multiplot <- function(..., plotlist = NULL, file, cols = 1, layout = NULL, title="", 
                      fontsize = 12, fontfamily = "Helvetica") {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (nchar(title)>0){
    layout <- rbind(rep(0, ncol(layout)), layout)
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), 
                                               ncol(layout), 
                                               heights = if (nchar(title)>0) {unit(c(0.5, rep(5,nrow(layout)-1)), "null")}
                                               else {unit(c(rep(5, nrow(layout))), "null")})))
    
    # Make each plot, in the correct location
    if (nchar(title)>0) {
      grid.text(title, 
                vp = viewport(layout.pos.row = 1, layout.pos.col = 1:ncol(layout)),
                gp = gpar(fontsize = fontsize, fontfamily = fontfamily))
    }
    
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}