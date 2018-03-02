
PSI_geneEvents <- function(events, type, showReplicates = TRUE){
  
  as_types <- c("A3SS", "A5SS", "SE", "RI", "MXE")
  if (!type %in% as_types){
    stop(cat("\"type\" should be one of the following: ", as_types))
  }
  
  annot <- events[[paste0(type,"_","events")]]
  if (length(unique(annot$geneSymbol)) > 1){
    stop(cat("Multiple genes found. Use geneEvents() to select AS events."))
  }
  
  PSI <- events[[paste0(type,"_","PSI")]]
  PSI_long <- reshape2::melt(PSI)
  colnames(PSI_long) <- c("ID", "Sample", "PSI")
  
  Condition <- rep("NA",nrow(PSI_long))
  idx.cond1 <- grep(paste0("^", events$conditions[1]), x = PSI_long$Sample,
                    perl = T )
  idx.cond2 <- grep(paste0("^", events$conditions[2]), x = PSI_long$Sample,
                    perl = T )
  Condition[idx.cond1] <- events$conditions[1]
  Condition[idx.cond2] <- events$conditions[2]
  
  PSI_long <- cbind(PSI_long, Condition)
  
  if (showReplicates){
    
    ggplot(PSI_long, aes(x = Condition, y = PSI, fill = Condition, color = Condition)) +
      #geom_boxplot() +
      geom_violin(trim = F, alpha = 0.6) +
      geom_jitter(position=position_jitter(0.05), size = 2) +
      theme_bw() +
      theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1),
            axis.text.y = element_text(size=12),
            axis.title.x = element_text(face="plain", colour="black", size=12),
            axis.title.y = element_text(face="plain", colour="black", size=12),
            legend.title=element_blank(),
            legend.text = element_text(face="plain", colour="black", size=12)) +
      ylab(paste(type, "PSI")) +
      #xlab("Sample") +
      scale_y_continuous(limits=c(-0.1, 1.05)) +
      scale_fill_manual(values = c("blue", "red") ) +
      scale_color_manual(values = c("blue", "red") ) +
      facet_grid(. ~ ID)

  } else{  
  
  ggplot(PSI_long, aes(x = Condition, y = PSI, fill = Condition, color = Condition)) +
    #geom_boxplot() +
    geom_violin(trim = F) +
    stat_summary(fun.y=mean, geom="point", size=2, color="black") +
    theme_bw() +
    theme(axis.text.x = element_text(size=12, angle = 45, hjust = 1),
          axis.text.y = element_text(size=12),
          axis.title.x = element_text(face="plain", colour="black", size=12),
          axis.title.y = element_text(face="plain", colour="black", size=12),
          legend.title=element_blank(),
          legend.text = element_text(face="plain", colour="black", size=12)) +
    ylab(paste(type, "PSI")) +
    #xlab("Sample") +
    scale_y_continuous(limits=c(-0.1, 1.05)) +
    scale_fill_manual(values = c("blue", "red") ) +
    scale_color_manual(values = c("blue", "red") ) +  
    facet_grid(. ~ ID)
    
  }
  
}