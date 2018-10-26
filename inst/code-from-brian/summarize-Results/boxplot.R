
 pi <- c( H$PP1[,1], H$PP2[,1],
          H$PP1[,2], H$PP2[,2],
          H$PP1[,3], H$PP2[,3],
          H$PP1[,4], H$PP2[,4],
          H$PP1[,5], H$PP2[,5],
          H$PP1[,6], H$PP2[,6] )

 x <- c(  rep("1 NSCLC", length(pi)/12 ), rep("1 NSCLC", length(pi)/12 ),
          rep("2 CRC.v", length(pi)/12 ), rep("2 CRC.v", length(pi)/12 ), 
          rep("3 CRC.vc", length(pi)/12 ), rep("3 CRC.vc", length(pi)/12 ), 
          rep("4 BD", length(pi)/12 ), rep("4 BD", length(pi)/12 ), 
          rep("5 ED.LH", length(pi)/12 ), rep("5 ED.LH", length(pi)/12 ), 
          rep("6 ATC", length(pi)/12 ), rep("6 ATC", length(pi)/12 ) )

 lab <- c( rep("MEM", length(pi)/12 ), rep("SEM", length(pi)/12 ),
          rep("MEM", length(pi)/12 ), rep("SEM", length(pi)/12 ), 
          rep("MEM", length(pi)/12 ), rep("SEM", length(pi)/12 ), 
          rep("MEM", length(pi)/12 ), rep("SEM", length(pi)/12 ), 
          rep("MEM", length(pi)/12 ), rep("SEM", length(pi)/12 ), 
          rep("MEM", length(pi)/12 ), rep("SEM", length(pi)/12 ) )

 PP <- as.data.frame( list(pi=pi, x=x, Method=factor(lab)))

# p <- ggplot(PP, aes(x=factor(x), y=pi, fill=Method)) + 
#      geom_boxplot(size=1, position=position_dodge(0.82)) + theme_classic() +
#      labs(title="",x="Basket", y = "Posterior Probability") +
#      theme( axis.title.x = element_text(size=19), axis.title.y = element_text(size=19), text = element_text(size=18),
#             legend.text = element_text(size = 25),
#             legend.key.size = unit(1, "cm"),
#             legend.title=element_text(face = "bold", size = 10),
#             axis.text.x = element_text(colour="black", size = 10) )
#
# p

 p <- ggplot(PP, aes(x=factor(x), y=pi, fill=Method)) + 
      geom_boxplot(size=1, position=position_dodge(0.82), show.legend=F) + theme_classic() +
      labs(title="",x="Basket", y = "Posterior Probability") +
      theme( axis.title.x = element_text(size=19, margin=margin(t = 18)), 
             axis.title.y = element_text(size=19, margin=margin(r = 15)), 
             text = element_text(size=20),
             legend.text = element_text(size = 25),
             legend.key.size = unit(1, "cm"),
             plot.margin = unit(c(1,0.9,1,1), "cm"),
             legend.title=element_text(face = "bold", size = 10),
             axis.text.x = element_text(colour="black", size = 10) ) + 
       scale_y_continuous(limits = YLIM)


