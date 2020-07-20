#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

#options(error=NULL)
library(shiny)
library(ggplot2)
library(fst)
library(ldmap)
library(plotly)
library(dplyr)
library(ldmap)

gene_df <- fst::read_fst("data/gene_df.fst") %>%
  as_tibble() %>%
  mutate(ldmr=as.character(as_ldmap_region(ldmr)))

anno_tm_df <- fst::read_fst("data/anno_tm_df.fst") %>%
    mutate(
      ld_snp = as_ldmap_snp(ld_snp),
      snp_pos = as_ldmap_snp(snp_pos),
      gwas_snp=as_ldmap_snp(gwas_snp),
           pos=positions(gwas_snp),
           z=abs(beta)/se,
           p=pnorm(z,lower.tail=F)) %>%
  as_tibble()


# Define server logic required to draw a histogram
s_ldmr <- gene_df$ldmr[1]

subset_gdf <- function(inp){
#  browser()
  tgdf <- dplyr::filter(gene_df,symbol==inp)
  tldmr <-tgdf$ldmr[1]
  wgdf <-dplyr::filter(gene_df,ldmr==tldmr)
  val_x <- c(min(tgdf$start),max(tgdf$end))
  scale_x <- c(min(wgdf$start),max(wgdf$end))
  return(list(tldmr=tldmr,value=c(max(val_x[1]-5000,scale_x[1]),min(val_x[2]+5000,scale_x[2])),scale_x=scale_x))
}


gen2_sin <-function(n){
  sin_a <-sin(seq(from=-2*pi,to = 2*pi,length.out = n/2))
  sin_b <-sin(seq(from=-2*pi,to = 2*pi,length.out = n/2))
  purrr::flatten_dbl(purrr::map2(sin_a,sin_b,~c(.x,.y)))[1:n]
}

gen_plots <-function(sv,s_ldmr,mods,annos){

  rsldmr  <- ldmap::as_ldmap_region(s_ldmr)
  rstarts <- starts(rsldmr)
  rends   <- ends(rsldmr)
  sv[1] <- max(rstarts,sv[1])
  sv[2] <- min(rends,sv[2])
  mods <- unique(c(mods,"null"))
  gdf <- dplyr::filter(gene_df,
                       ldmr==s_ldmr,
                       start>=sv[1] | end <= sv[2]) %>%
    mutate(mid = start+(end-start)/2) %>%
    arrange(mid) %>%
    mutate(flip = rep(c(1.5, -1.5),
                      length.out = n())+gen2_sin(n()))
  xdf <- dplyr::filter(anno_tm_df,
                       ldmr == s_ldmr,
                       model %in% mods,
                       dplyr::between(pos,sv[1],sv[2]),
                       anno %in% annos)
  pva <-ggplot(xdf,aes(x=pos,y=-log10(p),col=anno,label=rsid)) +
    geom_point() +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    coord_cartesian(xlim =sv)+theme(legend.position = "none")
  ppva <- ggplotly(pva)

  pvb <-(ggplot(xdf,aes(x=pos,y=pip,col=anno,label=rsid)) +
       geom_point() +
       theme(axis.title.x = element_blank(),
             axis.text.x = element_blank(),
             axis.ticks.x = element_blank())+
       coord_cartesian(xlim =sv) +
       facet_grid(model~.,labeller = "label_value"))
  ppvb <- ggplotly(pvb)
#  ppvb <- subplot(pvbl,nrows = length(pvbl))

  pvc <- (ggplot(gdf,
                 aes(xmin=start,xmax=end,fill=symbol,label=symbol,y=0.0)
                 ) +
            geom_rect(ymin=-0.05,ymax=.05) + geom_segment(aes(x=mid,xend=mid,y=0,yend=flip,col=symbol))+
            geom_text(aes(x=mid,y=0.0,col=symbol),nudge_y=gdf$flip) +
            coord_cartesian(xlim =sv) +
            theme(axis.title.y=element_blank(),
                  legend.position = "none",
                  plot.margin = margin(0, 0, 0, 0)))+xlab("position")+ylim(-2.5,2.5)
  ppvc <- ggplotly(pvc)
  sp <-subplot(ppva,ppvb,ppvc,nrows=3,shareX=TRUE)
  twgl <- suppressWarnings(toWebGL(sp))
  # ret_plot <- wrap_plots(pva,pvb,pvc,ncol=1,nrow=3,heights = c(1,3,1),guides = 'collect') +
  #   theme(plot.margin = margin(0, 0, 0, 0)) +
  #   plot_annotation(s_ldmr)
  # return(ret_plot)
  return(twgl)
}


shinyServer(function(input, output,session) {



  observe({

    # xdf <- dplyr::filter(anno_tm_df,
    #                      ldmr == s_ldmr,
    #                      model %in% mods,
    #                      anno %in% annos)


    gene_inp <-input$gene
    sub_res <-subset_gdf(gene_inp)
    updateSelectInput(session,"regions",selected=sub_res$tldmr)
    updateSliderInput(session, "receive", value = sub_res$value,
                      min = sub_res$scale_x[1], max = sub_res$scale_x[2])
  })

  observe({
    s_ldmr <-input$regions
    tgdf <- dplyr::filter(gene_df,ldmr==s_ldmr)
    scale_x <- c(min(tgdf$start),max(tgdf$end))

    updateSliderInput(session, "receive", value = scale_x,
                      min = scale_x[1], max = scale_x[2])
  })



  output$postPlotA <- renderPlotly({
    sv <-input$receive
    s_ldmr <-input$regions
    mods <-input$models
    annos <-input$annos
    gen_plots(sv,s_ldmr,mods,annos)
  }
  )

})
