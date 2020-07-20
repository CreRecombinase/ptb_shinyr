#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(ldmap)
library(shiny)
library(plotly)
library(dplyr)

sel_models <- c("atac",
            "ataccdecctl",
            "atacctcmctl",
            "utme3treatedhic")
ldmrs <- as.character(as_ldmap_region(unique(fst::read_fst("data/gene_df.fst",columns = "ldmr")$ldmr)))
anno_df <- fst::read_fst("data/anno_tm_df.fst") %>%
    mutate(ld_snp=as_ldmap_snp(ld_snp),
           snp_pos=as_ldmap_snp(snp_pos),
           ldmr = as_ldmap_region(ldmr),
           gwas_snp=as_ldmap_snp(gwas_snp),
           pos=positions(gwas_snp),
           z=abs(beta)/se,
           p=pnorm(z,lower.tail=F)) %>% filter(model==stringr::str_to_lower(model))


ldmrs <- distinct(anno_df,ldmr,model) %>% filter(model=="null") %>% pull(ldmr) %>% unique() %>% as.character()
#ldmrs <- c("chr14:103012102_105001723", "chr1:21736588_23086883", "chr2:73174848_75630086",
# "chr3:121974097_123517768", "chr3:126214943_128194861", "chr3:137371083_139954597",
# "chr3:139954597_141339097", "chr3:154714218_156008700", "chr4:174264132_176570716",
# "chr5:156628700_158825698", "chr9:79471208_81079055", "chr3:4431647_6083183"
# )
ldmrs <- as_ldmap_region(ldmrs)
annotations <- unique(anno_df$anno)
models <- unique(anno_df$model[anno_df$model!="null"])

ildmr <- fst::read_fst("data/top_pip_gene.fst") %>% mutate(ldmr=as.character(as_ldmap_region(ldmr)))
ldrm <- magrittr::set_names(as_ldmap_region(ildmr$ldmr),paste0(ildmr$ldmr,":",as.character(ildmr$null)))

mmr <- c(103058996L, 104647235L)
# Define UI for application that draws a histogram
shinyUI(fluidPage(

    # Application title
    titlePanel("Susie Plots"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
        selectInput("regions",
                    "locus",choices=ldrm),
        sliderInput("receive","Region",min=mmr[1],max=mmr[2],value=mmr),
        checkboxGroupInput("models",label = "Models",choices=models,selected=sel_models),
        checkboxGroupInput("annos",label = "Annotatons",choices=annotations,selected=annotations),
        textInput(inputId = "gene",label="gene",value="HAND2",placeholder = "HAND2")
    ),
    mainPanel(
            plotlyOutput("postPlotA",height = "800px"),

    )
    )
))
