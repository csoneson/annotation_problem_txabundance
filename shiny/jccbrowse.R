suppressPackageStartupMessages({
  library(shiny)
  library(ggplot2)
  library(Gviz)
  library(dplyr)
  library(tidyr)
  library(scatterpie)
  library(rtracklayer)
  library(ggrepel)
})

quant_methods <- c("hera", "kallisto", "RSEM", "Salmon", "SalmonCDS",
                   "SalmonSTAR", "SalmonKeepDup", "StringTie")

basedir <- "/Volumes/charlotte/annotation_problem_txabundance"

bigwig_files <- list(HAP1 = paste0(basedir, "/STARbigwig/20170918.A-WT_4_Aligned.sortedByCoord.out.bw"),
                     Cortex = paste0(basedir, "/STARbigwig/20151016.A-Cortex_RNA_Aligned.sortedByCoord.out.bw"))

gene_models <- readRDS(paste0(basedir,
                              "/reference/Gviz/Homo_sapiens.GRCh38.90_gviz_genemodels.rds"))$genemodels_exon

score_files <- list(HAP1 = readRDS(paste0(basedir,
                                          "/output/20170918.A-WT_4_combined_coverages_with_scores.rds")),
                    Cortex = readRDS(paste0(basedir,
                                            "/output/20151016.A-Cortex_RNA_combined_coverages_with_scores.rds")))

junction_tables <- lapply(score_files, function(f) {
  f$junctions %>% dplyr::filter(method %in% quant_methods)
})

transcript_tables <- lapply(score_files, function(f) {
  f$transcripts %>% dplyr::filter(method %in% quant_methods)
})

gene_tables <- lapply(score_files, function(f) {
  f$genes %>% dplyr::filter(method %in% quant_methods) %>%
    dplyr::select(gene, method, covOKfraction, intron_exon_ratio,
                  uniqjuncreads, mmjuncreads, uniqjuncfraction, score) %>%
    tidyr::spread(method, score)
})

all_genes <- unique(gene_models$gene)

options(ucscChromosomeNames=FALSE, envir=.GlobalEnv)

jccbrowse <- function(bigwig_files, score_files, all_genes, gene_models, 
                      junction_tables, transcript_tables, gene_tables) {
  options(ucscChromosomeNames=FALSE, envir=.GlobalEnv)
  p_layout <- shinydashboard::dashboardPage(
    skin = "purple", 
    
    shinydashboard::dashboardHeader(title = "JCC score"),
    
    shinydashboard::dashboardSidebar(
      radioButtons(inputId = "study", label = "Library", 
                   choices = c("HAP1", "Cortex"), selected = "HAP1"),
      
      selectizeInput(inputId = "gene", label = "Gene (Ensembl ID)",
                     choices = all_genes, selected = all_genes[1])
    ),
    
    shinydashboard::dashboardBody(
      shinydashboard::tabBox(
        width = 12,
        tabPanel(
          "Gene summary plot",
          fluidRow(plotOutput("gviz_plot", width = "100%", height = "400px")),
          fluidRow(
            column(width = 6, plotOutput("tpms_plot", width = "100%", height = "400px")),
            column(width = 6, plotOutput("junctions_plot", width = "100%", height = "400px"))
          )
        ),
        tabPanel(
          "HAP1 gene table",
          DT::dataTableOutput("hap1_gene_table")
        ),
        tabPanel(
          "Cortex gene table",
          DT::dataTableOutput("cortex_gene_table")
        )
      )
    )
  )
  
  server_function <- function(input, output, session) {
    options(ucscChromosomeNames=FALSE, envir=.GlobalEnv)

    muted <- c("#DC050C","#E8601C","#7BAFDE","#1965B0","#B17BA6",
               "#882E72","#F1932D","#F6C141","#F7EE55","#4EB265",
               "#90C987","#CAEDAB","#777777")
    
    values <- reactiveValues()
    
    observe({
      validate(need(
        tolower(input$gene) %in% tolower(all_genes),
        sprintf("Invalid gene")
      ))
      gm <- subset(gene_models, tolower(gene) == tolower(input$gene) | 
                     tolower(gene_name) == tolower(input$gene))
      gm <- subset(gm, gene == gene[1])  ## Select only one gene if there are many with the same name
      id <- unique(gm$gene_name)
      idshow <- paste0(id, " (", unique(gm$gene), ")")
      show_chr <- unique(seqnames(gm))[1]
      gm <- subset(gm, seqnames == show_chr)
      min_coord <- min(start(gm)) - 0.2*(max(end(gm)) - min(start(gm)))
      max_coord <- max(end(gm)) + 0.05*(max(end(gm)) - min(start(gm)))
      gm$transcript <- factor(gm$transcript, levels = unique(gm$transcript))
      
      ## Other features in the considered region
      gmo <- gene_models[overlapsAny(gene_models,
                                     GRanges(seqnames = show_chr,
                                             ranges = IRanges(start = min_coord,
                                                              end = max_coord),
                                             strand = "*"))]
      gmo <- gmo[!(gmo %in% gm)]
      gmo <- reduce(gmo)
      
      txs <- levels(gm$transcript)
      ncols <- nlevels(gm$transcript)
      cols <- colorRampPalette(muted)(ncols)
      names(cols) <- txs
      
      values$gm <- gm
      values$gmo <- gmo
      values$cols <- cols
      values$min_coord <- min_coord
      values$max_coord <- max_coord
      values$show_chr <- show_chr
      values$idshow <- idshow
      values$id <- id
      values$txs <- txs
    })
    
    # =========================== Gviz plot ================================= ##
    output$gviz_plot <- renderPlot({
      validate(need(
        tolower(input$gene) %in% tolower(all_genes),
        sprintf("Invalid gene")
      ))
      grtr <- GeneRegionTrack(values$gm, showId = TRUE, col = NULL, 
                              fill = values$cols[values$gm$transcript],
                              name = "", col.title = "black", 
                              background.title = "transparent", min.height = 10)
      grtr2 <- GeneRegionTrack(values$gmo, showId = TRUE, col = "black", fill = "white",
                               name = "", col.title = "black", showId = FALSE,
                               background.title = "transparent", min.height = 10)
      
      gtr <- GenomeAxisTrack()
      
      rtr <- DataTrack(range = bigwig_files[[input$study]],
                       type = "histogram",
                       chromosome = unique(seqnames(values$gm)),
                       col.title = "black",
                       fill = "grey",
                       col = "grey",
                       col.histogram = "grey",
                       fill.histogram = "grey",
                       cex.title = 0,
      )
      
      tracks <- c(rtr, gtr, grtr, grtr2)
      suppressWarnings(Gviz::plotTracks(tracks, chromosome = values$show_chr, 
                                        from = values$min_coord, to = values$max_coord, 
                                        main = values$idshow, 
                                        min.width = 0, min.distance = 0, collapse = FALSE))
    })
    
    # ============================ TPM plot ================================= ##
    output$tpms_plot <- renderPlot({
      validate(need(
        tolower(input$gene) %in% tolower(all_genes),
        sprintf("Invalid gene")
      ))
      txtmp <- transcript_tables[[input$study]] %>% dplyr::filter(gene == input$gene) %>%
        dplyr::mutate(transcript = factor(transcript))
      validate(need(
        nrow(txtmp) > 0,
        sprintf("No transcripts")
      ))
      ggplot(txtmp,
             aes(x = method, y = TPM, fill = transcript)) + 
        geom_bar(stat = "identity", position = "fill") + xlab("") + 
        scale_fill_manual(values = values$cols, name = "") + 
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 9),
              legend.text = element_text(size = 7),
              legend.position = "bottom") + 
        guides(fill = guide_legend(nrow = length(unique(txtmp$transcript)) %/% 6 + 1))
    })
    
    # ========================= Junction plot =============================== ##
    output$junctions_plot <- renderPlot({
      validate(need(
        tolower(input$gene) %in% tolower(c(values$gm$gene)),
        sprintf("Invalid gene")
      ))
      jl <- junction_tables[[input$study]] %>%
        dplyr::filter(gene == input$gene)
      for (tt in values$txs) {
        jl[[tt]] <- grepl(tt, jl$transcript)
      }
      jl[, txs] <- sweep(jl[, values$txs], 1, rowSums(jl[, values$txs]), "/")
      validate(need(
        nrow(jl) > 0,
        sprintf("No junctions")
      ))
      ggplot() + geom_abline(intercept = 0, slope = 1) + 
        geom_scatterpie(aes(x = scaled.cov, y = uniqreads, r = max(scaled.cov)/13), 
                        cols = values$txs, data = jl, color = NA) + 
        facet_wrap(~ methodscore, nrow = 2) + 
        coord_equal(ratio = 1) + 
        expand_limits(x = range(c(jl$scaled.cov, jl$uniqreads)), 
                      y = range(c(jl$scaled.cov, jl$uniqreads))) + 
        scale_fill_manual(values = values$cols, name = "") + 
        xlab("Scaled predicted coverage") + 
        ylab("Number of uniquely mapped reads") + 
        theme_bw() + theme(strip.text = element_text(size = 7),
                           legend.text = element_text(size = 7),
                           legend.position = "none")
    })
    
    # =========================== Gene tables =============================== ##
    output$cortex_gene_table <- DT::renderDataTable(
      DT::datatable(gene_tables[["Cortex"]] %>% 
                      dplyr::mutate(covOKfraction = signif(covOKfraction, 2),
                                    intron_exon_ratio = signif(intron_exon_ratio, 2),
                                    uniqjuncfraction = signif(uniqjuncfraction, 2)),
                    options = list(scrollX = TRUE))
    )
    output$hap1_gene_table <- DT::renderDataTable(
      DT::datatable(gene_tables[["HAP1"]] %>% 
                      dplyr::mutate(covOKfraction = signif(covOKfraction, 2),
                                    intron_exon_ratio = signif(intron_exon_ratio, 2),
                                    uniqjuncfraction = signif(uniqjuncfraction, 2)),
                    options = list(scrollX = TRUE))
    )
    
  }
  
  shinyApp(ui = p_layout, server = server_function)
}

