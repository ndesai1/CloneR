source("functions.R")
source("config.R")
library(shiny)
library(DT)


options(shiny.maxRequestSize = 2000*1024^2)
options(scipen=3)


cloneR= function(  sample_filename
                 , mutation_filename
                 , cnv_filename
                 , gene_list_filename
                 , ann_filename
                 , outdir){
  sample_list = gene_list = mutations = cnvs = snps = samples = sample_list = mutation_list = cnv_list = genes_of_interest = NULL

  # ==== Load Configuration File =====
  samples = read.config.file( sample_filename )
  sample_list = dlply( samples, .(sample) )

  # ==== Genes of interest File =====
  gene_list = read.gene.file( gene_list_filename )

  # ==== Load mutations =====
  if(!is.null(mutation_filename)) mutations = read.mutation.file(mutation_filename)

  # ==== Load CNVs =====
  if(!is.null(cnv_filename)) cnvs = read.cnv.file(cnv_filename)

  if(is.null(mutations) & is.null(cnvs)) stop("\nYou MUST provide Mutation file and/or CNV file")
  if(!is.null(mutations)) mutations = mutations[which(mutations[,1]%in%samples[,1]),]
  if(!is.null(cnvs)) cnvs = cnvs[which(cnvs[,1]%in%samples[,1]),]
  cat("Number of samples:\t",nrow(samples))
  cat("\nNumber of samples with mutation data:\t",length(unique(mutations$sample)))
  cat("\nNumber of samples with CNV data:\t",length(unique(cnvs$sample)))
  cat("\n")

  # --------------------------------------------
  # CREATING RESULT FOLDER
  # --------------------------------------------

  analysis = create.output.folder(samples, outdir, title )

  # --------------------------------------------
  # DATA PROCESSING AND ASSESMENT OF CLONALITY
  # --------------------------------------------

  # incProgress(0.2)

  if(!is.null(cnvs)) cnv_list = dlply( cnvs, .(sample))

  # ==== Clonality assessment - MUTATIONS =====

  if(!is.null(mutations)){
    cat("Assessing clonality of somatic mutations...")
    mutation_list = dlply( mutations, .(sample) )

    # Assign CNV
    if(!is.null(mutation_list) & !is.null(cnv_list)){
      mutation_list = mapply( get.CNV.status , mutation_list, cnv_list[ names(mutation_list)], SIMPLIFY = F)
    }


    # Frequency Correction by TC and CN
    mutation_list = mapply( get.cor.tumor.content, mutation_list, sample_list[ names( mutation_list ) ], SIMPLIFY = F)

    mutation_list = mapply( get.CNV.clonality.for.SNVs, mutation_list, sample_list[ names( mutation_list ) ], SIMPLIFY = F)

    mutation_list = lapply( mutation_list, get.SNV.clonality)

    mutation_list        = mutation_list[ names(sample_list) ]
    names(mutation_list) = names(sample_list)

    cat("Done\n")
  }



  # ==== Clonality assessment - CNV =====

  # incProgress(0.4)

  if(!is.null(cnvs)){
    cat("Assessing clonality of CNVs...\n")

    cnv_list = cnv_list[ names(sample_list) ]
    cnv_list = mapply( get.CNV.clonality, cnv_list, sample_list, SIMPLIFY = F)

    cnv_list_region = cnv_list
    cat("\nDone\n")
  }


  # === Assign genes to CNV regions ====

  if(!is.null(cnv_list)){
    cat("Annotating genes in CNV regions...\n")
    if(file.exists(ann_filename)){
      cnv_list = lapply_pb(cnv_list, get.genes.in.CNV.regions, annotation_filename = ann_filename )
    }else{
      cnv_list = lapply(cnv_list, get.genes.in.CNV.regions )
    }
    cat("\nDone\n")
  }

  # ==== Merge Mutation and CNV datasets =====

  # incProgress(0.6)

  if(!is.null(mutation_list) & !is.null(cnv_list)){
    mutation_list = mutation_list[ names(sample_list) ]
    cnv_list      = cnv_list[ names(sample_list) ]
    dataset_list  = mapply( prepare.dataset, mutation_list, cnv_list, SIMPLIFY = F)
  }else if(!is.null(mutation_list) & is.null(cnv_list)){
    dataset_list = lapply(mutation_list, prepare.dataset.1set)
  }else{
    dataset_list = lapply(cnv_list, prepare.dataset.1set)
  }

  # ==== Set genes of interest =====

  dataset_list = lapply( dataset_list, set.gene.category, gl=gene_list )

  # --------------------------------------------
  # GENERATING OUTPUT and REPORT
  # --------------------------------------------

  # ==== Clone Composition =====

  cat("Plotting clone composition...")
  clone_comp_list   = lapply(dataset_list, get.clone.composition)
  clone_composition = lapply(clone_comp_list, clone.composition.plot)
  cat("Done\n")

  # ==== Density plot =====

  cat("Plotting clonality ditributions...\n")
  density_plot_list = lapply(dataset_list, density.plot)
  cat("\nDone\n")

  # ==== Clonality heatmap =====

  cat("Plotting clonality of gene of interest...")
  if(!is.null(gene_list)){
    heatmap_genes = lapply_pb( dataset_list, heatmap.genes )
    genes_of_interest = lapply( dataset_list, get.clonality.gene.of.interest )
    genes_of_interest = do.call("rbind",genes_of_interest[ !sapply(genes_of_interest,is.null) ])
  }
  cat("Done\n")

  # ==== Exporting plots =====

  # incProgress(0.8)

  cat("Exporting plots...")

  pb = txtProgressBar(min = 0, max = length(samples$sample), style = 3)
  i = 0
  for( i_sample in samples$sample ){
    setTxtProgressBar(pb, i)
    savePlot(analysis, i_sample,
             clone_composition[[ i_sample ]],
             density_plot_list[[ i_sample ]],
             heatmap_genes[[i_sample]],
             gene_list)
    i = i+1
  }

  cat("\nDone\n")

  # ==== REPORTS =====

  cat("Generating reports...")

  for( i_sample in samples$sample ){
    write.table(dataset_list[[i_sample]][,c('sample','alt.type','id','cell','gname','assignedCNV','category')],
                file = paste0(analysis, "/", i_sample,"/",i_sample,".composition.tsv" ), row.names = F, col.names = T, quote = F, sep = "\t")
  }

  dataset = do.call( 'rbind', dataset_list )
  write.table(dataset[,c('sample','alt.type','id','cell','gname','assignedCNV','category')],
              file = paste(analysis, "clonality.tsv", sep='/' ), row.names = F, col.names = T, quote = F, sep = "\t")

  composition = do.call("rbind",clone_comp_list)
  write.table(composition, file = paste(analysis, "composition.tsv", sep='/'), row.names = F, col.names = T, quote = F, sep = "\t")

  global_report(analysis, composition)

  if(!is.null(gene_list)){
    if(nrow(genes_of_interest)>0){
      write.table(genes_of_interest, file = paste( analysis, "clonality_genes_of_interest.tsv", sep='/'), row.names = F, col.names = T, quote = F, sep = "\t")
    }
  }
  cat("Done\n")

  # ==== IMAGE =====

  cat("Saving R image...")


  sInfo = sessionInfo()
  save(list = ls(all = TRUE),file=paste( analysis, "cloneR_image.RData",sep="/"))
  cat("Done\n")

  composition$p_monoclonal =  round(composition$monoclonal,3)*100
  composition$p_biclonal   =  round(composition$biclonal,3)*100
  composition$p_polyclonal =  round(composition$polyclonal,3)*100


  ret = list(comp=composition, dataset=dataset[,c('sample', 'alt.type','gname','id','cell','cell',"category")])
  ret$comp    = unrowname(ret$comp)


  ret$dataset         = unrowname(ret$dataset)
  # ret$dataset$CN_raw  = round(ret$dataset$CN_raw, 2)
  # ret$dataset$freq    = round(ret$dataset$freq,2)
  # ret$dataset$freq.tc = round(ret$dataset$freq.tc,2)
  ret$dataset$cell    = round(ret$dataset$cell,2)

  # colnames(ret$dataset) = c('Sample','Alteration','HUGO symbol','id','Allele Frequency (AF)', 'estimated AF','Copy Number (CN)','CN_raw','Clonality','Gene category')
  colnames(ret$dataset) = c('Sample','Alteration','HUGO symbol','id','Clonality','Gene category')



}



# Define server logic required to draw a histogram
shinyServer(function(input, output) {

  current = reactiveValues(dataset =NULL, composition =NULL)


      # CLONER = eventReactive(input$sumbit,{
      #   # withProgress( message = "Running CloneR", value = 0.1,
      #   #               {
      #       # })
      # })


      CLONER = eventReactive(input$submit, {
              title = sample_filename = mutation_filename = cnv_filename = snp_folder = NULL

              title              = input$title
              outdir             = input$outdir
              sample_filename    = input$sample_file$datapath
              mutation_filename  = input$mutation_file$datapath
              cnv_filename       = input$cnv_file$datapath
              # snp_folder         = input$snp_folder

              gene_list_filename = ifelse(is.null(input$gene_coordinate_file), "GeneLists/Actionable_genes.tsv", input$gene_coordinate_file$datapath)
              ann_filename       = ifelse(is.null(input$list_gene_file), "GeneLists/Agilent_genes.tsv", input$list_gene_file$datapath)

              print(title)
              print(outdir)
              print(snp_folder)
              print(gene_list_filename)
              print(ann_filename)

              if(!is.null(outdir) & !is.null(sample_filename) & ( !is.null(mutation_filename) | (!is.null(cnv_filename)) ) ){

                cloneR(sample_filename
                       , mutation_filename
                       , cnv_filename
                       # , snp_folder
                       , gene_list_filename
                       , ann_filename
                       , outdir)

              }

      })


      # TAB Clone Composition
      output$global_report = renderPlotly({
        results = CLONER()

        current$dataset = results[["dataset"]]
        current$composition = results[["comp"]]

        browser()
        m = clone.composition.plot.overall.ply(current$composition)
        print(m)
        # plot_ly(m, y=sample, x=percentage,  type='bar', orientation = "h", color=composition, colors=color_clone_composition)
        plot_ly( y=m$sample, x=m$percentage,  type='bar', orientation = "h", color=m$composition, colors=color_clone_composition_pl) %>%
        layout(p, barmode = 'stack',
                   yaxis=list(title='', tickfont=list(family = "Arial, sans-serif")),
                   xaxis=list(title='Alterations (%)', titlefont=list(family = "Arial, sans-serif"), tickfont=list(family = "Arial, sans-serif"))
        )



      })

      # # TAB Summary
       output$composition   = DT::renderDataTable({
         comp = current$composition[,c("sample","n","n_monoclonal","n_biclonal","n_polyclonal","p_monoclonal","p_biclonal","p_polyclonal","composition" )]
         colnames(comp) = c('Sample','Alterations','Monoclonal (n)','Biclonal (n)','Polyclonal (n)','Monoclonal (%)','Biclonal (%)', 'Polyclonal (%)', 'Compostion')

         DT::datatable(comp, options = list(orderClasses = TRUE, pageLength = 100))

       })

      # TAB Alteration Clonalities
       output$dataset       = DT::renderDataTable({


         DT::datatable(current$dataset, options = list(orderClasses = TRUE, pageLength = 100))
       })

})
