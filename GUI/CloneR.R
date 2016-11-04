sample_filename='../Example/test.sample.txt'
mutation_filename='../Example/test.mutations.txt'
cnv_filename='../Example/test.cnvs.txt'
gene_list_filename='../GeneLists/Actionable_genes.tsv'
outdir='~/tmp'

source("../functions.R")
source("../config.R")


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


