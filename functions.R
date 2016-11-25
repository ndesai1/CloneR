# LOADING DATA functions =================================================================================================

# This function read the feature of each sample
# Tab-separeted file
# columns: sample, gender (man=T, female=F), tumor content (i.e. 0.7)
read.config.file=function(filename){
  if(!is.null(filename)){
    x = read.table(filename, h=F, sep='\t')
    colnames(x) =c('sample','man','tc')
    x$tc[ which(is.na(x$tc)) ] = 1
    a = is.na(x$sample)
    b = is.na(x$man)
    if( sum(a)>0 | sum(b)>0 ){
      warning(paste("NA values found in 'Sample', 'Gender' fields in file ", filename," /nCorresponding samples will be removed from the analysis", sep="", coll=""))
      if( sum(a)>0 ){
        if(sum(a)!=nrow(x)){
          x = x[ which(!is.na(x$sample)),]
        }else{
          stop(paste("No sample correctly provided in file ", filename," /nCheck the format of the input file", sep="", coll=""))
        }
      }
      if( sum(b)>0 ){
        if(sum(b)!=nrow(x)){
          x = x[ which(!is.na(x$man)),]
        }else{
          stop(paste("No sample correctly provided in file ", filename," /nCheck the format of the input file", sep="", coll=""))
        }
      }
    }
    return(x)
  } else {
    stop("Sample_list file not found")
  }
}

# This function read the gene list
# Tab-separeted file
# columns: HUGO gene symbol, category
read.gene.file=function(filename){
  if(!is.null(filename)){
    x = read.table(filename, h=F, sep='\t')
    colnames(x) =c('gname','category')
    return(x)
  } else {
    stop("Gene_list file not found")
  }
}

# This function read the mutations of each sample
# Tab-separeted file
# columns: sample, chr, position , reference, variant ,  allele frequency, gene symbol, mutation type
read.mutation.file=function(filename){
  if(!is.null(filename)){
    x = read.table(filename, h=F, sep='\t')
    if(ncol(x)>=8){
      colnames(x) =c('sample','chrom','position','ref','varallele','freq','gname', 'alt.type') #'covRef', 'covMut',
      x$type="Mutation"
      if(length(grep("chr",x$chrom))==0) x$chrom=paste('chr',x$chrom, sep="", coll="")
      x = x[which(!is.na(x$freq)),]
      # if(sum(x$freq<=1)==nrow(x)) x$freq=x$freq*100
      chrm_no = sapply(x$chrom, toString)
      chrm_no = sapply(strsplit(chrm_no, split='r'), function(x) (x[2]))
      x$id = paste0(x$sample, ":" , chrm_no, ":", x$position,":", x$ref) #changed id from previous, need to change back later
      
      x$assignedCNV=NA
      x$CN_major  = 1
      x$CN_minor  = 1
      x$CN_a =2
      x$CN_t    = 2
      x$CN_normal = 2
    } else {
      stop("read.mutation.file: missing fields in mutation file ")
    }
  } else {
    stop("read.mutation.file: Mutation file not found")
  }

  x$assignedCNV = NA
  return(x)
}

# This function read the cnvs of each sample
# Tab-separeted file
# columns: sample, chr, start , end, cnv type, CN, CN raw
read.cnv.file=function(filename){
  if(!is.null(filename)){
    x = read.table(filename, h=F, sep='\t')
    if(ncol(x)==9){
      colnames(x) =c('sample','chrom','start','end','assignedCNV', 'CN_major', 'CN_minor', 'CN_a', 'CN_t')
      x$CN_normal = 2
      x$type="CNV"
      x$alt.type=x$assignedCNV
      x$id = with(x, paste(chrom,".",start,".",end, sep = "", coll = ""))
      return(x)
    }else if(ncol(x)==10){
      colnames(x) =c('sample','chrom','start','end','assignedCNV', 'CN_major', 'CN_minor', 'CN_a',  'CN_t', 'CN_normal')
      x$type="CNV"
      x$alt.type=x$assignedCNV
      x$id = with(x, paste(chrom,".",start,".",end, sep = "", coll = ""))
      return(x)
    }else{
      stop("read.cnv.file: missing fields in cnv file ")
    }
  } else {
    stop("read.CNV.file: CNV file not found")
  }
}

# This function determine whether a mutation is a SNV or an InDel
set.mutation.type = function(m){
  m$alt.type = "InDel"
  m$alt.type[which(nchar(m$varallele)==1 & nchar(m$ref)==1 & m$ref%in%c("A","C","G","T") & m$varallele%in%c("A","C","G","T"))] = "SNV"
  return(m)
}

# This function exclude samples without the associated SNP file from the CNV analysis
exclude.samples.without.SNPs = function(c, s){
  id = which( !is.na(s$snp_filename) )
  if( length(id)>0){
    ic = which( c$sample%in%s$sample[id] )
    if(length(ic)>0){
      c =  c[ ic,  ]
    }else{
      c = NULL
    }
  }
  return(c)
}

# This function read the cnvs of each sample
# Tab-separeted file
# columns: id, chrom, position, freq.N, freq.T
read.SNP.file=function(filename){
  if(!is.null(filename)){
    x = read.table(filename, h=F, sep='\t')
    if(ncol(x)==5 & nrow(x)>0){
      colnames(x) =c('id','chrom','position','freq.N','freq.T')
      id.na = which(is.na(x$freq.T))
      if(length(id.na)>0) x = x[-id.na,]
      if(sum(x$freq.T<=1)==nrow(x)){
        x$freq.T = x$freq.T*100
        x$freq.N = x$freq.N*100
      }
      if(length(grep("chr",x$chrom))==0){
        x$chrom = paste( "chr", x$chrom, sep="", coll="")
      }
      return(x)
    }else{
      warning(paste("read.snp.file: ,",filename," missing fields in the file", sep="", coll=""))
    }
  } else {
    warning(paste("read.SNP.file: ,",filename," heterozygous SNP file not found", sep="", coll=""))
  }
}

# This function plot assign CNV status to those mutations falling in a CNV region
get.CNV.status = function(m, c=NULL, chr=chr_levels){
  require(GenomicRanges)
  if(!is.null(c) & !is.null(m)){
    im  = with(m, GRanges(chrom, IRanges(position,position) )); seqlevels(im, force=T)=chr
    ic  = with(c, GRanges(chrom, IRanges(start,end) )); seqlevels(ic, force=T)=chr
    ov  = findOverlaps(im,ic)
    if(length(ov)>0){
      m$assignedCNV[queryHits(ov)] = c$assignedCNV[subjectHits(ov)]
      m$CN_major[queryHits(ov)]    = c$CN_major[subjectHits(ov)]
      m$CN_minor[queryHits(ov)]    = c$CN_minor[subjectHits(ov)]
      m$CN_a[queryHits(ov)] = c$CN_a[subjectHits(ov)]
      m$CN_t[queryHits(ov)]      = c$CN_t[subjectHits(ov)]
      m$CN_normal[queryHits(ov)]      = c$CN_normal[subjectHits(ov)]

    }
    #     else{
    #       warning("No mutations in CNV regions")
    #     }
  }
  return(m)
}

# This function  assign genes to CNV regions
get.genes.in.CNV.regions = function(c, annotation_filename = "GeneLists/Agilent_genes.tsv"){
  require(GenomicRanges)
  if(!is.null(c)){
    ag = read.delim(annotation_filename,stringsAsFactors=F, h=F)
    if(nrow(ag)>0){
      ig = with(ag,GRanges(V1,IRanges(V2,V3) ) )
      ic  = with(c, GRanges(chrom, IRanges(start,end) ))
      ov  = findOverlaps(ig,ic)
      c$gname=NA
      if(length(subjectHits(ov))>0) c$gname[subjectHits(ov)] = ag$V4[queryHits(ov)]
      # c = as.data.frame( cbind('gname'=ag$V4[queryHits(ov)], c[subjectHits(ov),]), stringsAsFactors=F )
    }else{
      warning("NO gene annotations provided. Excluding CNVs from the analysis")
      c = NULL
    }
  }
  return(c)
}


# CORRECTED FREQUENCY functions  =================================================================================

# SNVS *********************************************
get.tc.correction.somatic = function( obs, tc, CNt, CNn=2){
  return( min(

    obs * ( 1 + (  ( CNn*(1-tc) )/( CNt * tc) ) )  ,

    1))
  # 100))
}

get.cor.tumor.content = function(x,y){

  for(i in 1:nrow(x))  x$freq.tc[i] = get.tc.correction.somatic(x$freq[i], y$tc, x$CN_t[i], 2 )

  if(y$man){
    idX = which(x$chrom=="chrX")
    idY = which(x$chrom=="chrY")
    if(length(idX)>0) for(i in idX) x$freq.tc[i] = get.tc.correction.somatic( x$freq[i], y$tc, x$CN_t[i], CNn=1)
    if(length(idY)>0) for(i in idY) x$freq.tc[i] = get.tc.correction.somatic( x$freq[i], y$tc, x$CN_t[i], CNn=1)
  }
  return(x)
}

get.CNV.proportion = function(CN_t, CN_a, CN_n  = 2 ){
  CN_prop = ( CN_t - CN_n ) / ( CN_a - CN_n )
  is_na = which(is.na(CN_prop))
  if(length(is_na)>0) CN_prop[is_na]=0
  return( CN_prop )
}

get.CNV.clonality.for.SNVs=function(x,y){
  tc = y$tc
  #get total absolute number of alleles in copy number aberrant cells
  #x$CN_a = with(x, CN_major + CN_minor )
  #get the average copy in tumor cells corrected for tumor content
  #x$CN_t = with(x, (CN_raw - (CN_normal * ( 1-tc )) ) / tc )

  #calculate the clonality of cells with CNVs (CNV_prop)
  x$CNV_prop = mapply(get.CNV.proportion,x$CN_t, x$CN_a, x$CN_normal)

  #get a list of possible n-values that will be used to determine clonality in subclonal populations
  for(i in 1:nrow(x)) {
    if ((x$CN_major[i]+x$CN_minor[i]) == x$CN_a[i])
    {
      x$n_vals[i] = paste(unique(c(0,1,x$CN_minor[i], x$CN_major[i], max(x$CN_a[i] -1, 0), x$CN_a[i]) ), collapse = ",")
    }else{
      x$n_vals[i]=paste(unique(c(0,1, x$CN_minor[i], max(x$CN_a[i]-1, 0) , x$CN_a[i])), collapse=",")
    }
  }  
  return(x)
}


# CLOANLITY functions  =================================================================================

# SNVS *********************************************

#This function calculates clonality for SNVs in subclonal copy number regions
get.clonality.subclonal.CNV = function(freq.tc, CNV_prop, CN_t, n_vals, CN_h = 2){
  
  n_values = as.integer( unlist(strsplit(n_vals, ",") ) )
  clonal_vals = F2 = NULL;
  n_vals = NULL;
  y = 1 - CNV_prop;
  
  for (n in n_values){
    
    f_2 = ((CN_t * freq.tc) - (n*CNV_prop) )/(CN_h * y )
    
    F2 = c(F2,
           f_2)
    
    if (f_2 <=1.05 & f_2 >0.01) # if n = 0, it means there are no SNV hits in CNV aberrant cells,
      # and therefore clonality is equal to y*f_2*2
    {
      m = ifelse(n==0, 0, 1)
      clonality = (m*CNV_prop) + (y * f_2 * CN_h)
      if (clonality<=1.05 & clonality >0){
        clonal_vals = c(clonal_vals, clonality)
        n_vals = c(n_vals, n)
      }
    }
    
    
  }
  if(length(clonal_vals)==0){
    clonal_vals = freq.tc*CN_t
  }
  
  cat("n\t",n_values, "\nF_2\t", F2, '\n')
  cat("n_accepted\t", n_vals, '\n')
  
  cat("\nclonalities\t", clonal_vals,'\n')
  clonality = mean(unique(clonal_vals))
  
  cat("clonality is", clonality, '\n')
  return(min(clonality, 1))
}

# I leave this function for you Nikita, please modify the one above
NIKITA.get.clonality.subclonal.CNV = function(freq.tc, CNV_prop, CN_t, n_vals, CN_h = 2){
  n_values = as.integer( unlist(strsplit(n_vals, ",") ) )
  clonal_vals = NULL;

  for (n in n_values){

    f_2 = ( (CN_t * freq.tc) - (n*CNV_prop) )/(CN_h * (1 - CNV_prop) )

    if (f_2 <=1 & f_2 >0) # if n = 0, it means there are no SNV hits in CNV aberrant cells,
      # and therefore clonality is equal to y*f_2*2
    {
      if (n==0) m = 0
      if (n>0) m = 1
      clonality = ( m*CNV_prop ) + ( f_2 * CN_h * (1 - CNV_prop) )

      if (clonality<=1 & clonality >0){
        clonal_vals = c(clonal_vals, clonality)
      }
    }


  }
  clonality = mean(unique(clonal_vals))
  return(clonality)
}


#reads information from dataframe y and calculates clonality
#y dataframe with CNV clonality information and CN information for each SNV
#columns: id, freq.tc, gname, alt.type, tc, CN_t, CNV_prop, n_vals, actual clonality
get.SNV.clonality = function(x){
  for(i in 1:nrow(x))
    if( (x$CNV_prop[i] <=0.9 & x$CNV_prop[i] >0.1)){
      x$cell[i] =  get.clonality.subclonal.CNV(   x$freq.tc[i]
                                                  , x$CNV_prop[i]
                                                  , x$CN_t[i]
                                                  , x$n_vals[i]
                                                  , x$CN_normal[i] )
    }else{
      # for all SNVs where copy number is not subclonal, clonality = freq.tc*CN_t
      x$cell[i] = min(x$freq.tc[i]*x$CN_t[i], 1)
    }
  return(x)
}

# CNVS *********************************************

# This function reads integrated information for CNV regions and outputs dataframe with CNV clonality
get.CNV.clonality=function(x,y){

  tc = y$tc
  #get total absolute number of alleles in copy number aberrant cells
  #x$CN_a = with(x, CN_major + CN_minor )
  #get the average copy in tumor cells corrected for tumor content
  #x$CN_t = with(x, (CN_t - (CN_normal * ( 1-tc )) ) / tc )

  #calculate the clonality of cells with CNVs (CNV_prop)
  x$cell = get.CNV.proportion(x$CN_t, x$CN_a, x$CN_normal)

  return(x)
}


# Output ===================================================================================================

# This function generate the dataset for the clonality.plot
prepare.dataset = function(m=NULL, c=NULL){
  mcol = c('sample','alt.type','cell','gname','type','id','assignedCNV','CN_major', 'CN_minor', 'CN_normal','CN_a','CN_t') #,'freq','freq.tc'
  if(!is.null(m) & is.null(c)){
    return(m[,mcol])
  }else if(is.null(m)  & !is.null(c)){
    return(c[,mcol])
  }else{
    return( as.data.frame(rbind(m[,mcol],c[,mcol]), stringsAsFactors=F) )
  }
}

prepare.dataset.1set = function(c){
  mcol = c('sample','type','cell','alt.type','id','assignedCNV','CN_major', 'CN_minor', 'CN_normal','CN_a','CN_t') #,'freq','freq.tc'
  if(!is.null(c)){
    return( c[,mcol] )
  } else {
    return(NULL)
  }
}


# This function marks the gene of interest
set.gene.category = function(m, gl=NULL){
  if(!is.null(m)){
    m$category = NA
    if(!is.null(gl)){
      m$category = gl$category[ match( m$gname, gl$gname) ]
    }
  }
  return(m)
}

# This function create the outout folder
create.output.folder = function (s, odir, title=NULL) {
  if(!dir.exists(odir)) dir.create(odir)

  job     = paste0('Job',length(grep(".Job", list.dirs(path = odir, recursive =F)))+1)
  timer    = format(Sys.time(), "%Y-%m-%d.%H_%M_%S")

  # if(!is.null(title)) job=paste(title,job,sep='.')

  analisys = paste(job,timer, sep=".")
  analisys = paste0(odir,"/", analisys )

  dir.create( analisys )
  for(i in s$sample ) dir.create( paste0( analisys, "/", i))
  return(analisys)
}


#this function get the clone composition
get.clone.composition = function(x, upper=.80, lower=.35){
  x$upper=upper
  x$lower=lower
  y = z = NULL
  if(!is.null(x)){
    if(sum(x$type=="Mutation")>0){
      y = ddply(subset(x, type=="Mutation"), .(sample), summarise,
                n = length(cell),
                n_monoclonal = sum(cell>=unique(upper)),
                n_biclonal   = sum(cell<unique(upper) & cell>=unique(lower)),
                n_polyclonal = sum(cell<unique(lower)))
    }

    if(sum(x$type=="CNV")>0){

      z = ddply(unique(subset(x, type=="CNV")[,c(1:4,6,11,12)]), .(sample), summarise,
                n = length(cell),
                n_monoclonal = sum(cell>=unique(upper)),
                n_biclonal   = sum(cell<unique(upper) & cell>=unique(lower)),
                n_polyclonal = sum(cell<unique(lower)))
    }

    if(!is.null(y) & !is.null(z)){
      y$n = y$n + z$n
      y$n_monoclonal = y$n_monoclonal + z$n_monoclonal
      y$n_biclonal = y$n_biclonal + z$n_biclonal
      y$n_polyclonal = y$n_polyclonal + z$n_polyclonal
    }else if(is.null(y) & !is.null(z)){
      y = z
    }else if(is.null(y) & is.null(z)){
      return(NULL)
    }

    y$monoclonal = y$n_monoclonal / y$n
    y$biclonal   = y$n_biclonal / y$n
    y$polyclonal = y$n_polyclonal / y$n

    code = c("M","B","P"); names(code)=c('monoclonal','biclonal','polyclonal')
    y$composition = code[names(which.max(y[,c('monoclonal','biclonal','polyclonal')]))]
  }
  return(y)
}


#this function plot a bar chart representing the clone composition
clone.composition.plot= function(x, cl=color_clone_composition){
  names(cl)=c("monoclonal",'biclonal','polyclonal')
  x$composition = factor(x$composition, levels=c("M","B","P"))

  if(!is.null(x)){
    mapper = as.list(x$composition)
    names(mapper) = x$sample
    map_labeller <- function(variable,value){
      return(mapper[value])
    }

    m = melt(x[,c('sample','polyclonal','biclonal','monoclonal')], id.vars = c("sample")) #1,4:2
    colnames(m)[2] = 'composition'
    p=ggplot(m, aes(x=sample,y=value,fill=composition))+
      geom_bar(width=0.5, stat="identity")+
      geom_segment(aes(x=-Inf,xend=-Inf,y=0,yend=1),col="black")+
      ylab("Alterations (%)")+xlab("")+
      scale_fill_manual(values=cl,
                        guide = guide_legend(title = NULL),
                        labels=c("Clonality<35%",'35%<Clonality<80%','Clonality>80%'))+
      scale_y_continuous(labels=c("0","25","50","75","100"))+
      facet_grid(sample ~ ., labeller=map_labeller)+
      theme_cloneR()+theme(panel.background = element_blank(),
                           legend.position  = "top",
                           legend.key=element_rect(size=1.5, color='white'),
                           legend.text      = element_text(color="black",size=10),
                           axis.text        = element_text(color="black",size=10),
                           axis.text.y      = element_blank(),
                           axis.ticks.y     = element_blank(),
                           strip.background = element_rect(fill=NA, colour=NA),
                           strip.text.y     = element_text(angle=0,size=12, colour="black")
      )+
      coord_equal(1/0.1)+
      coord_flip()+
      geom_bar(width=0.5, stat="identity",color="black", show.legend=FALSE)

    return(p)
  } else{
    return(NULL)
  }
}

clone.composition.plot.overall= function(x, cl = color_clone_composition){
  if(!is.null(x)){
    mapper = as.list(x$composition)
    names(mapper) = x$sample
    map_labeller <- function(variable,value){
      return(mapper[value])
    }
    names(cl)=c("monoclonal",'biclonal','polyclonal')
    x$composition = factor(x$composition, levels=c("M","B","P"))
    m = melt(x[,c('sample','polyclonal','biclonal','monoclonal')], id.vars = c("sample")) #1,4:2
    m$sample = factor(as.character(m$sample), x$sample[order(x$composition,x$monoclonal)] ) # ,x$biclonal,x$polyclonal
    colnames(m)[2] = 'composition'
    p=ggplot(m, aes(x=sample,y=value,fill=composition))+
      geom_bar(width=0.9, stat="identity")+
      geom_segment(aes(x=-Inf,xend=-Inf,y=0,yend=1),col="black")+
      ylab("Alterations (%)")+xlab("")+
      scale_fill_manual(values=cl,
                        guide = guide_legend(title = NULL),
                        labels=c("Clonality<35%",'35%<Clonality<80%','Clonality>80%'))+
      scale_y_continuous(labels=c("0","25","50","75","100"))+
      theme_cloneR()+theme(panel.background = element_blank(),
                           legend.position  = "left",
                           legend.key       = element_rect(size=.5, color='white'),
                           legend.text      = element_text(color="black",size=10),
                           axis.text        = element_text(color="black",size=10),
                           axis.text.y      = element_text(size=rel(.75))
      )+
      geom_bar(width=0.9, stat="identity",color="black", show.legend = FALSE)+
      coord_equal(1/0.05)+coord_flip()
    return(p)
  } else{
    return(NULL)
  }
}

clone.composition.plot.overall.ply= function(x, cl = color_clone_composition){
  if(!is.null(x)){
    mapper = as.list(x$composition)
    names(mapper) = x$sample
    map_labeller <- function(variable,value){
      return(mapper[value])
    }
    names(cl)=c("monoclonal",'biclonal','polyclonal')
    x$composition = factor(x$composition, levels=c("M","B","P"))
    m = melt(x[,c('sample','polyclonal','biclonal','monoclonal')], id.vars = c("sample")) #1,4:2
    m$sample = factor(as.character(m$sample), x$sample[order(x$composition,x$monoclonal)] ) # ,x$biclonal,x$polyclonal
    m$value = m$value*100
    colnames(m)[2:3] = c('composition','percentage')


    return(m)
  } else{
    return(NULL)
  }
}

bar.composition = function( x , cl = color_clone_composition){
  nn = c('Monoclonal','Biclonal','Polyclonal'); names(nn) = c('M','B','P')

  y = table(x$composition)
  y = melt(y)
  y$Var1  = factor(y$Var1, levels = names(nn))
  y$p = 100*(y$value/sum(y$value))
  y$label = paste0(round(y$p, 1), '% (',y$value,")")
  ggplot(y, aes(y=p, x=Var1, fill=Var1))+
    geom_bar(stat="identity",width = .9, color="black")+
    geom_text(aes(y=p+2.5, label=label))+
    geom_segment(aes(x=-Inf,xend=-Inf,y=0,yend=100),col="black")+
    scale_x_discrete(labels=nn[as.character(y$Var1)])+
    scale_fill_manual(values=cl)+
    ylab("Samples (%)")+xlab("")+
    theme_classic()+
    theme( legend.position  = "none",
           legend.key       = element_rect(size=1.5, color='white'),
           legend.text      = element_text(color="black",size=12),
           axis.text        = element_text(color="black",size=10),
           axis.text.x      = element_text(angle=90))+coord_equal(0.05)
}

# These functions generate axes
base_breaks_x = function(x, br, la, reverse=T){
  d <- data.frame(y=-Inf, yend=-Inf, x=0, xend=1)
  if(reverse){
    return(list(geom_segment(data=d, aes(x=x, y=y, xend=xend, yend=yend), inherit.aes=FALSE), scale_x_reverse(breaks=br, labels=la)))
  }else{
    return(list(geom_segment(data=d, aes(x=x, y=y, xend=xend, yend=yend), inherit.aes=FALSE), scale_x_continuous(breaks=br, labels=la)))
  }
}

base_breaks_x2 = function(x, br, la, reverse=T){
  d <- data.frame(y=-Inf, yend=-Inf, x=0, xend=x)
  if(reverse){
    return(list(geom_segment(data=d, aes(x=x, y=y, xend=xend, yend=yend), inherit.aes=FALSE), scale_x_reverse(breaks=br, labels=la)))
  }else{
    return(list(geom_segment(data=d, aes(x=x, y=y, xend=xend, yend=yend), inherit.aes=FALSE), scale_x_continuous(breaks=br, labels=la)))
  }
}

base_breaks_y = function(m, br, la ){
  d <- data.frame(x=Inf, xend=Inf, y=0, yend=m)
  list(geom_segment(data=d, aes(x=x, y=y, xend=xend, yend=yend), inherit.aes=FALSE), scale_y_continuous(breaks=br, labels=la))
}

base_breaks_y2 = function(m, br, la ){
  d <- data.frame(x=-Inf, xend=-Inf, y=0, yend=m)
  list(geom_segment(data=d, aes(x=x, y=y, xend=xend, yend=yend), inherit.aes=FALSE), scale_y_continuous(breaks=br, labels=la))
}

check.alterations.per.category = function( y ){
  # duplicates
  iy = paste(y$cell,'.',y$id, sep="", coll="")
  iy = which(!duplicated(iy))
  if(length(iy)>0){
    y = y[iy, ]
  }

  is_na = which(is.na(y$gname))
  if(length(is_na)>0) y$gname[is_na] = y$id[is_na]


  # cat        = rep(0,4); names(cat) = c("SNV","InDel","Gain","Loss");
  # alternative = c("InDel","SNV","Loss","Gain"); names(alternative) = c("SNV","InDel","Gain","Loss");
  #
  # cat[ names(table(y$Alteration_type)) ] = table(y$Alteration_type)
  #
  # y$original_category = y$Alteration_type
  #
  # id = names(which(cat<=5 & cat>0))
  #
  # if( sum(id%in%c("SNV","InDel"))==2){
  #   y$Alteration_type[ which(y$original_category %in%c("SNV","InDel"))] = 'SNV'
  # } else if (sum(id%in%c("Gain","Loss"))==2){
  #   y$Alteration_type[ which(y$original_category %in%c("Gain","Loss"))] = 'Gain'
  # } else {
  #   for(i in id){
  #     y$Alteration_type[ which(y$original_category == i)] =  alternative[ i  ]
  #   }
  # }
  # y$original = c("No","Yes")[as.numeric( y$Alteration_type==y$original_category )+1 ]
  return(y)
}

prepare.dataset.to.density.plot = function(y){
  tmp = subset(y, Alteration_type%in%c("Gain",'Loss'))

  if(nrow(tmp)>0){
    cx = ddply(tmp, .(id,Alteration_type,cell), summarise, gname=paste(gname,collapse=",") )
  }

  mut = subset(y, Alteration_type%in%c("SNV","InDel"))

  dd = NULL

  if(nrow(mut)>0 & nrow(tmp)>0){
    dd = as.data.frame( rbind(mut[,c('Alteration_type','cell')], cx[,c('Alteration_type','cell')]) , check.names=F)
  }else if(is.null(tmp)){
    dd = as.data.frame( mut[,c('Alteration_type','cell')] , check.names=F)
  }else if(is.null(mut)){
    dd = as.data.frame( cx[,c('Alteration_type','cell')] , check.names=F)
  }

  dd$Alteration_type=factor(dd$Alteration_type, levels=  c("SNV","InDel","Gain","Loss"))

  m = with(dd, table(Alteration_type, cell))
  ix = apply(m,1, function(x) sum(x!=0)==1)
  if(sum(ix)>0){
    for(i in names(ix)[ix]){
      a = subset(dd, Alteration_type==i)
      b = a[1,]
      b$cell=b$cell-1
      if(nrow(a)==1){
        dd = as.data.frame(rbind(dd,a,b),stringsAsFactors=F)
      }else{
        dd = as.data.frame(rbind(dd,b),stringsAsFactors=F)
      }
    }
  }
  return(dd)
}

# This function plot derives and plots the density distributions of clonality for SNVs, InDels, amplifications, and deletions separately
density.plot=function(x,adj=5, names=T, fill_cols=color_density_plot){
  if(!is.null(x)){
    colnames(x)[2]='Alteration_type'
    px = NULL

    # SELECT ONLY MUTATIONS THAT DO NOT UNDERGO CNVS

    px = subset(x, type=="Mutation" | type=="CNV")[,c('cell','Alteration_type','gname','category','id')]

    if(nrow(px)<2){
      cat("density.plot: LESS THEN 2 ALTERATIONS TO DRAW THE PLOT")
      return(NULL)
    }

    px = check.alterations.per.category( px )

    # PREPARE Dataset

    dd = prepare.dataset.to.density.plot(px)

    # Density Plot
    p =
      ggplot(dd, aes(x=cell,fill=Alteration_type)) + geom_density(alpha=.5, adjust=4) + aes(y = ..count..) +
      geom_vline(xintercept=.80, linetype="dashed", color="grey")+
      geom_vline(xintercept=.35, linetype="dashed", color="grey")+
      scale_fill_manual(values=fill_cols,
                        guide = guide_legend(title = NULL)
      )

    # leg = g_legend(p)

    # PLOT POINTS FOR DRIVERS

    gg  = ggplot_build(p)$data[[1]]
    gga = unique(ggplot_build(p)$plot[[1]]$Alteration_type) #     ptmp = ggplot(dd, aes(x=cell,fill=Alteration_type)) + geom_histogram() #     gtmp =  ggplot_build(ptmp)$data[[1]]
    ymax = max(gg$y)
    ymed = ymax/2
    lmax = ifelse( ymax>=0.5, round(ymax), 1)
    lmed = ifelse(lmax==1, 0.5, ifelse( ymed>0.5, round(ymed), 0.5))

    #drivers=subset(px, !is.na(category) )
    drivers = NULL

    if(!is.null(drivers)>0){

      p.drivers = plot.drivers(drivers, gg, gga )

      drvrs=ddply(drivers, .(cell,Alteration_type), summarise, gname=paste(gname, collapse=","))

      plvs=c("SNV","InDel","Gain","Loss"); names(plvs)=c("SNV","InDel","Gain","Loss")
      plvs = plvs[which(names(plvs)%in%gga)]

      p.drvrs=rbind()
      for(j in unique((drvrs$Alteration_type))){
        dd1 = subset(gg,group == which(plvs == j) )
        y1  =  subset(drvrs,(Alteration_type)==j )
        yy1=xxx1=c()
        for (i in y1$cell){
          xxx1 = c(xxx1, dd1$x[which(abs(dd1$x-i)==min(abs(dd1$x-i)))][1])
          yy1  = c(yy1,  dd1$y[which(abs(dd1$x-i)==min(abs(dd1$x-i)))][1])
        }
        if(!is.null(xxx1)){
          y1$cell=xxx1
          y1$y=yy1
          p.drvrs=rbind(p.drvrs,y1)
        }
      }

      p.drvrs      = as.data.frame(p.drvrs, stringsAsFactors=F)

      p +
        geom_point(data=p.drvrs, aes(x = cell, y=y), colour="black", show.legend = F)+


        ishape = c(19,17); names(ishape)=c("Yes","No")

      if(names==T){
        p= p + geom_point(data=p.drvrs, aes(x = cell, y=y, shape=original), colour="black", show.legend = F)+
          geom_text_repel(data=p.drvrs, aes(x=cell, y=y, label = gname))#        (data=p.drivers, aes(x = cell,y=y,label=gname, family=""),  size=rel(3.5), hjust=rel(.75), vjust=rel(-.25), show.legend = F)
        # scale_shape_manual(values=ishape[levels(p.drivers$original)], guide = 'none')
        if(sum(p.drivers$nx>0)>0){
          p = p +   geom_segment(data=subset(p.drivers,nx>0), aes(x = cell,y=y,xend=nx,yend=ny))
        }
      }

    }
    #----

    p= p +theme_cloneR()+
      theme(    legend.position  = "top",
                legend.text      = element_text(color="black",size=10),
                axis.text        = element_text(size=rel(0.85)),
                axis.title.x     = element_text(size=rel(1)),
                axis.title.y     = element_text(size=rel(1)),
                rect             = element_blank())+
      base_breaks_y(ymax, c(0,ymed,ymax), as.character(c(0,lmed,lmax)))+
      base_breaks_x(gg$x, c(0,.35,.80,1), c("0","35","80","100")) +
      ylab('Expected Alterations')+xlab('Alteration clonality (%)')
    return(p)
  } else {
    return(NULL)
  }
}

get.clonality.gene.of.interest = function(x){
  drivers = NULL
  drivers=subset(x, !is.na(category) )
  if(nrow(drivers)>0){
    hc = ddply( drivers, .(sample, gname), summarise, cell=cell[which.max(cell)])
    r = unique(hc[,2])
    c = unique(hc[,1])
    m = matrix(0, nr=length(r), nc=length(c), dimnames=list(r,c))
    for(i in 1:nrow(hc)) m[as.character(hc$gname[i]), as.character(hc$sample[i]) ] = hc$cell[i]
    m = melt(m)
    colnames(m) = c('Gene','Sample','Alteration_clonality')
    return(m)
  }else{
    return(NULL)
    warning("NO gene of interest provided: upload your gene list")
  }
}

#this function get the heatmap of select genes
heatmap.genes = function( x ){
  p = NULL
  drivers = NULL
  drivers=subset(x, !is.na(category) )
  if(nrow(drivers)>0){
    hc = ddply( drivers, .(sample, gname), summarise, cell=cell[which.max(cell)])
    r = unique(hc[,2])
    c = unique(hc[,1])
    m = matrix(0, nr=length(r), nc=length(c), dimnames=list(r,c))
    for(i in 1:nrow(hc)) m[as.character(hc$gname[i]), as.character(hc$sample[i]) ] = hc$cell[i]
    m = melt(m)
    colnames(m)[3] ='Alteration_clonality'
    myPalette <- colorRampPalette(c("white",brewer.pal(9, "Reds")), space="Lab")
    col = myPalette(100)
    mx = max(m$Alteration_clonality)
    mn = min(m$Alteration_clonality)
    if(mn==mx) mn=0
    p = ggplot(m, aes(y=Var2,x=Var1, fill=Alteration_clonality))+
      geom_tile(color="black")+
      coord_equal()+
      # scale_fill_gradientn(colours = myPalette(100), name='Alteration clonality', limits = c(mn,mx), breaks =c(mn,mx), labels=as.character(c(mn,mx)))+  #,title.position = "top",label.position="bottom"
      scale_fill_gradientn(colours = myPalette(100), name='Alteration clonality', limits = c(0,100), breaks =c(0,50,100), labels=as.character(c(0,50,100)))+  #,title.position = "top",label.position="bottom"
      theme_cloneR()+
      theme(legend.position  = 'top',
            legend.title     = element_text(color="black",size=12),#rel(1.5)
            legend.text      = element_text(color="black",size=12),
            axis.text        = element_text(color="black",size=10),
            axis.text.x      = element_text( hjust=1,angle=45), #vjust=0.5,angle=45,
            axis.text.y      = element_blank(),
            axis.ticks.y     = element_blank(),
            strip.background = element_rect(fill=NA, color="black"),
            strip.text       = element_text(size=12)) +
      xlab("") +ylab("")
  }else{
    warning("NO gene of interest provided: upload your gene list")
    return(NULL)
  }
  return(p)
}

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  legend
}
savePlot = function (analisys, sample, composition, density_plot, ht_genes, gene_list = NULL) {

  ccR='cR_report'
  ccR_version='CloneR v.1 2016'

  pdf(file=paste( analisys, "/", sample, "/", ccR, ".pdf",sep="", coll="" ),width = unit(8.27, "inches"), height = unit(11.69, "inches"))

  pushViewport(viewport(layout=grid.layout(8, 3, widths=c(0.05, 1, 0.1), heights=c(.05,.01,.15,.01, 0.4,.01,.2 ,.05))))

  grid.text( sample , vp = viewport(layout.pos.row = 1, layout.pos.col = 2), gp=gpar(fontface='bold'), just="center")

  grid.text("A.", vp = viewport(layout.pos.row = 2, layout.pos.col = 1), gp=gpar(fontface='bold'), just="left")

  if(!is.null(composition)){
    print( composition, vp=viewport(layout.pos.row=3,layout.pos.col=2))
  }

  grid.text("B.", vp = viewport(layout.pos.row = 4, layout.pos.col = 1), gp=gpar(fontface='bold'), just="left")

  if(!is.null(density_plot)){
    print(density_plot, vp=viewport(layout.pos.row=5,layout.pos.col=2:3))
    #     density_plot$vp = viewport(layout.pos.row=5,layout.pos.col=2:3)
    #     grid.draw(density_plot)
  }

  grid.text("C.", vp = viewport(layout.pos.row = 6, layout.pos.col = 1), gp=gpar(fontface='bold'), just="left")

  if(!is.null(gene_list)){
    if(!is.null(ht_genes)){
      print(ht_genes,    vp=viewport(layout.pos.row=7,layout.pos.col=2))
    }
  }

  grid.text(ccR_version, vp = viewport(layout.pos.row = 8, layout.pos.col = 2))
  dev.off()
}


global_report = function(analysis, composition ){
  if(nrow(composition)>1){
    ccR='cR_overall_clone_composition'
    ccR_version='CloneR v.1 2016'

    pdf(file=paste(analysis, "/", ccR, ".pdf",sep="", coll="" ),width = unit(8.27, "inches"), height = unit(11.69, "inches"))

    pushViewport(viewport(layout=grid.layout(4, 5, widths=c(0.01, .48, 0.02, .48, 0.01), heights=c( .05, .45,.45, .05))))

    grid.text("A.", vp = viewport(layout.pos.row = 1, layout.pos.col = 1), gp=gpar(fontface='bold'), just="left")

    p = clone.composition.plot.overall(composition)
    leg = g_legend(p)
    print(p+theme(legend.position='none') , vp=viewport(layout.pos.row=2:3,layout.pos.col=2))

    grid.text("B.", vp = viewport(layout.pos.row = 1, layout.pos.col = 3), gp=gpar(fontface='bold'), just="left")

    print( bar.composition(composition), vp=viewport(layout.pos.row=2,layout.pos.col=4))

    leg$vp =viewport(layout.pos.row=3,layout.pos.col=4)
    grid.draw(leg)
    grid.text(ccR_version, vp = viewport(layout.pos.row = 4, layout.pos.col = 3), just="left")

    dev.off()
  }
}
