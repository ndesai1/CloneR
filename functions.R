# Loading Configuration #####
# print.logo = function(){
#   cat("\n")
#   cat("██████╗██╗      ██████╗ ███╗   ██╗███████╗██████╗\n")
#   cat("██╔════╝██║     ██╔═══██╗████╗  ██║██╔════╝██╔══██╗\n")
#   cat("██║     ██║     ██║   ██║██╔██╗ ██║█████╗  ██████╔╝\n")
#   cat("██║     ██║     ██║   ██║██║╚██╗██║██╔══╝  ██╔══██╗\n")
#   cat("╚██████╗███████╗╚██████╔╝██║ ╚████║███████╗██║  ██║\n")
#   cat("╚═════╝╚══════╝ ╚═════╝ ╚═╝  ╚═══╝╚══════╝╚═╝  ╚═╝\n")
#   cat(  "\n")
# }

# print.logo()

# cat("Loading configuration.R ...")
source("config.R")
# cat("Done\n")

# Loading Functions #####

# cat("Loading Functions ...")
# Loading data functions =======

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
      if(sum(x$freq<=1)==nrow(x)) x$freq=x$freq*100
      x$assignedCNV=NA
      x$CN = 2
      x$CN_raw = 2

    } else {
      stop("read.mutation.file: missing fields in mutation file ")
    }
  } else {
    stop("read.mutation.file: Mutation file not found")
  }

  x$assignedCNV = NA
  return(x)
}

# This function determine whether a mutation is a SNV or an InDel
set.mutation.type = function(m){
  m$alt.type = "InDel"
  m$alt.type[which(nchar(m$varallele)==1 & nchar(m$ref)==1 & m$ref%in%c("A","C","G","T") & m$varallele%in%c("A","C","G","T"))] = "SNV"
  return(m)
}

# This function read the cnvs of each sample
# Tab-separeted file
# columns: sample, chr, start , end, cnv type, CN, CN raw
read.cnv.file=function(filename){
  if(!is.null(filename)){
    x = read.table(filename, h=F, sep='\t')
    if(ncol(x)==7){
      colnames(x) =c('sample','chrom','start','end','assignedCNV', 'CN', 'CN_raw')
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

# Tumour content correction =======

# This function corrects for the tumour content of the sample
# get.cor.tumor.content = function(x,y){
#   x$freq.tc = round(x$freq/y$tc)
#   if(sum(x$freq.tc>100)>0) x$freq.tc[ which(x$freq.tc>100) ] = 100
#   return(x)
# }

get.tc.correction.somatic = function( obs, tc, CNt, CNn=2){
  return( min(

    obs * ( 1 + (  ( CNn*(1-tc) )/( CNt * tc) ) )  ,

    100))
}

get.cor.tumor.content = function(x,y){

  for(i in 1:nrow(x))  x$freq.tc[i] = get.tc.correction.somatic(x$freq[i], y$tc, x$CN_raw[i], 2 )

  if(y$man){
    idX = which(x$chrom=="chrX")
    idY = which(x$chrom=="chrY")
    if(length(idX)>0) for(i in idX) x$freq.tc[i] = get.tc.correction.somatic( x$freq[i], y$tc, x$CN_raw[i], CNn=1)
    if(length(idY)>0) for(i in idY) x$freq.tc[i] = get.tc.correction.somatic( x$freq[i], y$tc, x$CN_raw[i], CNn=1)
  }
  return(x)
}


# Clonality assessment =======

# This function prepare the mutation dataset for the clonal assessment
prepare.mutation.dataset = function(x){
  x$id=with(x,paste(chrom,".",position,".",ref,".",varallele,sep="",coll=""))  # mm = min(x$freq.tc)
  y  = subset(x,freq.tc<=50)
  x  = z = subset(x,freq.tc>50)
  if(nrow(x)>0){
    x$freq.tc = x$freq.tc - 50    # x$cov.tc = round( x$freq.tc * x$cov.tc,0)
    x$position = x$position+1
    z$freq.tc = 50    # z$cov.tc = round( z$freq.tc * z$cov.tc,0)
    tmp = rbind(y,x,z)
    tmp = subset(tmp,freq.tc>=1 & freq.tc<=50)
  }else{
    tmp=y
  }
  tmp
}

# This function prepare the mutation dataset for the clonal assessment in case of male samples
prepare.mutation.dataset.male=function(x){
  tmp = tmp2 = NULL
  tmp = subset(x, !chrom%in%c("chrX","chrY") )
  if(nrow(tmp)>0)  tmp = prepare.mutation.dataset(tmp)
  tmp2    = subset(x, chrom%in%c("chrX","chrY"))
  if(nrow(tmp2)>0){
    tmp2$id = with(tmp2,paste0(chrom,".",position,".",ref,".",varallele))
    tmp2$freq.tc = tmp2$freq.tc/2
    if(!is.null(tmp)){
      tmp = rbind(tmp,tmp2)
    }else{
      tmp = tmp2
    }
  }
  tmp
}



# This function calculate the clonality
get.clonality = function(x,y){
  if(!is.null(x)){
    clonality = x$freq.tc*2
    if(!y$man){
      clonality[which(x$chrom=="chrX")] =  x$freq.tc[which(x$chrom=="chrX")]
    }
    x$cell = clonality
  }
  return(x)
}

# Auxialiary functions
up = function(x) {
  x = x[ which(x>50) ]
  return(median(x,na.rm=T))
}

dw = function(x) {
  x = x[ which(x<50) ]
  return(median(x,na.rm=T))
}

# This function assign the SNP frequncies to CNVs
get.hetero.SNPs.freq.in.CNV.regions = function(c,s){
  if(!is.null(c) & !is.null(s)){
    snp = read.SNP.file(s$snp_filename)
    if(nrow(c)>0 & nrow(snp)>0){
      ig = with(snp, GRanges(chrom,IRanges(position,position) ) )
      ic = with(c, GRanges(chrom, IRanges(start,end) ))
      ov = NULL
      ov = findOverlaps(ig,ic)
      if(length(ov)>0 & subjectLength(ov)>0){
        tmp = as.data.frame( cbind( c[subjectHits(ov),c("id",'assignedCNV','CN','CN_raw')], freq = snp[queryHits(ov),c('freq.T')]), stringsAsFactors=F )
        tmp[,'freq'] = as.numeric(tmp[,'freq'])
        res = ddply(tmp, .(id), summarise, freq.T.up = up( freq ), freq.T.dw = dw( freq ))
        c = cbind(c, res[match(c$id, res$id),c('freq.T.up','freq.T.dw')])
        ix = which(is.na(c$freq.T.up) & is.na(c$freq.T.dw))
        if(length(ix)>0){
          if(length(ix)!=nrow(c)){
            c=c[-ix,]
          } else{
            warning("NO heterozygous SNPs in the regions")
            c=NULL
          }
        }
      }else{
        warning("NO heterozygous SNPs in CNV regions. Excluding CNVs from the analysis")
        c = NULL
      }
    }else{
      warning("cnv_list or snp_list empty")
    }
  }
  return(c)
}


# correzione per il tc (germinali)

get.tc.correction.germline = function( obs, tc, CNt, CNn=2, Fg = 50){

  if(obs<60 & obs>40){
    return(50)
  }else{

    return( min(

      ( obs * ( 1 + (  ( CNn*(1-tc) )/( CNt * tc) ) ) ) - (Fg * (  ( CNn*(1-tc) )/( CNt * tc) ) ) ,

      100))
  }
}

#================================================================
# chiamata di pyclone per capire se sono clonali o sottoclonali
#================================================================

# This function assign the clonality of CNVs

get.cor.tumor.content.CNV = function(x,y){
  if(is.null(x)) return(NULL)
  x$freq.tc = NA
  x$freq = NA
  w = rbind()
  # z = subset(x, assignedCNV=='Gain')
  z = x
  ix = which(is.na(z$freq.T.up))
  if(length(ix)>0){
    z$freq.T.up[ix] = 100 -  z$freq.T.dw[ ix ]
    ix = which(is.na(z$freq.T.up))
    if(length(ix)>0) z=z[-ix,]
    if (is.null(z)) return(x)
  }

  x$freq = z$freq.T.up

  if(nrow(z)>0){

    for(i in 1:nrow(z)) z$freq.tc[i] = get.tc.correction.germline( z$freq.T.up[i], y$tc, z$CN_raw[i])

    if(y$man){
      idX = which(z$chrom=="chrX")
      idY = which(z$chrom=="chrY")
      if(length(idX)>0) for(i in idX) z$freq.tc[i] = get.tc.correction.germline( z$freq.T.up[i], y$tc, z$CN_raw[i], CNn=1)
      if(length(idY)>0) for(i in idY) z$freq.tc[i] = get.tc.correction.germline( z$freq.T.up[i], y$tc, z$CN_raw[i], CNn=1)
    }
    w = rbind(w,z)
  }
   return(w)
}


# This function assign the clonality of CNVs
get.clonality.CNV = function(x, y ){
  if(!is.null(x)){
    if(nrow(x)>0){

      x$cell = ifelse( x$freq.tc<60 & x$freq.tc>40 ,
                       100 ,
                       (round(x$freq.tc, 0) - 50 )*2
      )

      ix = which(x$cell>100)
      if(length(ix)>0) x$cell[ix] = 100


      ix = which(is.na(x$cell))
      if(length(ix)>0){
        if(length(ix)!=nrow(x)){
          x=x[-ix,]
        } else{
          warning("NO heterozygous SNPs in the regions")
          x=NULL
        }
      }
    }else{
      warning("Impossible assessing CNV clonality")
      # x = NULL
    }
  }
  return(x)
}

# This function plot assign CNV status to those mutations falling in a CNV region
get.CNV.status = function(m, c=NULL, chr=chr_levels){
  if(!is.null(c) & !is.null(m)){
    im  = with(m, GRanges(chrom, IRanges(position,position) )); seqlevels(im, force=T)=chr
    ic  = with(c, GRanges(chrom, IRanges(start,end) )); seqlevels(ic, force=T)=chr
    ov  = findOverlaps(im,ic)
    if(length(ov)>0){
      m$assignedCNV[queryHits(ov)] = c$assignedCNV[subjectHits(ov)]
      m$CN[queryHits(ov)] = c$CN[subjectHits(ov)]
      m$CN_raw[queryHits(ov)] = c$CN_raw[subjectHits(ov)]
    }
    #     else{
    #       warning("No mutations in CNV regions")
    #     }
  }
  return(m)
}

# This function plot assign genes to CNV regions
get.genes.in.CNV.regions = function(c, annotation_filename = "GeneLists/Agilent_genes.tsv"){
  if(!is.null(c)){
    ag = read.delim(annotation_filename,stringsAsFactors=F, h=F)
    if(nrow(ag)>0){
      ig = with(ag,GRanges(V1,IRanges(V2,V3) ) )
      ic  = with(c, GRanges(chrom, IRanges(start,end) ))
      ov  = findOverlaps(ig,ic)
      c = as.data.frame( cbind('gname'=ag$V4[queryHits(ov)], c[subjectHits(ov),]), stringsAsFactors=F )
    }else{
      warning("NO gene annotations provided. Excluding CNVs from the analysis")
      c = NULL
    }
  }
  return(c)
}

get.CNV.without.mutations = function(c, m, chr=chr_levels){
  if(!is.null(c)){
    c$with_mutations=F
    if(!is.null(m)){
      im  = with(m, GRanges(chrom, IRanges(position,position) )); seqlevels(im, force=T)=chr
      ic  = with(c, GRanges(chrom, IRanges(start,end) )); seqlevels(ic, force=T)=chr
      ov  = findOverlaps(ic,im)
      if(length(ov)>0){
        c$with_mutations[queryHits(ov)] = T
      }
    }
  }
  return(c)
}


# This function generate the dataset for the clonality.plot
prepare.dataset = function(m=NULL, c=NULL){
  mcol = c('sample','type','cell','alt.type','gname','id','assignedCNV','CN','CN_raw','freq','freq.tc')
  if(!is.null(m) & is.null(c)){
    return(m[,mcol])
  }else if(is.null(m)  & !is.null(c)){
    return(c[,mcol])
  }else{
    return( as.data.frame(rbind(m[,mcol],c[,mcol]), stringsAsFactors=F) )
  }
}

prepare.dataset.1set = function(c){
  mcol = c('sample','type','cell','alt.type','gname','id','assignedCNV','CN','CN_raw','freq','freq.tc')
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


# Output =======

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
get.clone.composition = function(x, upper=80, lower=35){
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

#this function get the clone composition
get.clone.composition.OLD = function(x, upper=80, lower=35){
  x$upper=upper
  x$lower=lower
  y = NULL
  if(!is.null(x)){
    y = ddply(x, .(sample), summarise,
              monoclonal = sum(cell>=unique(upper))/length(cell),
              biclonal   = sum(cell<unique(upper) & cell>=unique(lower))/length(cell),
              polyclonal = sum(cell<unique(lower))/length(cell))
    code = c("M","B","P"); names(code)=c('monoclonal','biclonal','polyclonal')
    y$composition = code[names(which.max(y[,2:4]))]
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
  d <- data.frame(y=-Inf, yend=-Inf, x=0, xend=100)
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

word_cloud = function (pdr) {
  # browser()
  if(nrow(pdr)>1){
    pdf(file.path(tempdir(), "Rplots.pdf"))
    plot(y = pdr$y, x=pdr$cell, type="n")
    pd             = as.data.frame(wordlayout(x=pdr$cell, y=pdr$y, words = pdr$code, xlim=c(0,Inf), ylim=c(0,Inf), cex=1.5),  stringsAsFactors=F, check.names =F)
    dev.off()
    tmp            = strsplit(rownames(pd), "\\:")
    pd$gname       = sapply(tmp, function(x) x[3] )
    pd$Alteration_type    = sapply(tmp, function(x) x[2] )
    rm(tmp)
    for (i in 1:nrow(pd)) {
      xl <- pd[i, 1]
      yl <- pd[i, 2]
      w  <- pd[i, 3]
      h  <- pd[i, 4]
      if (pdr$cell[i] < xl || pdr$cell[i] > xl + w || pdr$y[i] < yl || pdr$y[i] > yl + h) {
        pdr$nx[i] <- xl + 0.45 * w
        pdr$ny[i] <- yl + .5 * h
      }
    }
  }
  return(pdr)
}

plot.drivers = function (drvrs, gg_builder, a ) {
  # browser()
  drvrs=ddply(drvrs, .(cell,Alteration_type, original), summarise, gname=paste(gname, collapse=","))

  plvs=c("SNV","InDel","Gain","Loss"); names(plvs)=c("SNV","InDel","Gain","Loss")
  plvs = plvs[which(names(plvs)%in%a)]

  p.drvrs=rbind()
  for(j in unique((drvrs$Alteration_type))){
    dd1 = subset(gg_builder,group == which(plvs == j) )
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

  p.drvrs$code = with( p.drvrs, paste(1:nrow(p.drvrs),":", Alteration_type, ":",gname, sep="",coll=""))
  p.drvrs$nx=0
  p.drvrs$ny=0

  p.drvrs = word_cloud(p.drvrs)
  return(p.drvrs)
}

check.alterations.per.category = function( y ){
  iy = paste(y$cell,'.',y$id, sep="", coll="")
  iy = which(!duplicated(iy))
  if(length(iy)>0){
    y = y[iy, ]
  }
  cat        = rep(0,4); names(cat) = c("SNV","InDel","Gain","Loss");
  alternative = c("InDel","SNV","Loss","Gain"); names(alternative) = c("SNV","InDel","Gain","Loss");
  cat[ names(table(y$Alteration_type)) ] = table(y$Alteration_type)
  y$original_category = y$Alteration_type
  id = names(which(cat<=5 & cat>0))
  if( sum(id%in%c("SNV","InDel"))==2){
    y$Alteration_type[ which(y$original_category %in%c("SNV","InDel"))] = 'SNV'
  } else if (sum(id%in%c("Gain","Loss"))==2){
    y$Alteration_type[ which(y$original_category %in%c("Gain","Loss"))] = 'Gain'
  } else {
    for(i in id){
      y$Alteration_type[ which(y$original_category == i)] =  alternative[ i  ]
    }
  }
  y$original = c("No","Yes")[as.numeric( y$Alteration_type==y$original_category )+1 ]
  return(y)
}

prepare.dataset.to.density.plot = function(y){
  tmp = subset(y, Alteration_type%in%c("Gain",'Loss'))
  cx = ddply(tmp, .(id,Alteration_type,cell), summarise, gname=paste(gname,collapse=",") )
  dd = NULL
  dd = as.data.frame( rbind(subset(y, Alteration_type%in%c("SNV","InDel"))[,c('Alteration_type','cell')], cx[,c('Alteration_type','cell')]) , check.names=F)
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
    colnames(x)[4]='Alteration_type'
    px = NULL

    # SELECT ONLY MUTATIONS THAT DO NOT UNDERGO CNVS

    # if(sum(x$type=="Mutation")>0){
    # px = subset(x, (type=="Mutation" & is.na(assignedCNV)) | type=="CNV")[,c('cell','Alteration_type','gname','category','id')]
    px = subset(x, type=="Mutation" | type=="CNV")[,c('cell','Alteration_type','gname','category','id')]
    # }else{
    # px = x[,c('cell','Alteration_type','gname','category','id')]
    # }

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
      geom_vline(xintercept=80, linetype="dashed", color="grey")+
      geom_vline(xintercept=35, linetype="dashed", color="grey")+
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

    drivers=subset(px, !is.na(category) )

    if(nrow(drivers)>0){
      p.drivers = plot.drivers(drivers, gg, gga )

      ishape = c(19,17); names(ishape)=c("Yes","No")

      if(names==T){
        p= p + geom_point(data=p.drivers, aes(x = cell, y=y, shape=original), colour="black", show.legend = F)+
          geom_text(data=p.drivers, aes(x = cell,y=y,label=gname, family=""),  size=rel(3.5), hjust=rel(.75), vjust=rel(-.25), show.legend = F)
        scale_shape_manual(values=ishape[levels(p.drivers$original)], guide = 'none')
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
      base_breaks_x(gg$x, c(0,35,80,100), c("0","35","80","100")) +
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

# cat("Done\n")


