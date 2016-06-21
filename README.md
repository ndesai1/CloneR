![CloneR](www/Cloner_logo.png =250x)
=======

[CloneR] is a **R** application that evaluates the clone composition of a tumour. Starting from molecular data (i.e. results od whole exome/genome sequences and/or genome-wide SNP arrays) and sample specific information (i.e. gender, tumour purity), it estimates:
- the allele frequency of somatic alterations (i.e. point mutations and Copy Number Variant (CNV) regions) correcting for tumour purity and CN status;
- the alteration clonality of each somatic alteration;
- the clone compostion of the tumour (i.e. Monoclonal, Biclonal and Polyclonal)
For each sample, CloneR estimetes the density distribution of the alterations clonalities. If specifed, CloneR produces an heatmap of the alteration clonality of a list of selected genes. 

## Usage
The easiest way to run CloneR is to install **shiny** package (and the required dependencies) in R, and use the function `runGithub()`. See the example below,
```
install.packages("shiny")
install.packages("shinyjs")

install.packages("plotly")

runGitHub("CloneR","your_username")
```

## Contributors

CloneR has been designed by Dr **Matteo Cereda** and Prof **Francesca Ciccarelli**. 

Main developer: Matteo Cereda. 

Contributions are always welcome!

## License

Please read the [Licence](LICENSE) first