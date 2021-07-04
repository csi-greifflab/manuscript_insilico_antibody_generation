# Title     : TODO
# Objective : TODO
# Created by: rahmadakbar
# Created on: 2020-09-11

# import stuff
library(themeakbar)
library(VennDiagram)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(data.table)
library(Rtsne)
library(stringr)
library(wesanderson)
library(factoextra)
library(gtools)
library(ggridges)
library(RCy3)
library(kebabs)

theme_set(themeakbar())
my_spectral <- colorRampPalette(brewer.pal(8,'Spectral'))(14)
outfigdir = 'figures'
darjeer = colorRampPalette(wes_palette('Darjeeling1', n=5))(14)
set2 = colorRampPalette(brewer.pal(8,'Set2'))(14)
my_spectral=darjeer
my_spectral=set2

# my_base_color = c('#78C8A3','#95C4EA','#EE2A7B', '#F7BF16')
my_base_color = c('#78C8A3','#EE2A7B','#95C4EA', '#F7BF16')
my_base_color2 = colorRampPalette(my_base_color)(14)
my_spectral=my_base_color2

nativecol = my_spectral[10]
igemcol = my_spectral[14]
igemcol_large= 'orange'


outpdf = function(infile, tag, width=8, height=8){
    # assumes infile has extensins (.csv etc)
    # uses the first bit of the file as name and add the tag
    # uses ggsave, opens pdf plot after saving
    inname = strsplit(tail(strsplit(infile, '/')[[1]], n=1), '\\.')[[1]][1]
    outname = sprintf('%s/%s_%s.pdf', outfigdir,inname,tag)
    ggsave(outname, width=width, height = height)
    system(sprintf('open %s', outname))
    print(sprintf('opening %s', outname))

}


outpng = function(infile, tag, width=8, height=8){
    # assumes infile has extensins (.csv etc)
    # uses the first bit of the file as name and add the tag
    # uses ggsave, opens pdf plot after saving
    inname = strsplit(tail(strsplit(infile, '/')[[1]], n=1), '\\.')[[1]][1]
    outname = sprintf('%s/%s_%s.png', outfigdir,inname,tag)
    ggsave(outname, width=width, height = height)
    system(sprintf('open %s', outname))
    print(sprintf('opening %s', outname))

}

fixfive =function(df){
    # fix 5E94
    df$antigen[df$antigen=='FiveE94'] = '5E94'
    return(df)
}

fixfive2 =function(df){
    # fix 5E94
    df$antigen[df$antigen=='FiveE94'] = '5E94'
    df$antigen2[df$antigen2=='FiveE94'] = '5E94'
    return(df)
}


fixfive3 =function(df){
    # fix 5E94
    df$antigen[df$antigen=='FiveE94'] = '5E94'
    df$crossantigen[df$crossantigen=='FiveE94'] = '5E94'
    return(df)
}

# check cross transfer affinities

crosstransfer_affinity = function(){
    infile = 'eleven_outfiles/crosstransfer_native_generated_nsamples_opulent_sampled.csv'
    df = read_csv(infile)
    df = fixfive3(df)
    df$crossantigen_sample = sprintf('%s_%s', df$crossantigen, df$sample)
    # meddf = df %>% group_by(antigen, crossantigen, ntrainseq, crossantigen_sample) %>% summarize(medval = median(Energy))
    meddf = df %>% group_by(antigen, crossantigen, ntrainseq) %>% summarize(medval = median(Energy))
    print(meddf)
    df$crossstatus[df$crossstatus=='cross'] = '+T'
    df$crossstatus[df$crossstatus=='nocross'] = '-T'
    ggplot(data=df) +
        geom_density_ridges(mapping=aes(x=Energy, y=crossantigen, fill=crossstatus)) +
        geom_text(data=meddf, mapping=aes(x=medval, y=crossantigen, label=medval)) +
        facet_grid(ntrainseq~antigen, scale='free') +
        labs(y=' ', x='Affinity (Energy)') +
        theme(axis.title=element_text(size=20),
                legend.title=element_blank(),
                legend.text=element_text(size=20)) +
        scale_fill_manual(values=c('red', 'gray'))
    outpng('fig', 'S3', width=20, height=10)
}


# run stuff
crosstransfer_affinity()
