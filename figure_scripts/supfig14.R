# Title     : TODO
# Objective : TODO
# Created by: rahmadakbar
# Created on: 6/8/21


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
library(Biostrings)
library(Polychrome)
theme_set(themeakbar())
my_spectral <- colorRampPalette(brewer.pal(8,'Spectral'))(14)
outfigdir = 'figures'
darjeer = colorRampPalette(wes_palette('Darjeeling1', n=5))(14)
set2 = colorRampPalette(brewer.pal(8,'Set2'))(14)
set3 = colorRampPalette(brewer.pal(8,'Set2'))(20)
my_spectral=darjeer
my_spectral=set2
my_spectral2 = set3


# my_base_color = c('#78C8A3','#95C4EA','#EE2A7B', '#F7BF16')
my_base_color = c('#78C8A3','#EE2A7B','#95C4EA', '#F7BF16')
my_base_color3 = colorRampPalette(my_base_color)(20)
my_base_color2 = colorRampPalette(my_base_color)(14)
my_spectral=my_base_color3

nativecol = my_spectral[10]
igemcol = my_spectral[14]


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


get_overlap = function(infile1){
    # infile1  = 'mason/outfiles/mHER_H3_AgPos_11300_11300_corpus_geminio_113000000_len10.csv'
    infile2 = 'mason/outfiles/mason_30.csv'
    df1 = read_csv(infile1)[1:6000000,]
    df2 = read_csv(infile2)
    print(df1)
    print(df2)
    overlap = intersect(df1$gseq, df2$seq)
    print(overlap)
}


get_position_freq = function(infile){
    # infile = ('mason/outfiles/mHER_H3_AgPos_11300_11300_corpus_geminio_113000000_pbind.csv')
    df = read_csv(infile)
    df = df[1:6000000,]
    str_set = AAStringSet(df$gseq)
    pssm = consensusMatrix(str_set, as.prob=TRUE)
    print(str_set)
    pssmdf = as.data.frame(pssm)
    colnames(pssmdf) = sprintf('P%s', c(4:13))
    pssmcum = cumsum(pssmdf)
    pssmcum$aa = rownames(pssm)
    pssmdf$aa = rownames(pssm)
    gdf = gather(pssmdf, 'aapos', 'freq', -aa)
    gdf$aapos = factor(gdf$aapos, levels=mixedsort(unique(gdf$aapos)))
    gdfsum = gather(pssmcum, 'aapos', 'freq', -aa)
    gdfsum$aapos = as.factor(gdfsum$aapos)
    gdfsum$freq2 = gdf$freq
    # gdf = tibble(gdf)
    # gdfsum = tibble(gdfsum)
    print(gdfsum)
    print(gdf)
    ggplot(data=gdf) +
        geom_bar(mapping= aes(x=aapos,y=freq, fill=aa), stat='identity', position=position_fill(reverse=TRUE)) +
        geom_text(data=gdfsum, mapping=aes(x=aapos, y=freq, label=aa,size=freq2), nudge_y=-0.01) +
        labs(x='Amino acid position', y='Frequency', fill='Amino Acid') +
        scale_fill_manual(values=my_spectral) +
        guides(size=FALSE)
    outpng(infile, '')

}


blank_theme <- theme_minimal()+
  theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.border = element_blank(),
  panel.grid=element_blank(),
  axis.ticks = element_blank(),
  plot.title=element_text(size=14, face="bold")
  )

get_pie = function(infile){
    df = read_csv(infile)
    df = df[1:6000000,]
    binddf = df[df$pbind>0.5,]
    freqbind = dim(binddf)[1]/dim(df)[1]
    print(binddf)
    piedf = data.frame(seqtype = c('Binder (P > 0.5)', 'Non-binder (P \u2264 0.5)'), freq=c(freqbind, 1-freqbind))
    piedf$percent = sprintf('%s%%', round(piedf$freq,2)*100)
    print(piedf)
    ggplot(data=piedf) +
        geom_bar(mapping=aes(x='', y=freq, fill=seqtype), width=1, stat='identity',
        position=position_fill(reverse=TRUE)) +
        coord_polar('y', start=0) +
        labs(fill='') +
        blank_theme +
        theme(axis.text.x=element_blank(),
        legend.title = element_text(size=10),
        legend.text = element_text(size=20))+
        geom_text(aes(x='', y = cumsum(freq), label=percent), nudge_y=-0.07, size=10) +
        scale_fill_manual(values=c(nativecol, 'grey'))
    outpng(infile, 'pie')

}


get_position_freq2 = function(infile, seqvar){
    # infile = ('mason/outfiles/mHER_H3_AgPos_11300_11300_corpus_geminio_113000000_pbind.csv')
    df = read_csv(infile)
    # df = df[1:6000000,]
    aaseqs = df[[seqvar]]
    str_set = AAStringSet(aaseqs)
    pssm = consensusMatrix(str_set, as.prob=TRUE)
    pssmdf = as.data.frame(pssm)
    colnames(pssmdf) = sprintf('P%s', c(4:13))
    pssmcum = cumsum(pssmdf)
    pssmcum$aa = rownames(pssm)
    pssmdf$aa = rownames(pssm)
    gdf = gather(pssmdf, 'aapos', 'freq', -aa)
    gdf$aapos = factor(gdf$aapos, levels=mixedsort(unique(gdf$aapos)))
    gdfsum = gather(pssmcum, 'aapos', 'freq', -aa)
    gdfsum$aapos = as.factor(gdfsum$aapos)
    gdfsum$freq2 = gdf$freq
    # gdf = tibble(gdf)
    # gdfsum = tibble(gdfsum)
    print(gdfsum)
    print(gdf)
    ggplot(data=gdf) +
        geom_bar(mapping= aes(x=aapos,y=freq, fill=aa), stat='identity', position=position_fill(reverse=TRUE)) +
        geom_text(data=gdfsum, mapping=aes(x=aapos, y=freq, label=aa,size=freq2), nudge_y=-0.01) +
        labs(x='Amino acid position', y='Frequency', fill='Amino Acid') +
        scale_fill_manual(values=my_spectral) +
        guides(size=FALSE)
    outpng(infile, '')

}

get_mse = function(){
    genpos = read_csv('mason/outfiles/mHER_H3_AgPos_11300_11300_corpus_geminio_113000000_pbind.csv')[1:6000000,]
    genneg = read_csv('mason/outfiles/mHER_H3_AgNeg_27539_27539_corpus_geminio_275390000_pbind.csv')[1:6000000,]
    trainpos = read_csv('datasets/mason/mHER_H3_AgPos.csv')
    trainneg = read_csv('datasets/mason/mHER_H3_AgNeg.csv')
    print(trainneg)
    genaaseqpos = genpos[['gseq']]
    genaaseqneg = genneg[['gseq']]
    trainaaseqpos = trainpos[['AASeq']]
    trainaaseqneg = trainneg[['AASeq']]
    genstrpos_set = AAStringSet(genaaseqpos)
    genstrneg_set = AAStringSet(genaaseqneg)
    genpssmpos = consensusMatrix(genstrpos_set, as.prob=TRUE)
    genpssmneg = consensusMatrix(genstrneg_set, as.prob=TRUE)
    trainstrpos_set = AAStringSet(trainaaseqpos)
    trainstrneg_set = AAStringSet(trainaaseqneg)
    print(trainstrneg_set)
    trainpssmpos = consensusMatrix(trainstrpos_set, as.prob=TRUE)
    trainpssmneg = consensusMatrix(trainstrneg_set, as.prob=TRUE)
    gen_posneg = sum((genpssmpos-genpssmneg)^2)
    gentrain_pospos = sum((genpssmpos-trainpssmpos)^2)
    gentrain_negneg = sum((genpssmneg-trainpssmneg)^2)
    train_posneg = sum((trainpssmpos-trainpssmneg)^2)
    print(gen_posneg)
    print(train_posneg)
    print(gentrain_pospos)
    print(gentrain_negneg)
    stop()
}


get_distribution = function(infile){
    df = read_csv(infile)
    df = sample_n(df, 5000)
    print(df)
    median = median(df$pbind)
    ggplot(data=df) +
        geom_violin(mapping=aes(x='',y=pbind), outlier.shape = NA) +
        geom_jitter(mapping=aes(x='', y=pbind), color='#F5BF1A', alpha=0.07) +
        geom_text(mapping=aes(x='', y=median, label=sprintf('Median: %s',round(median,2))), nudge_y=0.01) +
        labs(x='Generated CDR-H3\nsequences', y= 'Binding probability')
    outpng(infile, 'proba_dist', 3,3)
}

# run stuff
# get_overlap()
get_overlap('mason/outfiles/mHER_H3_AgPos_11300_11300_corpus_geminio_113000000_pbind.csv')
get_overlap('mason/outfiles/mHER_H3_AgNeg_27539_27539_corpus_geminio_275390000_pbind.csv')
# get_position_freq('mason/outfiles/mHER_H3_AgPos_11300_11300_corpus_geminio_113000000_pbind.csv')
# get_position_freq('mason/outfiles/mHER_H3_AgNeg_27539_27539_corpus_geminio_275390000_pbind.csv')
# get_pie('mason/outfiles/mHER_H3_AgPos_11300_11300_corpus_geminio_113000000_pbind.csv')
# get_pie('mason/outfiles/mHER_H3_AgNeg_27539_27539_corpus_geminio_275390000_pbind.csv')
# get_position_freq2('datasets/mason/mHER_H3_AgPos.csv','AASeq')
# get_position_freq2('datasets/mason/mHER_H3_AgNeg.csv','AASeq')
# get_mse()
# get_distribution('mason/outfiles/mHER_H3_AgPos_11300_11300_corpus_geminio_113000000_pbind.csv')
# get_distribution('mason/outfiles/mHER_H3_AgNeg_27539_27539_corpus_geminio_275390000_pbind.csv')

