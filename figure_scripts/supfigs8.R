# Title     : TODO
# Objective : TODO
# Created by: rahmadakbar
# Created on: 2020-10-07



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

get_epitope_count_distr = function(){
    ninfile = 'eleven_outfiles/antigens_feature_native_count.csv'
    ginfile= 'eleven_outfiles/antigens_feature_generated_nsamples_10k_count.csv'
    ndf = read_csv(ninfile)
    ndf = fixfive(ndf)
    gdf = read_csv(ginfile)
    gdf = fixfive(gdf)
    nantigens = unique(ndf$antigen)
    gantigens = unique(gdf$antigen)
    print(nantigens)
    print(gantigens)
    print(ndf)
    sagpcount = ndf %>% group_by(antigen) %>% summarize(sumval = sum(agpcount))
    print(sagpcount)
    print(gdf)
    gdf = gdf %>% group_by(antigen, nsample, ntrainseq) %>%
        mutate(agpcountfrac = agpcount/sum(agpcount))
    ndf = ndf %>% group_by(antigen, nsample, ntrainseq) %>%
        mutate(agpcountfrac = agpcount/sum(agpcount))
    print(gdf)
    print(ndf)
    mdf = rbind(gdf, ndf)
    mdf$ntrainseq = as.character(mdf$ntrainseq)
    print(unique(mdf$ntrainseq))
    # stop()
    mdf = mdf[mdf$antigen!='2DD8',]
    mdf = mdf[mdf$agpcountfrac>0.01,]
    mdf$ntrainseq[mdf$ntrainseq==0] = 'native'
    mdf$ntrainseq = factor(mdf$ntrainseq, levels=c(700,7000,10000,20000,
    30000,40000,50000,60000,70000,'native'))
    ggplot(data=mdf, mapping = aes(x=reorder(AGbindPositions,agpcountfrac), y=agpcountfrac)) +
        geom_boxplot(mapping=aes(fill=ntrainseq,color=ntrainseq), position = 'dodge2',outlier.shape=NA,lwd=1) +
        # geom_point(aes(color = origin), size = 5, shape = 21, position = position_jitterdodge()) +
        # geom_bar(mapping=aes(x=reorder(AGbindPositions,agpcountfrac), y=agpcountfrac,
        # fill=ntrainseq,color=ntrainseq), stat='identity') +
        facet_wrap(~antigen, nrow=2, scale='free') +
        theme(legend.title = element_blank())+
        theme(axis.text.x = element_blank())+
        scale_fill_manual(values=my_base_color2)+
        scale_color_manual(values=my_base_color2) +
        labs(x='Epitopes', y= 'Fractions of epitopes in models')
    outpng('fig', 'S9', width=20,height=10)
}

plot_epitope_cor = function(){
    infile = 'eleven_outfiles/antigens_feature_native_generated_nsamples_count_cor.csv'
    df = read_csv(infile)
    df = na.omit(df)
    df$ntrainseq = factor(df$ntrainseq, levels=c('ntrain700',
                                                        'ntrain7000',
                                                        'ntrain10000',
                                                        'ntrain20000',
                                                        'ntrain30000',
                                                        'ntrain40000',
                                                        'ntrain50000',
                                                        'ntrain60000',
                                                        'ntrain70000'))
    meddf = df %>% group_by(antigen, ntrainseq) %>% summarize(medval= median(corval))
    meddf$medval = round(meddf$medval,2)
    print(df)
    print(meddf)
    ggplot(data=df, mapping=aes(x=ntrainseq, y=corval)) +
        geom_boxplot(mapping=aes(color='black'), color=igemcol) +
        geom_point(mapping=aes(x=ntrainseq, y=corval, color=igemcol),size = 5, shape = 21, position =
        position_jitterdodge()) +
        geom_text(data= meddf, mapping=aes(x=ntrainseq, y=1.1, label=medval), size=3) +
        facet_wrap(~antigen, nrow=2) +
        labs(x='Models', y='Correlation (Pearson) between generated and native epitopes') +
        theme(axis.text.x = element_text(angle=90),
                legend.position='NA')
    outpng('fig', 'S6', width=15)
}



get_epitope_count_overlap = function(){
    ninfile = 'eleven_outfiles/antigens_feature_native_count.csv'
    ginfile= 'eleven_outfiles/antigens_feature_generated_nsamples_10k_count.csv'
    ndf = read_csv(ninfile)
    ndf = fixfive(ndf)
    gdf = read_csv(ginfile)
    gdf = fixfive(gdf)
    nantigens = unique(ndf$antigen)
    gantigens = unique(gdf$antigen)
    print(nantigens)
    print(gantigens)
    print(ndf)
    print(gdf)
    antigens = unique(ndf$antigen)
    print(antigens)
    antigens2 = antigens[antigens!='2DD8']
    print(antigens2)
    outdf = tibble()
    for (antigen in antigens) {
        sgdf = gdf[gdf$antigen==antigen,]
        sndf = ndf[ndf$antigen==antigen,]
        print(sgdf)
        print(sndf)
        nsamples = unique(sgdf$nsample)
        ntrainseqs = unique(sgdf$ntrainseq)
        print(nsamples)
        print(ntrainseqs)
        for (nsample in nsamples){
            for(ntrainseq in ntrainseqs){
                fdf = sgdf[(sgdf$nsample==nsample) & (sgdf$ntrainseq==ntrainseq),]
                print(fdf)
                print(sndf)
                mdf = merge(fdf, sndf, by='AGbindPositions', all=TRUE)
                # mdf = merge(fdf, sndf, by='AGbindPositions')
                mdf[is.na(mdf)] = 0
                print(mdf)
                # print(dim(mdf))
                # print(dim(fdf))
                # print(dim(sndf))
                # print('stuff')
                mdf$origin.x = 'generated'
                mdf$nsample.x = nsample
                mdf$ntrainseq.x = ntrainseq
                mdf$antigen.x = antigen
                mdf$origin.y = 'native'
                mdf$nsample.y = nsample
                mdf$ntrainseq.y = 0
                mdf$antigen.y = antigen
                print(mdf)
                # stop()
                # corval = cor(mdf$agpcount.x, mdf$agpcount.y)
                # print(corval)
                # cordf = data.frame(antigen = antigen,corval = corval, nsample=nsample, ntrainseq=ntrainseq)
                # outdf = rbind(outdf, cordf)
                outdf = rbind(outdf, mdf)
            }
        }
    }
    outdf$ntrainseq.x = sprintf('ntrain%s', outdf$ntrainseq.x)
    # outdf$ntrainseq = factor(outdf$ntrainseq, levels=c('ntrain700',
    # 'ntrain7000',
    # 'ntrain10000',
    # 'ntrain20000',
    # 'ntrain30000',
    # 'ntrain40000',
    # 'ntrain50000',
    # 'ntrain60000',
    # 'ntrain70000'))
    outname = 'eleven_outfiles/antigens_feature_native_generated_nsamples_count_xy.csv'
    write.csv(outdf, outname)
}


plot_epitope_countxy = function(){
    infile = 'eleven_outfiles/antigens_feature_native_generated_nsamples_count_xy.csv'
    df = read_csv(infile)
    df = na.omit(df)
    df$ntrainseq.x = factor(df$ntrainseq.x, levels=c('ntrain700',
    'ntrain7000',
    'ntrain10000',
    'ntrain20000',
    'ntrain30000',
    'ntrain40000',
    'ntrain50000',
    'ntrain60000',
    'ntrain70000'))
    # meddf = df %>% group_by(antigen.x, ntrainseq.x) %>% summarize(medval= median(corval))
    # meddf$medval = round(meddf$medval,2)
    # df = df[df$agpcount.x>500,]
    # stop()
    # print(meddf)
    ggplot(data=df, mapping=aes(x=reorder(AGbindPositions,agpcount.x), y=agpcount.x)) +
        geom_boxplot(mapping=aes(color=ntrainseq.x, fill=ntrainseq.x), position='dodge2', outlier.shape=NA) +
        # geom_point(mapping=aes(x=ntrainseq, y=corval, color=igemcol),size = 5, shape = 21, position =
        # position_jitterdodge()) +
        # geom_text(data= meddf, mapping=aes(x=ntrainseq, y=1.1, label=medval), size=3) +
        facet_wrap(~antigen.x, nrow=2, scale='free') +
        labs(x=' ', y='Concordance/overlap between generated and native epitopes') +
        theme(axis.text.x = element_blank(),
        legend.position='bottom',
        legend.title=element_blank()) +
        scale_color_manual(values=my_base_color2) +
        scale_fill_manual(values=my_base_color2)
    outpng('tmp', 'tmp', width=20, height=10)
}


plot_epitope_countxy2 = function(){
    infile = 'eleven_outfiles/antigens_feature_native_generated_nsamples_count_xy.csv'
    df = read_csv(infile)
    df = na.omit(df)
    df$ntrainseq.x = factor(df$ntrainseq.x, levels=c('ntrain700',
    'ntrain7000',
    'ntrain10000',
    'ntrain20000',
    'ntrain30000',
    'ntrain40000',
    'ntrain50000',
    'ntrain60000',
    'ntrain70000'))
    # meddf = df %>% group_by(antigen.x, ntrainseq.x) %>% summarize(medval= median(corval))
    # meddf$medval = round(meddf$medval,2)
    # df = df[df$agpcount.x>500,]
    # stop()
    # print(meddf)
    print(df)
    stop()
    ggplot(data=df, mapping=aes(x=reorder(AGbindPositions,agpcount.x), y=agpcount.x)) +
        geom_boxplot(mapping=aes(color=ntrainseq.x, fill=ntrainseq.x), position='dodge2', outlier.shape=NA) +
    # geom_point(mapping=aes(x=ntrainseq, y=corval, color=igemcol),size = 5, shape = 21, position =
    # position_jitterdodge()) +
    # geom_text(data= meddf, mapping=aes(x=ntrainseq, y=1.1, label=medval), size=3) +
        facet_wrap(~antigen.x, nrow=2, scale='free') +
        labs(x=' ', y='Concordance/overlap between generated and native epitopes') +
        theme(axis.text.x = element_blank(),
        legend.position='bottom',
        legend.title=element_blank()) +
        scale_color_manual(values=my_base_color2) +
        scale_fill_manual(values=my_base_color2)
    outpng('tmp', 'tmp', width=20, height=10)
}


get_epitope_count_overlap2 = function(){
    ninfile = 'eleven_outfiles/antigens_feature_native_count.csv'
    ginfile= 'eleven_outfiles/antigens_feature_generated_nsamples_10k_count.csv'
    ndf = read_csv(ninfile)
    ndf = fixfive(ndf)
    gdf = read_csv(ginfile)
    gdf = fixfive(gdf)
    nantigens = unique(ndf$antigen)
    gantigens = unique(gdf$antigen)
    print(nantigens)
    print(gantigens)
    print(ndf)
    print(gdf)
    antigens = unique(ndf$antigen)
    print(antigens)
    antigens2 = antigens[antigens!='2DD8']
    print(antigens2)
    outdf = tibble()
    for (antigen in antigens) {
        sgdf = gdf[gdf$antigen==antigen,]
        sndf = ndf[ndf$antigen==antigen,]
        print(sgdf)
        print(sndf)
        nsamples = unique(sgdf$nsample)
        ntrainseqs = unique(sgdf$ntrainseq)
        print(nsamples)
        print(ntrainseqs)
        for (nsample in nsamples){
            for(ntrainseq in ntrainseqs){
                fdf = sgdf[(sgdf$nsample==nsample) & (sgdf$ntrainseq==ntrainseq),]
                print(fdf)
                print(sndf)
                mdf = merge(fdf, sndf, by='AGbindPositions', all=TRUE)
                # mdf = merge(fdf, sndf, by='AGbindPositions')
                mdf$agpcount.x[is.na(mdf$agpcount.x)] = 0
                print(mdf)
                nepitope = dim(sndf)[1]
                gepitope = dim(fdf)[1]
                noverlap = length(intersect(sndf$AGbindPositions, fdf$AGbindPositions))
                mdf$origin.x = 'generated'
                mdf$nsample.x = nsample
                mdf$ntrainseq.x = ntrainseq
                mdf$antigen.x = antigen
                mdf$origin.y = 'native'
                mdf$nsample.y = nsample
                mdf$ntrainseq.y = 0
                mdf$antigen.y = antigen
                print(mdf)
                overlapdf = data.frame(antigen = antigen,noverlap = noverlap, nepitope=nepitope,
                gepitope=gepitope, nsample=nsample, ntrainseq=ntrainseq)
                if (nsample!='sample0' & ntrainseq==70000){
                    print('passing this one')
                }
                else{outdf = rbind(outdf, overlapdf)}
            }
        }
    }
    outdf$ntrainseq = sprintf('ntrain%s', outdf$ntrainseq)
    # outdf$ntrainseq = factor(outdf$ntrainseq, levels=c('ntrain700',
    # 'ntrain7000',
    # 'ntrain10000',
    # 'ntrain20000',
    # 'ntrain30000',
    # 'ntrain40000',
    # 'ntrain50000',
    # 'ntrain60000',
    # 'ntrain70000'))
    print(outdf)
    outname = 'eleven_outfiles/antigens_feature_native_generated_nsamples_count_overlap.csv'
    write.csv(outdf, outname, row.names=FALSE)
}

plot_epitope_count_overlap2 = function(){
    infile = 'eleven_outfiles/antigens_feature_native_generated_nsamples_count_overlap.csv'
    df = read_csv(infile)
    df = na.omit(df)
    df$ntrainseq = factor(df$ntrainseq, levels=c('ntrain700',
    'ntrain7000',
    'ntrain10000',
    'ntrain20000',
    'ntrain30000',
    'ntrain40000',
    'ntrain50000',
    'ntrain60000',
    'ntrain70000'))
    # meddf = df %>% group_by(antigen.x, ntrainseq.x) %>% summarize(medval= median(corval))
    # meddf$medval = round(meddf$medval,2)
    # df = df[df$agpcount.x>500,]
    # stop()
    # print(meddf)
    print(df)
    df = gather(df, key='counttype', value='countval', noverlap, nepitope, gepitope)
    print(df)
    df$counttype = factor(df$counttype, levels=c('nepitope', 'gepitope', 'noverlap'))
    meddf = df %>% group_by(antigen, ntrainseq, counttype) %>% summarize(medval=median(countval))
    ggplot(data=df, mapping=aes(x=ntrainseq, y=countval)) +
        geom_boxplot(mapping=aes(color=counttype, fill=counttype), position='dodge2', outlier.shape=NA) +
    # geom_point(mapping=aes(color=counttype),size = 5, shape = 21, position =
    # position_jitterdodge()) +
    geom_text(data= meddf, mapping=aes(x=ntrainseq, y=1.1, label=medval, color=counttype), size=3,
        position=position_dodge2(width=1)) +
        facet_wrap(~antigen, nrow=2) +
        labs(x='Number of training CDR-H3 sequences',
        y='Overlap of epitopes recognized by\nnative and generated CDR-H3 sequences') +
        theme(axis.text.x = element_text(angle=90),
        legend.position='bottom',
        legend.title=element_blank(),
        axis.title=element_text(size=25),
        legend.text=element_text(size=25)) +
        scale_color_manual(values=c(nativecol, igemcol, 'black')) +
        scale_fill_manual(values=c(nativecol, igemcol, 'black'))
    outpng('fig', 'S9b', width=22, height=10)
}

# run stuff
# get_epitope_count_distr()
# plot_epitope_cor()
# get_epitope_count_overlap()
# plot_epitope_countxy()
# plot_epitope_countxy2()
get_epitope_count_overlap2()
plot_epitope_count_overlap2()
