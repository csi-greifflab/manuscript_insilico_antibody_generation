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
    tinfile = 'eleven_outfiles/antigens_feature_crosstransfer_generated_nsamples_10k_count.csv'
    ndf = read_csv(ninfile)
    ndf = fixfive(ndf)
    gdf = read_csv(ginfile)
    gdf = fixfive(gdf)
    gdf = gdf[(gdf$ntrainseq==700) | (gdf$ntrainseq==7000),]
    tdf = read_csv(tinfile)
    tdf = fixfive(tdf)
    nantigens = unique(ndf$antigen)
    gantigens = unique(gdf$antigen)
    gdf = gdf %>% group_by(antigen, nsample, ntrainseq) %>%
        mutate(agpcountfrac = agpcount/sum(agpcount))
    ndf = ndf %>% group_by(antigen, nsample, ntrainseq) %>%
        mutate(agpcountfrac = agpcount/sum(agpcount))
    tdf = tdf %>% group_by(antigen, nsample, ntrainseq) %>%
        mutate(agpcountfrac=agpcount/sum(agpcount))
    ndf$crossstatus = rep('-T', dim(ndf)[1])
    gdf$crossstatus = rep('-T', dim(gdf)[1])
    tdf$crossstatus = rep('+T', dim(tdf)[1])
    ndf$model_name = sprintf('%s_ntrain%s', ndf$crossstatus, ndf$ntrainseq)
    gdf$model_name = sprintf('%s_ntrain%s', gdf$crossstatus, gdf$ntrainseq)
    tdf$model_name = sprintf('%s_ntrain%s', tdf$crossstatus, tdf$ntrainseq)
    print(ndf)
    print(gdf)
    print(tdf)
    mdf = rbind(ndf,gdf,tdf)
    mdf = mdf[mdf$agpcountfrac>0.01,]
    print(unique(mdf$model_name))
    mdf$model_name[mdf$model_name=='-T_ntrain0'] = 'native'
    mdf = mdf[mdf$antigen != '2DD8',]
    mdf$model_name = factor(mdf$model_name, levels=c('-T_ntrain700',
                                                        '+T_ntrain700',
                                                        '-T_ntrain7000',
                                                        '+T_ntrain7000',
                                                        'native'))
    ggplot(data=mdf) +
        geom_boxplot(mapping=aes(x=reorder(AGbindPositions,agpcountfrac), y=agpcountfrac, fill=model_name,
        color=model_name), position='dodge2', outlier.shape=NA, lwd=1) +
        facet_wrap(~antigen, scale='free', nrow=2) +
        theme(axis.text.x=element_blank()) +
        scale_fill_manual(values=my_base_color2)+
        scale_color_manual(values=my_base_color2)+
        labs(x='Epitopes', y= 'Fractions of epitopes in models')
    outpng('fig', 'S11',widt=20, height=10)
}

plot_epitope_cor = function(){
    infile = 'eleven_outfiles/antigens_feature_native_transfer_generated_nsamples_count_cor.csv'
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
    df$xtag = sprintf('%s_%s', df$crossstatus, df$ntrainseq)
    df$xtag = factor(df$xtag,levels=c('-T_ntrain700',
        '+T_ntrain700',
        '-T_ntrain7000',
        '+T_ntrain7000'))
    meddf = df %>% group_by(antigen, ntrainseq,xtag) %>% summarize(medval= median(corval))
    meddf$medval = round(meddf$medval,2)
    print(df)
    print(meddf)
    ggplot(data=df, mapping=aes(x=xtag, y=corval)) +
        geom_boxplot(mapping=aes(color=crossstatus)) +
        geom_point(mapping=aes(color=crossstatus), size = 5, shape = 21, position = position_jitterdodge()) +
        geom_text(data= meddf, mapping=aes(x=xtag, y=1.1, label=medval), size=3) +
        facet_wrap(~antigen, nrow=2) +
        labs(x='Models', y='Correlation (Pearson) between generated and native epitopes') +
        theme(axis.text.x = element_text(angle=90),
                legend.position='NA') +
        scale_color_manual(values=c(my_base_color[2], 'black'))
    outpng('fig', 'S7', width=15)
}


get_epitope_count_overlap = function(){
    ninfile = 'eleven_outfiles/antigens_feature_native_count.csv'
    ginfile= 'eleven_outfiles/antigens_feature_generated_nsamples_10k_count.csv'
    tinfile = 'eleven_outfiles/antigens_feature_crosstransfer_generated_nsamples_10k_count.csv'
    ndf = read_csv(ninfile)
    ndf = fixfive(ndf)
    gdf = read_csv(ginfile)
    gdf = fixfive(gdf)
    gdf = gdf[(gdf$ntrainseq==700) | (gdf$ntrainseq==7000),]
    tdf = read_csv(tinfile)
    tdf = fixfive(tdf)
    nantigens = unique(ndf$antigen)
    nantigens = nantigens[nantigens!='2DD8']
    gantigens = unique(gdf$antigen)
    print(ndf)
    print(gdf)
    print(tdf)
    nsamples = unique(gdf$nsample)
    print(nsamples)
    ntrainseqs = unique(gdf$ntrainseq)
    print(ntrainseqs)
    outdf = tibble()
    for (antigen in nantigens){
        for (nsample in nsamples){
            for (ntrainseq in ntrainseqs){
                sndf = ndf[ndf$antigen==antigen,]
                sgdf = gdf[(gdf$antigen==antigen) & (gdf$nsample==nsample) & (gdf$ntrainseq==ntrainseq),]
                stdf = tdf[(tdf$antigen==antigen) & (tdf$nsample==nsample) & (tdf$ntrainseq==ntrainseq),]
                print(sndf)
                print(sgdf)
                print(stdf)
                nepitope = dim(sndf)[1]
                gepitope = dim(sgdf)[1]
                tepitope = dim(stdf)[1]
                gnoverlap = length(intersect(sndf$AGbindPositions, sgdf$AGbindPositions))
                tnoverlap = length(intersect(sndf$AGbindPositions, stdf$AGbindPositions))
                print(nepitope)
                print(gepitope)
                print(tepitope)
                print(gnoverlap)
                print(tnoverlap)
                overlapdf = data.frame(antigen = antigen,gnoverlap = gnoverlap, tnoverlap=tnoverlap,
                nepitope=nepitope, gepitope=gepitope, tepitope=tepitope, nsample=nsample, ntrainseq=ntrainseq)
                outdf=rbind(outdf, overlapdf)}
        }
    }
    print(outdf)
    outname = 'eleven_outfiles/antigens_feature_native_crosstransfer_generated_nsamples_count_overlap.csv'
    write.csv(outdf, outname, row.names=FALSE)
}


plot_epitope_count_overlap = function() {
    infile = 'eleven_outfiles/antigens_feature_native_crosstransfer_generated_nsamples_count_overlap.csv'
    df = read_csv(infile)
    print(unique(df$ntrainseq))
    df$ntrainseq=sprintf('ntrain%s', df$ntrainseq)
    print(unique(df$ntrainseq))
    df$ntrainseq=factor(df$ntrainseq, levels=c('ntrain700', 'ntrain7000'))
    print(df)
    gdf = gather(df,key='counttype', value='countvalue', gnoverlap, tnoverlap, nepitope, gepitope, tepitope)
    print(gdf)
    gdf$ntrainseq = factor(gdf$ntrainseq, levels=c('ntrain700',
    'ntrain7000',
    'ntrain10000',
    'ntrain20000',
    'ntrain30000',
    'ntrain40000',
    'ntrain50000',
    'ntrain60000',
    'ntrain70000'))
    print(gdf)
    gdf$counttype[gdf$counttype=='gepitope'] = '-T_epitope'
    gdf$counttype[gdf$counttype=='tepitope'] = '+T_epitope'
    gdf$counttype[gdf$counttype=='gnoverlap'] = '-T_overlap'
    gdf$counttype[gdf$counttype=='tnoverlap'] = '+T_overlap'
    gdf$counttype[gdf$counttype=='nepitope'] = 'native_epitope'
    print(gdf)
    gdf$counttype = factor(gdf$counttype, levels=c('-T_epitope', '+T_epitope', 'native_epitope','-T_overlap',
    '+T_overlap'))
    meddf = gdf %>% group_by(antigen, ntrainseq, counttype) %>% summarize(medval = median(countvalue))
    print(meddf)
    ggplot(data=gdf, mapping=aes(x=ntrainseq, y=countvalue)) +
        geom_boxplot(mapping=aes(color=counttype, fill=counttype),
        position='dodge2', outlier.shape=NA) +
        geom_text(data=meddf, mapping=aes(x=ntrainseq, y=-10, color=counttype,label=medval),
        position=position_dodge2(width=1)) +
    # geom_point(mapping=aes(color=counttype),size = 3, shape = 21, position =
    #     position_jitterdodge()) +
        facet_wrap(~antigen, nrow=2) +
        scale_color_manual(values=c('black', 'orange', nativecol, 'black', 'orange'))+
        scale_fill_manual(values=c('black', 'orange', nativecol, 'black', 'orange'))+
        labs(x='Number of training CDR-H3 sequences',
        y='Overlap of epitopes recognized by\nnative and generated CDR-H3 sequences') +
        theme(legend.title=element_blank())
    outpng('fig', 'S11b', width=15)

}

# run stuff
# get_epitope_count_distr()
# plot_epitope_cor()
get_epitope_count_overlap()
plot_epitope_count_overlap()
