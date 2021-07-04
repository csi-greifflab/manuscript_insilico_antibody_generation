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

get_epitope_cor = function(){
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
                mdf = merge(fdf, sndf, by='AGbindPositions')
                mdf[is.na(mdf)] = 0
                print(mdf)
                print(dim(mdf))
                print('stuff')
                corval = cor(mdf$agpcount.x, mdf$agpcount.y)
                print(corval)
                cordf = data.frame(antigen = antigen,corval = corval, nsample=nsample, ntrainseq=ntrainseq)
                outdf = rbind(outdf, cordf)
            }
        }
    }
    outdf$ntrainseq = sprintf('ntrain%s', outdf$ntrainseq)
    outdf$ntrainseq = factor(outdf$ntrainseq, levels=c('ntrain700',
                                                        'ntrain7000',
                                                        'ntrain10000',
                                                        'ntrain20000',
                                                        'ntrain30000',
                                                        'ntrain40000',
                                                        'ntrain50000',
                                                        'ntrain60000',
                                                        'ntrain70000'))
    outname = 'eleven_outfiles/antigens_feature_native_generated_nsamples_count_cor.csv'
    write.csv(outdf, outname)
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
    outpng('fig', 'S5', width=15)
}

# run stuff
get_epitope_cor()
plot_epitope_cor()