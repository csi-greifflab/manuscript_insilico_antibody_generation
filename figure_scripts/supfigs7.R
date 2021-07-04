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
    tinfile = 'eleven_outfiles/antigens_feature_crosstransfer_generated_nsamples_10k_count.csv'
    ndf = read_csv(ninfile)
    ndf = fixfive(ndf)
    gdf = read_csv(ginfile)
    gdf = fixfive(gdf)
    tdf = read_csv(tinfile)
    tdf = fixfive(tdf)
    nantigens = unique(ndf$antigen)
    gantigens = unique(gdf$antigen)
    print(nantigens)
    print(gantigens)
    print(ndf)
    print(gdf)
    print(tdf)
    antigens = unique(ndf$antigen)
    print(antigens)
    antigens2 = antigens[antigens!='2DD8']
    print(antigens2)
    crossantigens = unique(tdf$crossantigen)
    print(crossantigens)
    outdf = tibble()
    for (antigen in antigens) {
        sgdf = gdf[gdf$antigen==antigen,]
        sndf = ndf[ndf$antigen==antigen,]
        stdf = tdf[tdf$antigen==antigen,]
        print(sgdf)
        print(sndf)
        nsamples = unique(sgdf$nsample)
        ntrainseqs = unique(stdf$ntrainseq)
        print(nsamples)
        print(ntrainseqs)
        for (nsample in nsamples){
            for(ntrainseq in ntrainseqs){
                fgdf = sgdf[(sgdf$nsample==nsample) & (sgdf$ntrainseq==ntrainseq),]
                ftdf = stdf[(stdf$nsample==nsample) & (stdf$ntrainseq==ntrainseq),]
                print(fgdf)
                print(ftdf)
                print(sndf)
                mgdf = merge(fgdf, sndf, by='AGbindPositions')
                mgdf[is.na(mgdf)] = 0
                gcorval = cor(mgdf$agpcount.x, mgdf$agpcount.y)
                gcordf = data.frame(antigen = antigen, corval = gcorval, nsample=nsample, ntrainseq=ntrainseq,
                crossstatus='-T', crossantigen=antigen)
                outdf = rbind(outdf, gcordf)
                for (crossantigen in crossantigens){
                    ftdfcross = ftdf[ftdf$crossantigen==crossantigen,]
                    mtdf = merge(ftdfcross, sndf, by='AGbindPositions')
                    mtdf[is.na(mtdf)] = 0
                    print(mtdf)
                    print(dim(mtdf))
                    print(mgdf)
                    print(dim(mgdf))
                    print('stuff')
                    tcorval = cor(mtdf$agpcount.x, mtdf$agpcount.y)
                    print(tcorval)
                    print(gcorval)
                    tcordf = data.frame(antigen = antigen, corval = tcorval, nsample=nsample, ntrainseq=ntrainseq,
                    crossstatus='+T', crossantigen=crossantigen)
                    outdf = rbind(outdf, tcordf)

                }
            }
        }
    }
    outdf = outdf %>% drop_na()
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
    outname = 'eleven_outfiles/antigens_feature_native_crosstransfer_generated_nsamples_count_cor.csv'
    write.csv(outdf, outname, row.names=FALSE)
}

plot_epitope_cor = function(){
    infile = 'eleven_outfiles/antigens_feature_native_crosstransfer_generated_nsamples_count_cor.csv'
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
        labs(x='Models', y='Correlation of epitope occupancy of native\nand generated CDR-H3 sequences') +
        theme(axis.text.x = element_text(angle=90),
                legend.position='NA') +
        scale_color_manual(values=c('orange', 'black'))
    outpng('fig', 'S7', width=15)
}

# run stuff
get_epitope_cor()
plot_epitope_cor()