# Title     : TODO
# Objective : TODO
# Created by: rahmadakbar
# Created on: 2020-09-15


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





fixfive2 =function(df){
    # fix 5E94
    df$antigen[df$antigen=='FiveE94'] = '5E94'
    df$antigen2[df$antigen2=='FiveE94'] = '5E94'
    return(df)
}

fixfive3 =function(df){
    # fix 5E94
    df$antigen[df$antigen=='FiveE94'] = '5E94'
    df$antigen2[df$crossantigen=='FiveE94'] = '5E94'
    return(df)
}

 freq2med = function(ld, freq){
        d2 = rep(ld, freq)
        med = median(d2)
    }

ld_native_generated_distr = function(){
    # infile = 'eleven_outfiles/antigens_native_generated_ldcount_merged_ldmerged.csv'
    infile = 'eleven_outfiles/antigens_native_generated_ldcount_rsample_merged.csv'
    df = read_csv(infile)
    df = fixfive2(df)
    df = df[df$origin!='generated_native',]
    df = df[df$antigen==df$antigen2,]
    print(df)
    df$origin = factor(df$origin, levels=c('native_native', 'generated_generated'))
    # meddf = df %>%  group_by(antigen, antigen2, origin) %>% summarize(medval=freq2med(ld, count))
    meadf = df %>% group_by(antigen, antigen2, ld,origin) %>% summarize(mean_nld=mean(nld),
        max_nld = max(nld),sd_nld=sd(nld), mean_ld=mean(ld))
    meddf = meadf %>% group_by(antigen, antigen2,origin) %>% summarize(medval=freq2med(ld,mean_nld),
        maxval=max(mean_nld))
    print(meddf)
    print(meadf)
    meddf$medlabel  = sprintf('Median: %s', meddf$medval)
    nmeddf = meddf[meddf$origin=='native_native',]
    gmeddf = meddf[meddf$origin!='native_native',]
    ggplot(data=meadf) +
        geom_bar(mapping=aes(x=ld, y=mean_nld, fill=origin), stat='identity', position='dodge2') +
        geom_text(data=nmeddf, mapping=aes(x=30, y=150000, label=medlabel), color=nativecol, hjust=1) +
        geom_errorbar(mapping=aes(x=ld, y=mean_nld,ymin=mean_nld-sd_nld,ymax=mean_nld+sd_nld),
            position=position_dodge2(width=1), size=0.2) +
        geom_text(data=gmeddf, mapping=aes(x=30, y=125000, label=medlabel), color=igemcol, hjust=1) +
        facet_wrap(~antigen, scale='free', nrow=2) +
        labs(x='Levenshtein distance (LD)', y= '# of CDR-H3 sequences') +
        theme(legend.title=element_blank(),
                legend.text = element_text(size=15)) +
        scale_fill_manual(values=c(nativecol, igemcol))

    outpng('fig', 'S1', width=17)

}

# runs stuff
ld_native_generated_distr()
