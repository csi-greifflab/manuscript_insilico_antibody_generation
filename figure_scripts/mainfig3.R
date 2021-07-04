# Title     : TODO
# Objective : TODO
# Created by: rahmadakbar
# Created on: 2020-08-21



# viz plots in main figure 3

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

my_gradient= colorRampPalette(c('red', 'black'))(9)
print(my_gradient)

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
    df$antigen.x[df$antigen.x=='FiveE94'] = '5E94'
    df$antigen.y[df$antigen.y=='FiveE94'] = '5E94'
    return(df)
}

viz_binder_bins = function(){
    infile  = 'eleven_outfiles/antigens_affinity_binders_bin.csv'
    df = read_csv(infile)
    # df$xtag = sprintf('%s_%s_%s')
    print(df)
    # df = df[df$antigen=='1ADQ',]
    df$binder_class = factor(df$binder_class, levels=c('hyper_binder', 'ultimate_binder', 'penultimate_binder',
    'binder'))
    gdf = df[df$origin=='generated',]
    ndf = df[df$origin=='native',]
    ggplot(data=gdf) +
    # geom_line(mapping=aes(x=ngenerated, y=nbinder_class, color=origin), size=1.4) +
        geom_path(mapping=aes(x=ngenerated, y=nbinder_class, color=origin, shape=binder_class), size=0.5) +
        geom_line(data=ndf, mapping=aes(x=ngenerated, y=nbinder_class, color=origin, shape=binder_class,
        linetype='dashed'), size=1.5) +
        geom_point(mapping=aes(x=ngenerated, y=nbinder_class, color=origin, shape=binder_class), size=5) +
        facet_wrap(~antigen) +
        scale_y_log10() +
        scale_x_log10() +
        labs(x='N generated', y=' ') +
        theme(legend.title = element_blank(),
        legend.position='none',
        axis.text.x = element_text(size=20))+
        scale_color_manual(values=c(igemcol, nativecol))

    outpng(infile, 'linepath', width=20, height=10)
}


viz_binder_bins2 = function(){
    infile  = 'eleven_outfiles/antigens_affinity_binders_bin_permuted.csv'
    df = read_csv(infile)
    # print(df)
    df$binder_class = factor(df$binder_class, levels=c('hyper_binder', 'ultimate_binder', 'penultimate_binder',
    'binder'))
    gdf = df[df$origin=='generated',]
    gdf=fixfive(gdf)
    medgdf = gdf %>%
        group_by(antigen, binder_class, ngenerated, origin) %>%
        summarize('mean' = mean(nbinder_class), 'std' = sd(nbinder_class))
    print(medgdf)
    ndf = df[(df$origin=='native') & (df$npermutation==1),]
    ndf=fixfive(ndf)
    myxticks = unique(medgdf$ngenerated)
    myyticks = unique(gdf$nbinder_class)
    myyticks = c(10, 100,10000,100000,400000)
    print(myxticks)
    print(myyticks)
    ggplot(data=medgdf) +
        geom_path(mapping=aes(x=ngenerated, y=mean, color=origin, shape=binder_class), size=0.5) +
        geom_line(data=ndf, mapping=aes(x=ngenerated, y=nbinder_class, color=origin, shape=binder_class,
        linetype='dashed'), size=1.5) +
        geom_point(mapping=aes(x=ngenerated, y=mean, color=origin, shape=binder_class), size=7) +
        geom_errorbar(mapping=aes(x=ngenerated, ymin=mean-std, ymax=mean+std),color='grey',size=2) +
        scale_x_log10(
        breaks = myxticks,
        labels = myxticks
        ) +
        scale_y_log10(
        breaks = myyticks,
        labels = myyticks
        ) +
        facet_wrap(~antigen, nrow=2) +
        labs(x='# of generated CDR-H3 sequences (log10)', y='# of binders per binder type (log10)') +
        theme(axis.text.y = element_text(size=20),
        axis.ticks.y=element_blank(),
        axis.text.x = element_text(angle=90, size= 13),
        legend.position = 'None',
            axis.title = element_text(size=30),
            strip.text = element_text(size=25)
        ) +
        scale_color_manual(values=c(igemcol, nativecol))

    outpng(infile, 'linepath', width=20, height=10)
}


bin_bar2 = function(infile){
    # infile='eleven_outfiles/antigens_affinity_native_generated_take2_devel_bin_small.csv'
    df = read_csv(infile)
    print(df)
    df$mtag = sprintf('%s_%s_%s_%s_%s_%s_%s', df$charge_at_7_binned, df$mol_weight_binned, df$gravy_binned,
    df$instability_index_binned, df$length_binned, df$mhc2_binned, df$mhc_binned)
    cdf = df %>% count(mtag, antigen, origin)
    cdf$xtag= sprintf('%s_%s', cdf$antigen, cdf$origin)
    cdf$origin = factor(cdf$origin,levels=c('native', 'generated'))
    cdf$xcolor = ifelse(cdf$origin=='generated', igemcol,nativecol)
    cdf = cdf %>% group_by(xtag) %>% mutate(xtagcount = n())
    cdf = cdf %>% group_by(mtag) %>% mutate(mtagcount = n())
    cdf = cdf %>% group_by(mtag) %>% mutate(mtagsum = sum(n))
    cdf$mtag_fct = factor(cdf$mtag)
    cdf$mtag_fct=fct_reorder2(cdf$mtag_fct, cdf$mtagcount, cdf$mtagsum)
    print(cdf)
    ggplot(cdf) +
        geom_tile(mapping=aes(x=xtag, y= mtag_fct, fill=n, color=origin, height=0.6, width=0.6),
        size=0.6) +
    #     geom_point(mapping=aes(x=xtag, y=mtag, size=n, color=origin)) +
        theme(axis.text.x = element_text(angle=90, size=15, color=cdf$xcolor),
        legend.position='bottom',
        axis.title.y = element_text(size=20),
        legend.text=element_text(size=20),
        legend.title = element_blank(),
        legend.key.width = unit(2.5, 'cm')) +
        scale_color_manual(values=c(nativecol, igemcol)) +
        scale_fill_gradient(low='white', high='red') +
        labs(x=' ', y='Developability encoding')
    outpdf(infile, 'developability_encoding', width=12, height=15)
}


bin_bar_large_generation = function(){
    # infile='eleven_outfiles/antigens_affinity_generated_native_develbin_sampled.csv'
    infile='eleven_outfiles/antigens_affinity_generated_native_develbin.csv'
    infile_small='eleven_outfiles/antigens_affinity_native_generated_take2_devel_bin_small.csv'
    dfsmall = read_csv(infile_small)
    dfsmall = dfsmall[dfsmall$origin=='generated',]
    dfsmall$origin = sprintf('%s_small', dfsmall$origin)
    dfsmall$sampletag = rep('small', dim(dfsmall)[1])
    df = read_csv(infile)
    # df$sampletag = ifelse(df$origin=='generated', 'large', 'ntrain')
    df$sampletag = ifelse(df$origin=='generated', 'large', 'ntrain')
    df$origin = sprintf('%s_%s', df$origin, df$sampletag)
    df = rbind(df, dfsmall)
    df$origin = factor(df$origin,levels=c('native_ntrain', 'generated_large', 'generated_small'))
    # df$mtag = sprintf('%s_%s_%s_%s_%s_%s_%s', df$charge_at_7_binned, df$mol_weight_binned, df$gravy_binned,
    # df$instability_index_binned, df$length_binned, df$mhc2_binned, df$mhc_binned)
    df$mtag = sprintf('%s_%s_%s_%s_%s_%s', df$charge_at_7_binned, df$mol_weight_binned, df$gravy_binned,
    df$instability_index_binned, df$length_binned, df$mhc2_binned, df$mhc_binned)
    print(unique(df$origin))
    origin_colors = c('native_ntrain'=nativecol, 'generated_large' = igemcol_large, 'generated_small'=igemcol)
    print(origin_colors)
    df$xtag= sprintf('%s_%s', df$antigen, df$origin)
    print(unique(df$xtag))
    bardf = df %>% count(xtag, mtag, origin) %>% count(xtag,origin)
    print(bardf)
    ## barplot
    ggplot(data=bardf) +
        geom_bar(mapping=aes(x=xtag, y=n, fill= origin), stat='identity')+
        geom_text(mapping=aes(x=xtag, y=n, label=n), nudge_y=2, size=7) +
        scale_fill_manual(values=c(nativecol,igemcol_large, igemcol)) +
        labs(y='# of developability combinations') +
        theme(axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        legend.position='NA',
        axis.text.y= element_text(size=30),
        axis.title.y=element_text(size=30))
    outpng('antigen', 'developability_encoding_large_bar', width=30)
    ## end barplot
    # stop()
    cdf = df %>% count(xtag, mtag, antigen, origin, sampletag)
    # cdf$xcolor = ifelse(cdf$origin=='generated', igemcol,nativecol)
    cdf$xcolor = sapply(cdf$origin, function(x) origin_colors[x])
    cdf = cdf %>% group_by(xtag) %>% mutate(xtagcount = n())
    cdf = cdf %>% group_by(mtag) %>% mutate(mtagcount = n())
    cdf = cdf %>% group_by(mtag) %>% mutate(mtagsum = sum(n))
    cdf$mtag_fct = factor(cdf$mtag)
    cdf$mtag_fct=fct_reorder2(cdf$mtag_fct, cdf$mtagcount, cdf$mtagsum)
    get_devel_sum = function(x) {
        sumbin = sum(as.numeric(unlist(strsplit(x, split='_'))))
        sumbin = sumbin+1
        return(sumbin)
    }
    cdf$sumbin = sapply(cdf$mtag, function(x) get_devel_sum(x))
    umtag = cdf$mtag_fct
    umtag2 = cdf$mtag
    print(umtag)
    print(umtag2)
    umtaglevel = levels(umtag)
    umtagcol = sapply(umtaglevel, function(x) get_devel_sum(x))
    ggplot(data=cdf) +
        geom_tile(mapping=aes(x=xtag, y= mtag_fct, fill=n, color=origin, height=0.6, width=0.6),
        size=0.7) +
    #     geom_point(mapping=aes(x=xtag, y=mtag, size=n, color=origin)) +
        theme(axis.text.x = element_text(angle=90, size=25, color=cdf$xcolor),
        legend.position='NA',
        axis.title.y = element_text(size=30),
        axis.text.y = element_text(size=15,color=umtagcol, face='bold'),
        legend.text=element_text(size=20),
        legend.title = element_blank(),
        legend.key.width = unit(2.5, 'cm')) +
        scale_color_manual(values=c(nativecol, igemcol_large, igemcol)) +
        scale_fill_gradient(low='white', high='red') +
        labs(x=' ', y='Developability encoding')
    outpng('antigen', 'developability_encoding_large', width=30, height=15)
}


bin_bar_sln = function(){
    infile = 'eleven_outfiles/antigens_affinity_generated_devel_native_develbin_sln.csv'
    df = read_csv(infile)
    df$mtag = sprintf('%s_%s_%s_%s_%s_%s', df$charge_at_7_binned, df$mol_weight_binned, df$gravy_binned,
    df$instability_index_binned, df$length_binned, df$mhc2_binned, df$mhc_binned)
    df$xtag= sprintf('%s_%s', df$antigen, df$sample_tag)
    print(unique(df$xtag))
    print(unique(df$mtag))
    bardf= df %>% count(xtag,mtag,sample_tag) %>%count(xtag,sample_tag)
    ggplot(data=bardf) +
        geom_bar(mapping=aes(x=xtag, y=n, fill= sample_tag), stat='identity')+
        geom_text(mapping=aes(x=xtag, y=n, label=n), nudge_y=2, size=7) +
        scale_fill_manual(values=c(igemcol_large, igemcol,nativecol)) +
        labs(y='# of developability combinations') +
        theme(axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        legend.position='NA',
        axis.text.y= element_text(size=30),
        axis.title.y=element_text(size=30))
    outpng('antigen', 'developability_encoding_sln_bar', width=30)


}


bin_tile_sln = function(){
    infile = 'eleven_outfiles/antigens_affinity_generated_devel_native_develbin_sln.csv'
    df = read_csv(infile)
    df$mtag = sprintf('%s_%s_%s_%s_%s_%s', df$charge_at_7_binned, df$mol_weight_binned, df$gravy_binned,
    df$instability_index_binned, df$length_binned, df$mhc2_binned, df$mhc_binned)
    df = fixfive(df)
    df$xtag= sprintf('%s_%s', df$antigen, df$sample_tag)
    print(unique(df$xtag))
    print(unique(df$mtag))
    tdf= df %>% count(xtag,mtag,sample_tag)
    tdf = tdf %>% group_by(xtag) %>% mutate(xtagcount = n())
    tdf = tdf %>% group_by(mtag) %>% mutate(mtagcount = n())
    tdf = tdf %>% group_by(mtag) %>% mutate(mtagsum = sum(n))
    tdf$mtag_fct = factor(tdf$mtag)
    tdf$mtag_fct=fct_reorder2(tdf$mtag_fct, tdf$mtagcount, tdf$mtagsum)
    get_devel_sum = function(x) {
        sumbin = sum(as.numeric(unlist(strsplit(x, split='_'))))
        sumbin = sumbin+1
        return(sumbin)
    }
    get_devel_col = function(x) {
        sumbin = sum(as.numeric(unlist(strsplit(x, split='_'))))
        sumbin = sumbin+1
        return(my_gradient[sumbin])
    }
    tdf$sumbin = sapply(tdf$mtag, function(x) get_devel_sum(x))
    umtag = tdf$mtag_fct
    umtag2 = tdf$mtag
    umtaglevel = levels(umtag)
    umtagcol = sapply(umtaglevel, function(x) get_devel_col(x))
    umtagsum = sapply(umtaglevel, function(x) get_devel_sum(x))
    ggplot(data=tdf) +
        geom_tile(mapping=aes(x=xtag, y= mtag_fct, fill=n, color=sample_tag, height=0.6, width=0.6),
        size=0.7) +
    #     geom_point(mapping=aes(x=xtag, y=mtag, size=n, color=origin)) +
        theme(axis.text.x = element_text(angle=90, size=25),
        legend.position='NA',
        axis.title.y = element_text(size=30),
        axis.text.y = element_text(size=15,color=umtagcol, face='bold'),
        legend.text=element_text(size=20),
        legend.title = element_blank(),
        legend.key.width = unit(2.5, 'cm')) +
        scale_color_manual(values=c(igemcol_large, igemcol, nativecol)) +
        scale_fill_gradient(low='white', high='red') +
        labs(x=' ', y='Developability encoding')
    outpng('antigen', 'developability_encoding_sln', width=30, height=15)


}


bin_tile_sln_permuted = function(){
    infile = 'eleven_outfiles/antigens_affinity_generated_devel_native_develbin_sln_permuted.csv'
    df = read_csv(infile)
    print(df)
    df$mtag = sprintf('%s_%s_%s_%s_%s_%s', df$charge_at_7_binned, df$mol_weight_binned, df$gravy_binned,
    df$instability_index_binned, df$length_binned, df$mhc2_binned, df$mhc_binned)
    df = fixfive(df)
    df$xtag= sprintf('%s_%s', df$antigen, df$sample_tag)
    print(unique(df$xtag))
    print(unique(df$mtag))
    tdf= df %>% count(xtag,mtag,sample_tag)
    tdf = tdf %>% group_by(xtag) %>% mutate(xtagcount = n())
    tdf = tdf %>% group_by(mtag) %>% mutate(mtagcount = n())
    tdf = tdf %>% group_by(mtag) %>% mutate(mtagsum = sum(n))
    tdf$mtag_fct = factor(tdf$mtag)
    tdf$mtag_fct=fct_reorder2(tdf$mtag_fct, tdf$mtagcount, tdf$mtagsum)
    get_devel_sum = function(x) {
        sumbin = sum(as.numeric(unlist(strsplit(x, split='_'))))
        sumbin = sumbin+1
        return(sumbin)
    }
    get_devel_col = function(x) {
        sumbin = sum(as.numeric(unlist(strsplit(x, split='_'))))
        sumbin = sumbin+1
        return(my_gradient[sumbin])
    }
    tdf$sumbin = sapply(tdf$mtag, function(x) get_devel_sum(x))
    umtag = tdf$mtag_fct
    umtag2 = tdf$mtag
    umtaglevel = levels(umtag)
    umtagcol = sapply(umtaglevel, function(x) get_devel_col(x))
    umtagsum = sapply(umtaglevel, function(x) get_devel_sum(x))
    xtagcol=unique(tdf$xtag)
    xtagcol[grepl('large',xtagcol)]='orange'
    xtagcol[grepl('native_sized',xtagcol)]=igemcol
    xtagcol[grepl('native',xtagcol)]=nativecol
    ggplot(data=tdf) +
        geom_tile(mapping=aes(x=xtag, y= mtag_fct, fill=n, color=sample_tag, height=0.6, width=0.6),
        size=0.7) +
    #     geom_point(mapping=aes(x=xtag, y=mtag, size=n, color=origin)) +
        theme(axis.text.x = element_text(angle=90, size=25, color=xtagcol, face='bold'),
        legend.position='bottom',
        axis.title.y = element_text(size=30),
        axis.text.y = element_text(size=15,color=umtagcol, face='bold'),
        legend.text=element_text(size=20),
        legend.title = element_blank(),
        legend.key.size = unit(2.5, 'cm')) +
        scale_color_manual(values=c(igemcol_large, igemcol, nativecol)) +
        scale_fill_gradient(low='white', high='red') +
        labs(x=' ', y='Developability encoding')
    outpng('antigen', 'developability_encoding_sln_permuted', width=34, height=15)


}


bin_bar_sln_permuted = function(){
    infile = 'eleven_outfiles/antigens_affinity_generated_devel_native_develbin_sln_permuted.csv'
    df = read_csv(infile)
    df$mtag = sprintf('%s_%s_%s_%s_%s_%s', df$charge_at_7_binned, df$mol_weight_binned, df$gravy_binned,
    df$instability_index_binned, df$length_binned, df$mhc2_binned, df$mhc_binned)
    df$xtag= sprintf('%s_%s', df$antigen, df$sample_tag)
    print(unique(df$xtag))
    print(unique(df$mtag))
    bardf= df %>% count(xtag,mtag,sample_tag,nsample) %>%count(xtag,sample_tag,nsample)
    medbardf = bardf %>% group_by(xtag,sample_tag) %>% summarize(medval=median(n),sdval=sd(n))
    print(medbardf)
    # stop()
    ggplot(data=medbardf) +
        geom_bar(mapping=aes(x=xtag, y=medval, fill= sample_tag), stat='identity')+
        geom_text(mapping=aes(x=xtag, y=medval, label=medval), nudge_y=2, size=7) +
        geom_errorbar(mapping=aes(x=xtag, y=medval, ymin=medval-sdval,ymax=medval+sdval)) +
        scale_fill_manual(values=c(igemcol_large, igemcol,nativecol)) +
        labs(y='# of developability combinations') +
        theme(axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        legend.position='NA',
        axis.text.y= element_text(size=30),
        axis.title.y=element_text(size=30))
    outpng('antigen', 'developability_encoding_sln_bar_permuted', width=30)


}

get_devel_combinations_cor = function(){
    infile = 'eleven_outfiles/antigens_affinity_generated_devel_native_develbin_sln_permuted_sampled.csv'
    df = read_csv(infile)
    print(df)
    df$mtag = sprintf('%s_%s_%s_%s_%s_%s', df$charge_at_7_binned, df$mol_weight_binned, df$gravy_binned,
    df$instability_index_binned, df$length_binned, df$mhc2_binned, df$mhc_binned)
    df$xtag= sprintf('%s_%s', df$antigen, df$sample_tag)
    print(unique(df$xtag))
    print(unique(df$mtag))
    countdf = df %>% group_by(xtag,mtag, antigen, sample_tag, nsample) %>% summarize(countval = n())
    ndf = countdf[countdf$sample_tag=='native',]
    print(unique(ndf$nsample))
    print(unique(countdf$sample_tag))
    antigens = unique(countdf$antigen)
    nsamples = unique(countdf$nsample)
    print(countdf)
    outdf = data.frame()
    for (antigen in antigens){
        for (nsample in nsamples){
            adf = countdf[(countdf$antigen==antigen)&(countdf$nsample==nsample),]
            nadf = ndf[ndf$antigen==antigen,]
            nadf$nsample = rep(nsample,dim(nadf)[1])
            gadf = adf[adf$sample_tag=='generated_large',]
            sadf = adf[adf$sample_tag=='generated_native_sized',]
            mngdf = merge(nadf, gadf, by='mtag', all=TRUE)
            mnsdf = merge(nadf, sadf, by='mtag', all=TRUE)
            ngcorval = cor(mngdf$countval.x, mngdf$countval.y, use='pairwise.complete.obs')
            nscorval = cor(mnsdf$countval.x, mnsdf$countval.y, use='pairwise.complete.obs')
            ngcordf = data.frame(corvar='native-generated_large',corval=ngcorval,nsample=nsample,
                                antigen=antigen)
            nscordf = data.frame(corvar='native-generated_native_sized',corval=nscorval,nsample=nsample,
                                antigen=antigen)
            outdf= rbind(outdf,ngcordf)
            outdf = rbind(outdf,nscordf)
        }
    }
    outdf = drop_na(outdf)
    outname='eleven_outfiles/antigens_affinity_generated_devel_native_develbin_sln_permuted_cor.csv'
    write.csv(outdf, outname, row.names=FALSE)

}

plot_devel_combinations_cor = function(){
    infile = 'eleven_outfiles/antigens_affinity_generated_devel_native_develbin_sln_permuted_cor.csv'
    df = read_csv(infile)
    print(df)
    print(min(df$corval))
    print(max(df$corval))
    stop()
    ggplot(data=df) +
        geom_boxplot(mapping=aes(x=corvar, y=corval)) +
        geom_point(mapping=aes(x=corvar, y=corval), position=position_dodge2(), size=5, shape=21) +
        facet_wrap(~antigen, nrow=2) +
        theme(axis.text.x = element_text(angle=90))
    outpng('tmp', 'tmp')
}

# run stuff
#viz_binder_bins()
# viz_binder_bins2() # panel A
# bin_bar2('eleven_outfiles/antigens_affinity_generated_native_develbin.csv')
# bin_bar_large_generation()
# bin_bar_sln() # panel B top
# bin_tile_sln() # panel B bottom
#bin_tile_sln_permuted() # panel B bottom
# bin_bar_sln_permuted() # panel B top
# get_devel_combinations_cor()
# plot_devel_combinations_cor()
