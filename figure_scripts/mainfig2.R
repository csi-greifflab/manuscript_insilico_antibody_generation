# Title     : TODO
# Objective : TODO
# Created by: rahmadakbar
# Created on: 2020-05-08

# main figure 2 of the manuscript

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




get_top_bot_overlap = function(infile, rank_tag){
    df = read_csv(infile)
    pdbids = unique(df$pdbid)
    df = df[df$rank==rank_tag,]
    source='igem'
    if (grepl('native', infile)){
        df=df %>% rename(col2='CDR3', col3='Slide', col4='Energy', col7='Structure')
        source='native'
    }
    groupeddf = group_by(df, pdbid) %>% group_split()
    get_seq=function(x){
        seq = x$CDR3
        return(seq)
    }
    get_tag=function(x){
        pdbid =x$pdbid
        rank=x$rank
        source=x$source
        # tag = sprintf('%s_%s_%s', pdbid, source, rank)
        tag = pdbid
        return(unique(tag))
    }
    seq_list = lapply(groupeddf, get_seq)
    tags = sapply(groupeddf, get_tag)
    names(seq_list)= tags
    filename=sprintf('figures/overlap_%s_%s.tiff', rank_tag,source)
    print(filename)
    venn.diagram(seq_list,
    filename=filename,
    # fill = c(nativecol, igemcol),
    fill = my_spectral[1:length(pdbids)],
    cat.col = my_spectral[1:length(pdbids)],
    col='NA',
    cat.default.pos = 'outer',
    scaled=TRUE,
    margin =0.1)
    command = sprintf('open %s', filename)
    system(command)
}



ld_dist_cutoff = function(){
    infile = 'outfiles/antigens_native_delicab_rank_develparams_ld_1pc_merged.csv'
    # infile = 'outfiles/antigens_native_delicab_rank_develparams_ld_10pc_merged.csv'
    df = read_csv(infile)
    df$source[df$source=='igem'] = 'generated'
    # df = df[(df$pdbid1==df$pdbid2) & (df$rank1==df$rank2),]
    meddf = df %>% group_by(pdbid1, rank1, source) %>% summarize(medval=median(ld))
    meddf$mlabel = sprintf('med: %s', meddf$medval)
    meddf$ypos = rep(1,dim(meddf)[1])
    meddf$ypos[meddf$source=='igem'] = 0.9
    df = df %>% group_by(pdbid1, rank1, ld, source) %>% summarize(size=n())
    normalit = function(m){
        (m - min(m))/(max(m)-min(m))
    }
    cumul_percent=function(m){
        cumsum(m)/sum(m)
    }
    df = df %>% group_by(pdbid1, rank1, source) %>% mutate(size_norm=normalit(size), size_cum=cumul_percent(size))
    df$tag = sprintf('%s_%s_%s',df$pdbid1, df$rank1, df$source)
    print(df[1:20,])
    # ggplot(data=df) +
    #     geom_bar(mapping=aes(x=ld, y=size_cum, fill=source), stat='identity', position='dodge2') +
    #     geom_text(data=meddf, mapping=aes(x=0, y=ypos, label=mlabel, color=source), size=3, hjust=0) +
    #     facet_grid(rank1~pdbid1) +
    #     theme(axis.text=element_text(size=20),
    #     axis.title=  element_text(size=20),
    #     legend.title = element_blank(),
    #     legend.text= element_text(size=20))+
    #     scale_fill_manual(values = c(igemcol, nativecol)) +
    #     scale_color_manual(values = c(igemcol, nativecol))
    # labs(y='cumulative fraction', x='LD')
    #
    # outpng(infile, 'tmp', width=12)
    ld_cutoff = 15
    df = df[(df$ld<=ld_cutoff) & (df$rank1=='top'),]
    df$source = factor(df$source, levels = c('native', 'generated'))
    print(df)
    meddf = df %>% group_by(pdbid1, rank1, source) %>% summarize(medval=median(ld))
    meddf$mlabel = sprintf('median: %s', meddf$medval)
    meddf$ypos = rep(1,dim(meddf)[1])
    meddf$ypos[meddf$source=='generated'] = 0.9
    ggplot(data=df) +
        geom_bar(mapping=aes(x=ld, y=size_cum, fill=source), stat='identity', position='dodge2') +
        geom_text(data=meddf, mapping=aes(x=0, y=ypos, label=mlabel, color=source), size=5, hjust=0) +
        # facet_grid(rank1~pdbid1) +
        facet_wrap(~pdbid1, nrow=2) +
        theme(axis.text=element_text(size=20),
        axis.title=element_text(size=25),
        legend.title = element_blank(),
        legend.text = element_text(size=20))+
        labs(y='Cumlative fraction of distance\nbetween native and generated', x='Levenshtein distance (LD)') +
        scale_fill_manual(values = c(nativecol, igemcol)) +
        scale_color_manual(values = c(nativecol, igemcol))
    outpdf('native_igem_diversity_ld', ld_cutoff, width=12)
}



get_native_igem_kmercount_cor_m = function(){
    infile = 'outfiles/antigens_native_delicab_rank_develparams_gapped_kmer_trimmed_m.csv'
    infile2 = 'outfiles/antigens_geminio_absolute_rank_develparams_gapped_kmer_trimmed_m.csv'
    df = read_csv(infile)
    df2 = read_csv(infile2)
    df = df[df$rank != 'all',]
    df2 = df2[df2$rank != 'all',]
    print(df[df$m==2,])
    print(df2[df2$m==2,])
    mdf = merge(df, df2, by='kmer', all=TRUE)
    mdf = na.omit(mdf)
    print(unique(mdf$pdbid.x))
    mdf2 = mdf[(mdf$rank.x==mdf$rank.y) & (mdf$pdbid.x==mdf$pdbid.y) & (mdf$m.x==mdf$m.y),]
    mdf3 = mdf[(mdf$rank.x==mdf$rank.y) & (mdf$pdbid.x==mdf$pdbid.y),]
    print(dim(mdf2))
    print(dim(mdf3))
    cordf = mdf2 %>% group_by(pdbid.x, rank.x, m.x) %>% summarize(corval = cor(count.x, count.y))
    cordf$corlabel = sprintf('Cor (P): %s', round(cordf$corval,3))
    print(unique(cordf$pdbid.x))
    cordf = cordf[cordf$rank.x=='top',]
    ggplot(cordf) +
        geom_bar(mapping=aes(x=m.x, y=corval), fill=my_spectral[1], stat='identity') +
        geom_text(data=cordf, mapping=aes(x=m.x, y=corval, label=round(corval,3)), angle=90,
        position=position_nudge(y=-0.1), color='white') +
        # facet_grid(rank.x~pdbid.x, scale='free') +
        facet_grid(~pdbid.x, scale='free') +
        theme(axis.title = element_text(size=25)) +
        labs(x='m', y = 'Cor. coef. (Pearson)')
    outpdf('native', 'igem_cor_m_top_bar')
}


get_native_igem_kmercount_cor = function(){
    infile = 'outfiles/antigens_native_delicab_rank_develparams_gapped_kmer_trimmed.csv'
    infile2 = 'outfiles/antigens_geminio_absolute_rank_develparams_gapped_kmer_trimmed.csv'
    df = read_csv(infile)
    df2 = read_csv(infile2)
    df = df[df$rank != 'all',]
    df2 = df2[df2$rank != 'all',]
    print(df)
    print(df2)
    # mdf = tibble(merge(df, df2, by='kmer'))
    mdf = merge(df, df2, by='kmer', all=TRUE)
    column_names = colnames(mdf)
    # mdf = tibble(mdf)
    mdf = mdf[(mdf$rank.x==mdf$rank.y) & (mdf$pdbid.x==mdf$pdbid.y),]
    # ggplot(mdf[mdf$rank.x=='all',]) +
    #     geom_point(mapping=aes(x=count.x, y=count.y)) +
    #     facet_grid(rank.x~pdbid.x, scale='free') +
    #     theme(axis.text.x = element_text(angle=90))
    # outpng('tmp', 'tmp')
    # topbotdf = mdf[mdf$rank.x!='all', ]
    topdf = mdf[mdf$rank.x=='top', ]
    # print(head(topbotdf, 10))
    # cordf = topbotdf %>% group_by(pdbid.x, rank.x) %>% summarize(corval = cor(count.x, count.y))
    cordf = mdf %>% group_by(pdbid.x, rank.x) %>% summarize(corval = cor(count.x, count.y))
    cordf$corlabel = sprintf('Cor (P): %s', round(cordf$corval,3))
    print(cordf)
    cordf = cordf[cordf$rank.x=='top',]
    ggplot(topdf) +
        geom_point(mapping=aes(x=count.x, y=count.y), color=my_spectral[1]) +
        geom_text(data=cordf, mapping=aes(x=0, y=10000, label=corlabel), hjust=0) +
        # facet_grid(rank.x~pdbid.x, scale='free') +
        facet_grid(~pdbid.x, scale='free') +
        theme(axis.text.x = element_text(angle=90),
            axis.title = element_text(size=25)) +
        labs(x='# of k-mer in native (k=1, m=2)', y = '# of k-mer in generated (k=1, m=2)')
    outpng('native', 'igem_top_cor')
}


binding_affinity_igem_native_dist_grid = function(){
    igem_antigenfile = 'outfiles/antigens_geminio_absolute_rank_develparams.csv'
    native_antigenfile = 'outfiles/antigens_native_delicab_rank_develparams.csv'
    igemdf = read_csv(igem_antigenfile)
    natdf  = read_csv(native_antigenfile)
    mergednames= c(names(igemdf)[2],names(igemdf)[4:dim(igemdf)[2]])
    natdf=natdf %>% rename(col2='CDR3', col3='Slide', col4='Energy', col7='Structure')
    igemdf = igemdf[,mergednames]
    natdf = natdf[, mergednames]
    mdf=rbind(natdf, igemdf)
    mdf$source[mdf$source=='igem'] = 'generated'
    mdf$seq_origin = sprintf('%s_%s', mdf$pdbid, mdf$source)
    print(mdf)
    gdf = gather(mdf, key='devel_params', value = 'value', Energy, charge_at_7, mol_weight, gravy,
    instability_index, length, netmhc2pan_averank)
    meddf = gdf %>% group_by(pdbid,seq_origin, devel_params, rank) %>% summarize(med_value = median(value), max_value=max(value), min_value=min(value))
    gdf = gdf[(gdf$devel_params=='Energy')&(gdf$rank=='top'),]
    meddf = meddf[(meddf$devel_params=='Energy')&(meddf$rank=='top'),]
    print(gdf)
    # stop()
    ggplot(data=gdf) +
        geom_density_ridges(mapping = aes(x=value, y=seq_origin, color=source), alpha=0.5, size=2)+
        # facet_grid(rank~devel_params, scales = 'free_x') +
        geom_text(data=meddf, mapping = aes(x = med_value, y=seq_origin,
        label = sprintf('median: %s', round(meddf$med_value, 2))), color='black',size=6, nudge_y=0.3, hjust=0.5)+
        # scale_color_manual(values = my_spectral)+
        scale_color_manual(values = c(igemcol, nativecol))+
        theme(axis.text=element_text(size=25),
        axis.title=element_text(size=25),
        legend.text = element_text(size=25),
        legend.title = element_text(size=25),
        strip.text = element_text(size=20))+
        labs(x='Binding energy', y='') +
        theme(legend.position = 'none')
    outpng('native', 'igem_energy_density')
}


develparams_igem_native_dist_grid = function(){
    igem_antigenfile = 'outfiles/antigens_geminio_absolute_rank_develparams.csv'
    native_antigenfile = 'outfiles/antigens_native_delicab_rank_develparams.csv'
    igemdf = read_csv(igem_antigenfile)
    natdf  = read_csv(native_antigenfile)
    mergednames= c(names(igemdf)[2],names(igemdf)[4:dim(igemdf)[2]])
    natdf=natdf %>% rename(col2='CDR3', col3='Slide', col4='Energy', col7='Structure')
    igemdf = igemdf[,mergednames]
    natdf = natdf[, mergednames]
    mdf=rbind(natdf, igemdf)
    mdf$source[mdf$source=='igem'] = 'generated'
    mdf$seq_origin = sprintf('%s_%s', mdf$pdbid, mdf$source)
    print(mdf)
    gdf = gather(mdf, key='devel_params', value = 'value', Energy, charge_at_7, mol_weight, gravy,
    instability_index, length, netmhc2pan_averank)
    meddf = gdf %>% group_by(pdbid,seq_origin, devel_params, rank) %>% summarize(med_value = median(value), max_value=max(value), min_value=min(value))
    gdf = gdf[(gdf$devel_params!='Energy')&(gdf$rank=='top'),]
    meddf = meddf[(meddf$devel_params!='Energy')&(meddf$rank=='top'),]
    devel_labels = c('charge_at_7'='Charge pH (7)',
                    'gravy'='Gravy',
                    'instability_index'='Instability index',
                    'length'='CDR-H3 Length',
                    'mol_weight'='Mol. weight',
                    'netmhc2pan_averank'='NetMHCIIpan')
    print(unique(gdf$devel_params))
    print(devel_labels)
    ggplot(data=gdf) +
        geom_density_ridges(mapping = aes(x=value, y=seq_origin, color=source), alpha=0.7, size=2)+
        facet_wrap(~devel_params, scales = 'free_x', nrow=1, labeller=as_labeller(devel_labels)) +
        geom_text(data=meddf, mapping = aes(x = med_value, y=seq_origin,
        label = sprintf('median: %s', round(meddf$med_value, 2))), color='black',size=6, nudge_y=0.3, hjust=0.2)+
        # scale_color_manual(values = my_spectral)+
        scale_color_manual(values = c(igemcol, nativecol))+
        theme(axis.text=element_text(size=25),
        axis.title=element_text(size=25),
        legend.text = element_text(size=25),
        legend.title = element_text(size=25),
        strip.text = element_text(size=30))+
        labs(x='', y='') +
        theme(legend.position = 'none')
    outpng('native', 'igem_develparam_density', width = 27, height=15)
}


show_colors = function(){
    df = data_frame(colors = my_base_color2, fake_y = rep(1, length(my_base_color2)))
    # df = data_frame(colors = my_base_color, fake_y = rep(1, length(my_base_color)))
    print(df)
    ggplot(data=df) +
        geom_bar(mapping=aes(x=colors, y=fake_y, fill=colors), stat='identity') +
        scale_fill_manual(values=my_base_color2) +
        theme(axis.text.x = element_text(angle=90))
    outpdf('my_base_color', '2_palette')
}



get_pairwise_overlap = function(){

    igemfile = 'outfiles/antigens_geminio_absolute_rank_develparams.csv'
    natfile = 'outfiles/antigens_native_delicab_rank_develparams.csv'
    # tags = c('1OAZ_native_bot', '1OAZ_native_top', '1FBI_native_bot')
    get_pdbid_from_tag = function(tag){
        print(tag)
        pdbid = strsplit(tag, '_')[[1]][1]
        return(pdbid)
    }
    igemdf = read_csv(igemfile)
    # igemdf = read_csv(igemfile, n_max=100)
    natdf = read_csv(natfile)
    # natdf = read_csv(natfile, n_max=100)
    mergednames= c(names(igemdf)[2],names(igemdf)[4:dim(igemdf)[2]])
    natdf=natdf %>% rename(col2='CDR3', col3='Slide', col4='Energy', col7='Structure')
    igemdf = igemdf[,mergednames]
    natdf = natdf[, mergednames]
    mdf  = rbind(natdf,igemdf)
    # mdf = mdf[mdf$rank!= 'all',]
    mdf = mdf[mdf$rank=='top',]
    mdf$source[mdf$source=='igem'] ='generated'
    groupeddf = group_by(mdf, pdbid, source, rank) %>% group_split()
    get_seq=function(x){
        seq = x$CDR3
        return(seq)
    }
    get_tag=function(x){
        pdbid =x$pdbid
        rank=x$rank
        source=x$source
        # tag = sprintf('%s_%s_%s', pdbid, source, rank)
        tag = sprintf('%s_%s', pdbid, source)
        return(unique(tag))
    }
    listIntersect <- function(inList) {
        X <- crossprod(table(stack(inList)))
        #X[lower.tri(X)] <- NA
        print(diag(X))
        print(X[1,])
        print(dim(X))
        X2 = copy(X)
        X2 = t(apply(X2,1, function(x) round(x/max(x)*100,2)))
        print(X2[1,])
        X[lower.tri(X)] <- X2[lower.tri(X2)]
        print(X)
        print(class(X))
        # diag(X) <- NA
        out <- na.omit(data.frame(as.table(X)))
        out[order(out$ind), ]
    }
    seqs = lapply(groupeddf, get_seq)
    tags = sapply(groupeddf, get_tag)
    names(seqs)= tags
    print(seq)
    overlaps = listIntersect(seqs)
    overlaps = overlaps[order(overlaps$ind),]
    colnames(overlaps) = c('tag1', 'tag2', 'freq')
    overlaps$freq2 = round(overlaps$freq/10000,3)
    overlaps$xpdbid = sapply(as.character(overlaps$tag1), get_pdbid_from_tag)
    overlaps$ypdbid = sapply(as.character(overlaps$tag2), get_pdbid_from_tag)
    # print(overlaps)
    # overlaps = spread(overlaps, tag2, freq)
    # print(overlaps)
    pdbids = unique(overlaps$xpdbid)
    pdbcolors = my_spectral[1:length(pdbids)]
    for (i in 1:length(pdbids)){
        overlaps$xpdbid[overlaps$xpdbid==pdbids[i]] = pdbcolors[i]
        overlaps$ypdbid[overlaps$ypdbid==pdbids[i]] = pdbcolors[i]
    }
    overlaps$xsource = sapply(as.character(overlaps$tag1), function(x) strsplit(x, split='_')[[1]][2])
    overlaps$ysource = sapply(as.character(overlaps$tag2), function(x) strsplit(x, split='_')[[1]][2])
    sources = unique(overlaps$xsource)
    sourcecolors = c(igemcol, nativecol)
    for (i in 1:length(sources)){
        overlaps$xsource[overlaps$xsource==sources[i]] = sourcecolors[i]
        overlaps$ysource[overlaps$ysource==sources[i]] = sourcecolors[i]
    }
    print(overlaps)
    ggplot(data=overlaps)+
        geom_tile(mapping=aes(x=tag1, y=tag2, fill=freq2)) +
        theme(axis.text.y = element_text(angle=0, color=overlaps$ysource, face='bold', size=15),
        # axis.text.y=element_text(color=overlaps$ysource, face='bold'),
        axis.text.x=element_text(angle=90, color=overlaps$ysource,face='bold', size=15),
        # axis.text.x=element_blank(),
        legend.position='None',
        axis.title.x=element_text(size=20)) +
        guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
        title.position = "top", title.hjust = 0.5))+
        geom_text(mapping=aes(x=tag1, y=tag2, label=freq), color='white', fontface='bold', size=6) +
        labs(y=' ', x='Overlap')
    outpdf('igem_native', 'top_bot_overlap_heatmap2', width=10)
    # filename='figures/temp.tiff'
    # venn.diagram(seqs[1:4],
    #     filename=filename,
    #     fill = c(my_spectral[1], my_spectral[5], my_spectral[9], my_spectral[14]),
    #     cat.col = c(my_spectral[1], my_spectral[5], my_spectral[9], my_spectral[14]),
    #     col='NA',
    #     cat.default.pos = 'outer',
    #     scaled=TRUE,
    #     margin =0.1)
    # command = sprintf('open %s', filename)
    # system(command)
}


do_jensen_shannon = function(){
    infile ='outfiles/generator_seven_absolute_merged_jsval.csv'
    df = read_csv(infile)
    df$nsample = factor(df$nsample, level=c(700,7000,70000))
    df$label = sprintf('%s (%s)', df$nsample, round(df$jsval,2))
    print(df)
    ggplot(data=df, mapping=aes(x=pdbid, y=jsval, fill=nsample))+
        geom_bar(stat='identity', position='dodge2') +
        geom_text(mapping=aes(label=label), position=position_dodge(width=1), vjust=-0.3) +
        labs(y='Jensen Shannon distance (native vs generated)', x='Antigens') +
        scale_fill_manual(values=c(igemcol,igemcol,igemcol)) +
        theme(legend.position='None',
                axis.text=element_text(size=20),
                axis.title=element_text(size=20))
    outpdf(infile, 'jsdistance')


    # print(unique(df$))

}

binding_affinity_native_generated = function(){
    ginfile = 'eleven_outfiles/antigens_affinity_generated_native_sized.csv'
    ninfile = 'eleven_outfiles/antigens_affinity_native.csv'
    gdf = read_csv(ginfile)
    ndf = read_csv(ninfile)
    df = rbind(gdf, ndf)
    df = fixfive(df)
    df$antigen_origin = sprintf('%s_%s', df$antigen, df$origin)
    meddf = df %>% group_by(antigen_origin, origin, antigen) %>% summarize(medval=round(median(Energy),1))
    meddf2 = meddf %>% group_by(antigen) %>% summarize(medval = mean(medval))
    print(meddf2)
    (antigen_order = structure(meddf2$medval, names=meddf2$antigen))
    print(str(antigen_order))
    antigen_order2 = sapply(df$antigen, function(x) antigen_order[x])
    print(antigen_order2)
    meddf$medlabel = sprintf('median: %s', meddf$medval)
    df$antigen_order = antigen_order2
    axiscol = ifelse(meddf$origin=='generated', igemcol, nativecol)
    print(axiscol)
    ggplot(data=df) +
        geom_density_ridges(mapping=aes(x=Energy, y=reorder(antigen_origin, antigen_order), color=origin), size=1.3) +
        geom_text(data=meddf, mapping=aes(x=-60, y=antigen_origin, label=medlabel, color=origin),nudge_y=0.3, size=5) +
        scale_color_manual(values=c(igemcol, nativecol)) +
        labs(y='CDR-H3 origin', x='Affinity (Binding energy)') +
        theme(legend.position='none',
                axis.text = element_text(size=15),
                axis.title=element_text(size=20),
                axis.text.y = element_text(colour=axiscol, face='bold'))
    outpng('native_generated', 'binding_affinity_distribution')
}


binding_fold_cor = function(){
    infile = 'eleven_outfiles/antigens_native_generated_structure_count.csv'
    df  = read_csv(infile)
    df = fixfive(df)
    nstructdf = df %>% group_by(antigen) %>% summarize(count=n())
    print(nstructdf)
    df = na.omit(df)
    # print(df)
    cordf = df %>% group_by(antigen) %>% summarize(corval=cor(native_count, generated_count),
                                                    count=n(),
                                                    yloc= max(generated_count),
                                                    xloc=min(native_count))
    cordf$label = sprintf('Cor.(P): %s', round(cordf$corval,3))
    print(cordf)
    ggplot(data=df) +
        geom_point(mapping=aes(x=native_count, y=generated_count),size=5, shape=21,color=my_base_color2[5]) +
        geom_text(data=cordf, mapping=aes(x=xloc, y=yloc*0.9, label=label), hjust=0, fontface='bold', size=7) +
        facet_wrap(~antigen, nrow=2, scale = 'free') +
        labs(x='# of paratope folds in native', y='# of paratope folds\nin generated') +
        theme(axis.title=element_text(size=40),
            strip.text=element_text(size=20),
            axis.text.x= element_text(size=20,angle=90, color=nativecol),
            axis.text.y=element_text(size=20, color=igemcol))
    outpng('native_generated', 'structure_count_cor', width=17, height=7)

}



binding_fold_bar = function(){
    infile = 'eleven_outfiles/antigens_native_generated_structure_count.csv'
    df  = read_csv(infile)
    df = fixfive(df)
    nstructdf = df %>% group_by(antigen) %>% summarize(count=n())
    print(nstructdf)
    df = na.omit(df)
    # print(df)
    gdf = gather(df, 'origin', 'n_found', -Structure,-antigen)
    print(gdf)
    cdf = gdf[gdf$n_found!=0,] %>% group_by(antigen, origin) %>% summarize(n_structure=n())
    print(cdf)
    cdf$ypos = 25000
    cdf$ypos[cdf$origin=='native_count'] = 10000
    cdf$xpos = 10
    cdf$strlabs = sprintf('Total: %s', cdf$n_structure)
    ggplot(data=gdf) +
        geom_bar(mapping=aes(x=reorder(Structure,-n_found), y=n_found, fill=origin), stat='identity') +
        geom_text(data=cdf, mapping=aes(x=xpos, y=ypos, label=strlabs, color=origin), hjust=0, fontface='bold',
        size=7) +
        facet_wrap(~antigen, nrow=2, scale = 'free') +
        labs(x='Paratope folds', y='# of paratope folds') +
        theme(axis.title=element_text(size=40),
        strip.text=element_text(size=20),
        axis.text.x= element_text(size=2,angle=90, color='black'),
        axis.text.y=element_text(size=20, color='black'),
        legend.position='NA') +
        scale_color_manual(values=c(igemcol, nativecol)) +
        scale_fill_manual(values=c(igemcol, nativecol))
    outpng('native_generated', 'structure_count_bar', width=16, height=7)

}


kmer_native_generated_cor = function(){
    ginfile = 'eleven_outfiles/antigens_affinity_native_gapped_kmer_trimmed.csv'
    ninfile = 'eleven_outfiles/antigens_affinity_generated_native_sized_gapped_kmer_trimmed.csv'
    gdf = read_csv(ginfile)
    ndf = read_csv(ninfile)
    print(gdf)
    print(ndf)
    mdf =  merge(gdf, ndf, by='kmer', all=TRUE)
    mdf = mdf[mdf$antigen.x==mdf$antigen.y,]
    mdf = fixfive3(mdf)
    print(mdf)
    print(colnames(mdf))
    cordf = mdf %>% group_by(antigen.x) %>% summarize(corval = cor(count.x, count.y))
    cordf$corlabel = sprintf('Cor. (P): %s', round(cordf$corval,3))
    print(cordf)
    ggplot(mdf) +
        geom_point(mapping=aes(x=count.x, y=count.y), color=my_base_color2[5], size=5, shape=21) +
        geom_text(data=cordf, mapping=aes(x=0, y=50000, label=corlabel), hjust=0, fontface='bold', size=4) +
    # facet_grid(rank.x~pdbid.x, scale='free') +
        facet_wrap(~antigen.x, scale='free', nrow=2) +
        theme(axis.text.x = element_text(angle=90, color=nativecol),
            axis.text.y=element_text(color=igemcol),
            axis.title = element_text(size=25),
            strip.text = element_text(size=20),
            axis.text = element_text(size=20)) +
        labs(x='# of k-mer in native (k=1, m=5)', y = '# of k-mer in generated (k=1, m=5)')
    outpng('native_generated', 'kmer_decomposition_cor', width = 12)
}

devel_params_native_generated = function(){
    ginfile = 'eleven_outfiles/antigens_affinity_generated_native_sized_devel.csv'
    ninfile = 'eleven_outfiles/antigens_affinity_native_devel.csv'
    gdf = read_csv(ginfile)
    ndf = read_csv(ninfile)
    print(gdf)
    print(ndf)
    mdf = rbind(gdf,ndf)
    mdf = fixfive(mdf)
    print(mdf)
    gdf = gather(mdf, key='devel_params', value = 'value', charge_at_7, mol_weight, gravy,
    instability_index, length, mhc2, mhc)
    gdf = gdf[gdf$devel_params!='length',]
    meddf = gdf %>% group_by(antigen,origin, devel_params) %>% summarize(med_value = median(value),
    max_value=max(value), min_value=min(value))
    print(gdf)
    # print(meddf)
    gdf$antigen_origin = sprintf('%s_%s', gdf$antigen, gdf$origin)
    meddf$antigen_origin = sprintf('%s_%s', meddf$antigen, meddf$origin)
    devel_labels = c('charge_at_7'='Charge pH (7)',
                'gravy'='Gravy',
                'instability_index'='Instability index',
                # 'length'='CDR-H3 Length',
                'mol_weight'='Mol. weight',
                'mhc2'='NetMHCIIpan',
                'mhc'='NetMHCpan')
    print(gdf)
    print(meddf)
    antigen_origins =unique(meddf$antigen_origin)
    print(antigen_origins)
    yaxiscol = ifelse(grepl('generated',antigen_origins), igemcol, nativecol)
    print(yaxiscol)
    ggplot(data=gdf) +
        geom_density_ridges(mapping = aes(x=value, y=antigen_origin, color=origin), alpha=0.7, size=2)+
        facet_wrap(~devel_params, scales = 'free_x', nrow=1, labeller=as_labeller(devel_labels)) +
        geom_text(data=meddf, mapping = aes(x = min_value, y=antigen_origin,
        label = sprintf('median: %s', round(meddf$med_value, 3))), color='black',size=10, nudge_y=0.3, hjust=0)+
        # scale_color_manual(values = my_spectral)+
        scale_color_manual(values = c(igemcol, nativecol))+
        theme(axis.text=element_text(size=30),
            axis.text.y = element_text(color=yaxiscol),
            axis.title=element_text(size=25),
            legend.text = element_text(size=25),
            legend.title = element_text(size=25),
            strip.text = element_text(size=35))+
        labs(x='', y='') +
        theme(legend.position = 'none')
    outpng('native_generated', 'devel_distr', width = 30, height=15)
}

pairwise_overlap_native_generated = function(){
    ginfile = 'eleven_outfiles/antigens_affinity_generated_native_sized_devel.csv'
    ninfile = 'eleven_outfiles/antigens_affinity_native_devel.csv'
    gdf = read_csv(ginfile)
    ndf = read_csv(ninfile)
    get_pdbid_from_tag = function(tag){
        print(tag)
        pdbid = strsplit(tag, '_')[[1]][1]
        return(pdbid)
    }
    sgdf = gdf %>% group_by(antigen) %>% summarize(nseq=n())
    sndf = ndf %>% group_by(antigen) %>% summarize(nseq=n())
    print(sgdf)
    print(sndf)
    # stop()
    mdf = rbind(gdf, ndf)
    groupeddf = group_by(mdf, antigen, origin) %>% group_split()
    get_seq=function(x){
        seq = x$CDR3
        return(seq)
    }
    get_tag=function(x){
        antigens =x$antigen
        origins =x$origin
        tag = sprintf('%s_%s', antigens, origins)
        return(unique(tag))
    }
    listIntersect <- function(inList) {
        X <- crossprod(table(stack(inList)))
        #X[lower.tri(X)] <- NA
        print(diag(X))
        print(X[1,])
        print(dim(X))
        X2 = copy(X)
        X2 = t(apply(X2,1, function(x) round(x/max(x)*100,2)))
        print(X2[1,])
        X[lower.tri(X)] <- X2[lower.tri(X2)]
        print(X)
        print(class(X))
        # diag(X) <- NA
        out <- na.omit(data.frame(as.table(X)))
        out[order(out$ind), ]
    }
    seqs = lapply(groupeddf, get_seq)
    tags = sapply(groupeddf, get_tag)
    names(seqs)= tags
    # print(seqs)
        overlaps = listIntersect(seqs)
    overlaps = overlaps[order(overlaps$ind),]
    colnames(overlaps) = c('tag1', 'tag2', 'freq')
    overlaps$freq2 = round(overlaps$freq/70000,3)
    overlaps$xpdbid = sapply(as.character(overlaps$tag1), get_pdbid_from_tag)
    overlaps$ypdbid = sapply(as.character(overlaps$tag2), get_pdbid_from_tag)
    # print(overlaps)
    # overlaps = spread(overlaps, tag2, freq)
    # print(overlaps)
    pdbids = unique(overlaps$xpdbid)
    pdbcolors = my_spectral[1:length(pdbids)]
    for (i in 1:length(pdbids)){
        overlaps$xpdbid[overlaps$xpdbid==pdbids[i]] = pdbcolors[i]
        overlaps$ypdbid[overlaps$ypdbid==pdbids[i]] = pdbcolors[i]
    }
    overlaps$xsource = sapply(as.character(overlaps$tag1), function(x) strsplit(x, split='_')[[1]][2])
    overlaps$ysource = sapply(as.character(overlaps$tag2), function(x) strsplit(x, split='_')[[1]][2])
    sources = unique(overlaps$xsource)
    sourcecolors = c(igemcol, nativecol)
    for (i in 1:length(sources)){
        overlaps$xsource[overlaps$xsource==sources[i]] = sourcecolors[i]
        overlaps$ysource[overlaps$ysource==sources[i]] = sourcecolors[i]
    }
    print(overlaps)
    ggplot(data=overlaps)+
        geom_tile(mapping=aes(x=tag1, y=tag2, fill=freq2)) +
        theme(axis.text.y = element_text(angle=0, color=overlaps$ysource, face='bold', size=15),
        # axis.text.y=element_text(color=overlaps$ysource, face='bold'),
        axis.text.x=element_text(angle=90, color=overlaps$ysource,face='bold', size=25),
        # axis.text.x=element_blank(),
        legend.position='None',
        axis.title.x=element_text(size=30)) +
        guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
        title.position = "top", title.hjust = 0.5))+
        geom_text(mapping=aes(x=tag1, y=tag2, label=freq), color='white', fontface='bold', size=6) +
        labs(y=' ', x='Overlap')
    outpdf('naitive_generated', 'seq_overlap_heatmap', height=15,width=20)

}



get_pairwise_overlap_native_generated_boxplot = function(){
    ginfile = 'eleven_outfiles/antigens_affinity_generated_native_sized_devel.csv'
    ninfile = 'eleven_outfiles/antigens_affinity_native_devel.csv'
    gdf = read_csv(ginfile)
    ndf = read_csv(ninfile)
    get_pdbid_from_tag = function(tag){
        print(tag)
        pdbid = strsplit(tag, '_')[[1]][1]
        return(pdbid)
    }
    sgdf = gdf %>% group_by(antigen) %>% summarize(nseq=n())
    sndf = ndf %>% group_by(antigen) %>% summarize(nseq=n())
    print(sgdf)
    print(sndf)
    # stop()
    mdf = rbind(gdf, ndf)
    groupeddf = group_by(mdf, antigen, origin) %>% group_split()
    get_seq=function(x){
        seq = x$CDR3
        return(seq)
    }
    get_tag=function(x){
        antigens =x$antigen
        origins =x$origin
        tag = sprintf('%s_%s', antigens, origins)
        return(unique(tag))
    }
    listIntersect <- function(inList) {
        X <- crossprod(table(stack(inList)))
        #X[lower.tri(X)] <- NA
        print(diag(X))
        print(X[1,])
        print(dim(X))
        X2 = copy(X)
        X2 = t(apply(X2,1, function(x) round(x/max(x)*100,2)))
        print(X2[1,])
        X[lower.tri(X)] <- X2[lower.tri(X2)]
        print(X)
        print(class(X))
        # diag(X) <- NA
        out <- na.omit(data.frame(as.table(X)))
        out[order(out$ind), ]
    }
    seqs = lapply(groupeddf, get_seq)
    tags = sapply(groupeddf, get_tag)
    names(seqs)= tags
    # print(seqs)
    overlaps = listIntersect(seqs)
    overlaps = overlaps[order(overlaps$ind),]
    colnames(overlaps) = c('tag1', 'tag2', 'freq')
    overlaps$freq2 = round(overlaps$freq/70000,3)
    overlaps$xpdbid = sapply(as.character(overlaps$tag1), get_pdbid_from_tag)
    overlaps$ypdbid = sapply(as.character(overlaps$tag2), get_pdbid_from_tag)
    # print(overlaps)
    # overlaps = spread(overlaps, tag2, freq)
    # print(overlaps)
    pdbids = unique(overlaps$xpdbid)
    pdbcolors = my_spectral[1:length(pdbids)]
    for (i in 1:length(pdbids)){
        overlaps$xpdbid[overlaps$xpdbid==pdbids[i]] = pdbcolors[i]
        overlaps$ypdbid[overlaps$ypdbid==pdbids[i]] = pdbcolors[i]
    }
    overlaps$xsource = sapply(as.character(overlaps$tag1), function(x) strsplit(x, split='_')[[1]][2])
    overlaps$ysource = sapply(as.character(overlaps$tag2), function(x) strsplit(x, split='_')[[1]][2])
    sources = unique(overlaps$xsource)
    sourcecolors = c(igemcol, nativecol)
    for (i in 1:length(sources)){
        overlaps$xsource[overlaps$xsource==sources[i]] = sourcecolors[i]
        overlaps$ysource[overlaps$ysource==sources[i]] = sourcecolors[i]
    }
    outname = 'eleven_outfiles/pairwise_overlap_native_generated.csv'
    write.csv(overlaps, outname, row.names=FALSE)
    stop()
    ggplot(data=overlaps)+
        geom_tile(mapping=aes(x=tag1, y=tag2, fill=freq2)) +
        theme(axis.text.y = element_text(angle=0, color=overlaps$ysource, face='bold', size=15),
        # axis.text.y=element_text(color=overlaps$ysource, face='bold'),
        axis.text.x=element_text(angle=90, color=overlaps$ysource,face='bold', size=25),
        # axis.text.x=element_blank(),
        legend.position='None',
        axis.title.x=element_text(size=30)) +
        guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
        title.position = "top", title.hjust = 0.5))+
        geom_text(mapping=aes(x=tag1, y=tag2, label=freq), color='white', fontface='bold', size=6) +
        labs(y=' ', x='Overlap')
    outpdf('naitive_generated', 'seq_overlap_heatmap', height=15,width=20)

}



plot_pairwise_overlap_native_generated_boxplot = function(){
    infile = 'eleven_outfiles/pairwise_overlap_native_generated.csv'
    df = read_csv(infile)
    print(df)
    df = df[df$tag1!=df$tag2,]
    gdf = df[(grepl('generated',df$tag2)) & (grepl('generated',df$tag1)),]
    gdf$origin_tag = 'generated_generated'
    ndf = df[(grepl('native',df$tag2)) & (grepl('generated',df$tag1)),]
    ndf$origin_tag = 'native_generated'
    mdf = rbind(gdf, ndf)
    mdf$freq2 = mdf$freq2*100
    mdf$origin_tag[mdf$origin_tag=='native_generated'] = 'Native vs generated'
    mdf$origin_tag[mdf$origin_tag=='generated_generated'] = 'Generated vs generated'
    mdf$origin_tag = factor(mdf$origin_tag, levels=c('Native vs generated', 'Generated vs generated'))
    meddf = mdf %>% group_by(origin_tag) %>% summarize(medval=median(freq2))
    meddf$medval = sprintf('Median: %s', meddf$medval)
    ggplot(data=mdf, mapping=aes(x=origin_tag, y=freq2)) +
        geom_boxplot(mapping=aes(color=origin_tag), outlier.shape=NA) +
        geom_point(mapping=aes(color=origin_tag),size = 7, shape = 21, position = position_jitterdodge())+
        geom_text(data=meddf, mapping=aes(x=origin_tag,y=0, label=medval), size=7) +
        labs(y='CDR-H3 overlap (%)', x=' ') +
        theme(legend.title=element_blank(),
                axis.title = element_text(size=25),
                legend.text= element_text(size=25),
                axis.text.x = element_blank(),
                axis.text.y=element_text(size=25)) +
        scale_color_manual(values=c(nativecol,igemcol))
    outpng('native_generated', 'seq_overlap_boxplot', width=9)

}

ld_count_cor = function () {
    infile = 'eleven_outfiles/antigens_native_generated_ldcount_merged_ldmerged.csv'
    df = read_csv(infile)
    gdf = df[df$origin=='generated_generated',]
    ndf = df[df$origin=='native_native',]
    print(dim(gdf))
    print(dim(ndf))
    antigens = unique(df$antigen)
    print(antigens)
    outdf = tibble()
    for (antigen in antigens){
        for (antigen2 in antigens){
            andf=ndf[(ndf$antigen==antigen) & (ndf$antigen2==antigen2),]
            agdf = gdf[(gdf$antigen==antigen) & (gdf$antigen2==antigen2),]
            mdf = merge(agdf, andf, by='ld')
            corval = cor(mdf$count.x, mdf$count.y)
            print(corval)
            cordf = data_frame('antigen'=antigen, 'antigen2'=antigen2, 'corval'=corval)
            outdf = rbind(outdf, cordf)

        }
    }
    outdf$xtag = sprintf('%s_%s', outdf$antigen, outdf$antigen2)
    ggplot(data=outdf)+
        geom_bar(mapping=aes(x=xtag, y=corval), stat='identity', position='dodge2') +
        facet_wrap(~antigen2,nrow=2, scale='free') +
        labs(x='Antigen_pair', y='Correleation of LD distances between generated and native') +
        theme(axis.text.x = element_text(angle=90))
    outpng('native_generated', 'ldfcount', width=11)
}

ld_between_native_generated_distr = function(){
    infile = 'eleven_outfiles/antigens_native_generated_ldcount_merged_ldmerged.csv'
    df = read_csv(infile)
    df = df[df$origin=='generated_native',]
    df = df[df$antigen==df$antigen2,]
    df %>% print(n=100)
    ggplot(data=df) +
        geom_bar(mapping=aes(x=ld, y=count), stat='identity', position='dodge2') +
        facet_wrap(~antigen)
    outpng('tmp', 'tmp')
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

ld_native_generated_distr = function(){
    infile = 'eleven_outfiles/antigens_native_generated_ldcount_merged_ldmerged.csv'
    df = read_csv(infile)
    df = fixfive2(df)
    df = df[df$origin!='generated_native',]
    df = df[df$antigen==df$antigen2,]
    print(df)
    ggplot(data=df) +
        geom_bar(mapping=aes(x=ld, y=count, fill=origin), stat='identity', position='dodge2') +
        facet_wrap(~antigen, scale='free')
    outpng('tmp', 'tmp')

}


antigen_position_distr = function(){
    ginfile = 'eleven_outfiles/antigens_feature_generated_native_sized.csv'
    ninfile = 'eleven_outfiles/antigens_feature_native.csv'
    gdf = read_csv(ginfile)
    ndf = read_csv(ninfile)
    df = rbind(gdf, ndf)
    df = df[df$antigen!='2DD8', ]
    df = fixfive(df)
    df$origin = factor(df$origin, levels = c('native', 'generated'))
    # df = df %>% group_by()
    df = df %>% group_by(antigen, AGbindPositions) %>%mutate(countval = n())
    nepidf = df %>% group_by(antigen, origin) %>% summarize(nepival=length(unique(AGbindPositions)))
    print(nepidf)
    nepidf$xpos = 10
    nepidf$ypos = 40000
    nepidf$ypos[nepidf$origin=='generated'] = 20000
    nepidf$nepilabel = sprintf('Total: %s', nepidf$nepival)
    ggplot(data=df) +
        geom_bar(mapping= aes(x=reorder(AGbindPositions, -countval), fill=origin)) +
        facet_wrap(~antigen, nrow=2, scale='free') +
        geom_text(data=nepidf, mapping=aes(x=xpos,y=ypos, label=nepilabel, color=origin), fontface='bold',
                    size=7) +
        theme(axis.text.x = element_text(angle=90, size=4),
            legend.title = element_blank(),
            axis.title = element_text(size=40),
            strip.text= element_text(size=20),
            legend.text = element_text(size=30),
            axis.text.y=element_text(size=20)) +
        scale_fill_manual(values=c(nativecol, igemcol)) +
        scale_color_manual(values=c(nativecol, igemcol)) +
        labs(y='# of epitopes', x = 'Epitopes')
    outpng('native_generated', 'antigenbindingpos_distr', width=17, heigh=12)
}



ld_between_native_generated_distr = function(){
    # infile = 'eleven_outfiles/antigens_native_generated_ldcount_merged_ldmerged.csv'
    infile = 'eleven_outfiles/antigens_native_generated_ldcount_rsample_merged.csv'
    df = read_csv(infile)
    df = fixfive2(df)
    df = df[df$origin=='generated_native',]
    df = df[df$antigen==df$antigen2,]
    df %>% print(n=100)
    freq2med = function(ld, freq){
        d2 = rep(ld, freq)
        med = median(d2)
    }
    meadf = df %>% group_by(antigen, antigen2, ld) %>% summarize(mean_nld=mean(nld),
                                                            max_nld = max(nld),sd_nld=sd(nld),
                                                            mean_ld=mean(ld))
    meddf = meadf %>% group_by(antigen, antigen2) %>% summarize(medval=freq2med(ld,mean_nld),maxval=max(mean_nld))
    print(meddf)
    # stop()
    # meddf2 = meddf %>% group_by(antigen, antigen2) %>% summarize(medval=median(medval),sdval = sd(medval))
    # print(meddf2)
    # stop()
    print(meadf)
    meddf$medlabel = sprintf('Median: %s', meddf$medval)
    ggplot(data=meadf) +
        geom_bar(mapping=aes(x=ld, y=mean_nld), stat='identity', position='dodge2') +
        geom_errorbar(mapping=aes(x=ld, y=mean_nld,ymin=mean_nld-sd_nld,ymax=mean_nld+sd_nld)) +
        geom_text(data=meddf, mapping=aes(x=0, y=maxval*1.1, label=medlabel), hjust=0, size=7) +
        facet_wrap(~antigen, scale = 'free',  nrow=2) +
        labs(y='# of CDR-H3 sequences', x='Levenshtein distance (LD)') +
        theme(axis.title=element_text(size=25),
                strip.text=element_text(size=25),
                axis.text = element_text(size=20))
    outpng('ld_between', 'native_generated', width = 17, height=10)
}

# run stuff

# get_top_bot_overlap('outfiles/antigens_geminio_absolute_rank_develparams.csv', 'top')
# get_top_bot_overlap('outfiles/antigens_native_delicab_rank_develparams.csv', 'top')
# get_native_igem_kmercount_cor_m()
# get_native_igem_kmercount_cor()
# binding_affinity_igem_native_dist_grid()
# develparams_igem_native_dist_grid()
# show_colors()
#get_pairwise_overlap()
# do_jensen_shannon()
# ld_dist_cutoff()


## new panels
binding_affinity_native_generated() # panel A
# binding_fold_cor() # panel B top
# binding_fold_bar() # panel B top
# kmer_native_generated_cor() # panel C
# devel_params_native_generated() # panel F
# pairwise_overlap_native_generated() # panel E old
# get_pairwise_overlap_native_generated_boxplot() # data for panel E
# plot_pairwise_overlap_native_generated_boxplot() # panel E
# ld_count_cor() #
# ld_between_native_generated_distr() #panel D
# ld_native_generated_distr()
# antigen_position_distr() # panel B bottom