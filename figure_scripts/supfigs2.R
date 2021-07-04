# Title     : TODO
# Objective : TODO
# Created by: rahmadakbar
# Created on: 2020-09-02


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
library(provenance)
library(distrEx)
library(transport)
library(kSamples)


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


affinity_nsamples_transfer_ks= function(){
    generated_infile = 'eleven_outfiles/antigens_affinity_ntrainingsamples_generated.csv'
    transfer_infile = 'eleven_outfiles/antigens_affinity_ntrainingsamples_transfer_generated.csv'
    native_infile = 'eleven_outfiles/antigens_affinity_native.csv'
    gdf = read_csv(generated_infile)
    tdf = read_csv(transfer_infile)
    ndf = read_csv(native_infile)
    ndf$model = rep('native', dim(ndf)[1])
    df=rbind(gdf,tdf, ndf)
    df$model2 = sprintf('%s_%s', df$antigen, df$model)
    models = unique(df$model2)
    # models = models[!grepl('native',models)]
    outdfs = data.frame()
    for (model in models){
        mdf = df[df$model2==model,]
        antigen = unique(mdf$antigen)
        origin = unique(mdf$origin)
        native_antigendf = df[df$model2==sprintf('%s_native', antigen),]
        k = ks.test(mdf$Energy, native_antigendf$Energy)
        ks=k$statistic
        # print(k[['intrinsic.discrepancy']])
        # print(model)
        # print(mdf$model)
        model_parts  = strsplit(mdf$model[1], '_')[[1]]
        model_short = sprintf('%s_%s', model_parts[1], model_parts[3])
        if (grepl('impoverished', model)){
            model_short=sprintf('-T_ntrain%s', model_parts[3])
            crossstatus = '-T'
        }
        else if (grepl('interposed', model)){
            model_short=sprintf('+T_ntrain%s', model_parts[3])
            crossstatus = '+T'
        }
        # print(model)
        else if (grepl('opulent', model) | grepl('native', model)){
            model_short = sprintf('-T_%s',mdf$model[1])
            crossstatus = '-T'
        }
        print(model_short)
        outdf = data.frame('antigen' = antigen,
                            'origin' = origin,
                            'model'=model,
                            'model_short'=model_short,
                            'KSD'= ks,
                            'crossstatus'=crossstatus)
        outdfs = rbind(outdfs, outdf)
        # print(outdfs)
        # if (grepl('opulent', model) | grepl('native', model)){
        #     outdfs=rbind(outdfs,outdf)
        #     }
    }
    print(outdfs)
    # outdfs$model_short = factor(outdfs$model_short, levels=c('impoverished_700',
    # 'interposed_700', 'impoverished_7000', 'interposed_7000', 'opulent', 'native'))
    outname = 'eleven_outfiles/affinity_nsamples_transfer_ks.csv'
    write.csv(outdfs, outname, row.names=FALSE)
    # ggplot(data=outdfs) +
    # geom_bar(mapping = aes(x=model, y=KSD), stat='identity', position = 'dodge2') +
    #     geom_boxplot(mapping = aes(x=model_short, y=KSD)) +
    #     geom_point(mapping = aes(x=model_short, y=KSD), poistion=position_jitterdodge()) +
    #     facet_wrap(~antigen, nrow=2) +
    #     labs(x='TransferType_NumberOfTrainingSequences', y='KSD against native') +
    #     theme(axis.text.x = element_text(angle=90, size=25),
    #             axis.title = element_text(size=30))
    # outpng('native', 'nsamples_transfererd_KSD', width=17)
    # df$model = factor(df$model, levels=c('native', 'opulent', sprintf('impoverished_sample%s_700', c(0:4)),
    # sprintf('impoverished_sample%s_7000',c(0:4)),
    # sprintf('interposed_sample%s_700',c(0:4)),
    # sprintf('interposed_sample%s_7000',c(0:4))))
    # print(unique(df$model))
}


plot_crosstransfer = function(){
    infile = 'eleven_outfiles/crosstransfer_native_generated_nsamples_opulent_sampled.csv'
    df = read_csv(infile)
    df$crossantigen_tag = sprintf('%s_%s_%s_%s', df$crossantigen, df$sample, df$ntrainseq,df$crossstatus)
    print(df)
    antigens = unique(df$antigen)
    samples = unique(df$sample)
    print(antigens)
    print(samples)
    outdfs = data.frame()
    for (antigen in antigens){
        ndf = df[(df$antigen==antigen) & (df$sample=='native'),]
        ctags = unique(df[df$antigen==antigen,]$crossantigen_tag)
        for (ctag in ctags){
            sdf = df[df$crossantigen_tag==ctag,]
            # k = ks.test(ndf$Energy,sdf$Energy)
            # ksd=k$statistic
            # print(sprintf('ks.test: %s', ksd))
            ksd2 = KS.diss(ndf$Energy, sdf$Energy)
            # tvd = TotalVarDist(ndf$Energy, sdf$Energy)
            # wd = wasserstein(ndf$Energy, sdf$Energy)
            # print(sprintf('KS.diss: %s, tvd: %s', ksd2, wd))
            # print(unique(sdf$crossstatus))
            origin = unique(sdf$origin)
            crossstatus = unique(sdf$crossstatus)
            ntrainseq = unique(sdf$ntrainseq)
            sample = unique(sdf$sample)
            if (grepl('sample', sample)){sample='sample'}
            # print(sample)
            outdf = data.frame('antigen'=antigen, 'origin'=origin, 'crossstatus'=crossstatus,
            'ctag'=ctag,'KSD'=ksd2, 'ntrainseq'=ntrainseq, 'sample'=sample)
            outdfs = rbind(outdfs, outdf)
        }
    }
    outdfs=as_tibble(outdfs)
    print(outdfs)
    stop(0)
    transfer_type = c('nocross' = 'NO', 'cross' = 'Biased')
    outdfs$crosstatus[outdfs$crossstatus=='cross'] = 'Biased'
    outdfs$xtag= sprintf('%s_%s_%s', outdfs$sample,outdfs$crossstatus, outdfs$ntrainseq)
    print(outdfs$xtag)
    outdfs$xtagfact = factor(outdfs$xtag, levels = c('sample_nocross_700', 'sample_cross_700',
    'sample_nocross_7000', 'sample_cross_7000',
    'opulent_nocross_70000','native_nocross_0'))
    print(outdfs)
    outdfs$antigen = as.character(outdfs$antigen)
    print(outdfs)
    ggplot(data=outdfs) +
    # geom_density_ridges(mapping=aes(x=Energy, y=crossantigen_tag, fill=crossstatus)) +
        geom_boxplot(mapping=aes(x=xtagfact, y=KSD)) +
        geom_point(mapping = aes(x=xtagfact, y=KSD), poistion=position_jitterdodge()) +
        facet_wrap(~antigen, nrow=2) +
        labs(x='TransferType_numberOfTrainingSequences', y='KSD against native') +
        theme(axis.text.x=element_text(angle=90,size=25),
        axis.title = element_text(size=30))

    outpng('native_generated', 'crosstransfer',  width=17)

}


affinity_nsamples = function(){
    generated_infile = 'eleven_outfiles/antigens_affinity_ntrainingsamples_generated.csv'
    # transfer_infile = 'eleven_outfiles/antigens_affinity_ntrainingsamples_transfer_generated.csv'
    native_infile = 'eleven_outfiles/antigens_affinity_native.csv'
    gdf = read_csv(generated_infile)
    # tdf = read_csv(transfer_infile)
    # print(gdf)
    # change_model_name = function(x){
    #     if (grepl('sample', x)){
    #         return(sprintf('impoverished_%s', x))
    #     }
    #     else{return(sprintf('opulent'))}
    # }
    # model2 = sapply(gdf$model, function(x) change_model_name(x))
    # print(gdf)
    # gdf$model = model2
    ndf = read_csv(native_infile)
    ndf$model = rep('native', dim(ndf)[1])
    df=rbind(gdf, ndf)
    print(unique(df$model))
    # stop()
    df$model = factor(df$model, levels=c('native', 'opulent', sprintf('impoverished_sample%s_700', c(0:4)),
    sprintf('impoverished_sample%s_7000',c(0:4)),
    sprintf('interposed_sample%s_700',c(0:4)),
    sprintf('interposed_sample%s_7000',c(0:4))))
    print(unique(df$model))
    df$origin = factor(df$origin, levels=c('native', 'generated'))
    ggplot(data=df) +
        geom_density_ridges(mapping=aes(x=Energy, y=model, color=origin)) +
        facet_wrap(~antigen, nrow=2) +
        labs(y='Model_NumberOfTrainingSequences') +
        theme(axis.title = element_text(size=25)) +
        scale_color_manual(values=c(nativecol,igemcol))
    outpng('affinity','ntrainingsamples', width=15)
}


plot_crosstransfer_ad = function(){
    # anderson darling
    infile = 'eleven_outfiles/crosstransfer_native_generated_nsamples_opulent_sampled.csv'
    df = read_csv(infile)
    df$crossantigen_tag = sprintf('%s_%s_%s_%s', df$crossantigen, df$sample, df$ntrainseq,df$crossstatus)
    print(df)
    antigens = unique(df$antigen)
    samples = unique(df$sample)
    print(antigens)
    print(samples)
    outdfs = data.frame()
    for (antigen in antigens){
        ndf = df[(df$antigen==antigen) & (df$sample=='native'),]
        ctags = unique(df[df$antigen==antigen,]$crossantigen_tag)
        for (ctag in ctags){
            sdf = df[df$crossantigen_tag==ctag,]
            # k = ks.test(ndf$Energy,sdf$Energy)
            # ksd=k$statistic
            # print(sprintf('ks.test: %s', ksd))
            # ksd2 = KS.diss(ndf$Energy, sdf$Energy)
            adouts = ad.test(ndf$Energy, sdf$Energy)
            adc = adouts$ad[1,1]
            # tvd = TotalVarDist(ndf$Energy, sdf$Energy)
            # wd = wasserstein(ndf$Energy, sdf$Energy)
            # print(sprintf('KS.diss: %s, tvd: %s', ksd2, wd))
            # print(unique(sdf$crossstatus))
            origin = unique(sdf$origin)
            crossstatus = unique(sdf$crossstatus)
            ntrainseq = unique(sdf$ntrainseq)
            sample = unique(sdf$sample)
            if (grepl('sample', sample)){sample='sample'}
            # print(sample)
            outdf = data.frame('antigen'=antigen, 'origin'=origin, 'crossstatus'=crossstatus,
            'ctag'=ctag,'KSD'=adc, 'ntrainseq'=ntrainseq, 'sample'=sample)
            outdfs = rbind(outdfs, outdf)
        }
    }
    outdfs=as_tibble(outdfs)
    transfer_type = c('nocross' = 'NO', 'cross' = 'Biased')
    outdfs$crosstatus[outdfs$crossstatus=='cross'] = 'Biased'
    outdfs$xtag= sprintf('%s_%s_%s', outdfs$sample,outdfs$crossstatus, outdfs$ntrainseq)
    print(outdfs$xtag)
    outdfs$xtagfact = factor(outdfs$xtag, levels = c('sample_nocross_700', 'sample_cross_700',
    'sample_nocross_7000', 'sample_cross_7000',
    'opulent_nocross_70000','native_nocross_0'))
    print(outdfs)
    outdfs$antigen = as.character(outdfs$antigen)
    print(outdfs)
    ggplot(data=outdfs) +
    # geom_density_ridges(mapping=aes(x=Energy, y=crossantigen_tag, fill=crossstatus)) +
        geom_boxplot(mapping=aes(x=xtagfact, y=KSD)) +
        geom_point(mapping = aes(x=xtagfact, y=KSD), poistion=position_jitterdodge()) +
        facet_wrap(~antigen, nrow=2) +
        labs(x='Models', y='AD criterion against native') +
        theme(axis.text.x=element_text(angle=90,size=25),
        axis.title = element_text(size=30))

    outpng('native_generated', 'crosstransfer',  width=17)

}



affinity_nsamples_crosstransfer_ks = function(){
    infile = 'eleven_outfiles/crosstransfer_native_generated_nsamples_opulent_sampled.csv'
    df = read_csv(infile)
    df$crossantigen_tag = sprintf('%s_%s_%s_%s', df$crossantigen, df$sample, df$ntrainseq,df$crossstatus)
    print(df)
    antigens = unique(df$antigen)
    samples = unique(df$sample)
    print(antigens)
    print(samples)
    outdfs = data.frame()
    for (antigen in antigens){
        ndf = df[(df$antigen==antigen) & (df$sample=='native'),]
        ctags = unique(df[df$antigen==antigen,]$crossantigen_tag)
        for (ctag in ctags){
            sdf = df[df$crossantigen_tag==ctag,]
            # k = ks.test(ndf$Energy,sdf$Energy)
            # ksd=k$statistic
            # print(sprintf('ks.test: %s', ksd))
            ksd2 = KS.diss(ndf$Energy, sdf$Energy)
            # tvd = TotalVarDist(ndf$Energy, sdf$Energy)
            # wd = wasserstein(ndf$Energy, sdf$Energy)
            # print(sprintf('KS.diss: %s, tvd: %s', ksd2, wd))
            # print(unique(sdf$crossstatus))
            origin = unique(sdf$origin)
            crossstatus = unique(sdf$crossstatus)
            ntrainseq = unique(sdf$ntrainseq)
            sample = unique(sdf$sample)
            if (grepl('sample', sample)){sample='sample'}
            # print(sample)
            outdf = data.frame('antigen'=antigen, 'origin'=origin, 'crossstatus'=crossstatus,
            'ctag'=ctag,'KSD'=ksd2, 'ntrainseq'=ntrainseq, 'sample'=sample)
            outdfs = rbind(outdfs, outdf)
        }
    }
    outdfs=as_tibble(outdfs)
    print(outdfs)
    # transfer_type = c('nocross' = 'NO', 'cross' = 'Biased')
    # outdfs$crosstatus[outdfs$crossstatus=='cross'] = 'Biased'
    # outdfs$xtag= sprintf('%s_%s_%s', outdfs$sample,outdfs$crossstatus, outdfs$ntrainseq)
    # print(outdfs$xtag)
    # outdfs$xtagfact = factor(outdfs$xtag, levels = c('sample_nocross_700', 'sample_cross_700',
    # 'sample_nocross_7000', 'sample_cross_7000',
    # 'opulent_nocross_70000','native_nocross_0'))
    # print(outdfs)
    # outdfs$antigen = as.character(outdfs$antigen)
    # print(outdfs)
    outname = 'eleven_outfiles/affinity_nsamples_crosstransfer_ks.csv'
    write.csv(outdfs, outname, row.names=FALSE)
}


fixfive =function(df){
    # fix 5E94
    df$antigen[df$antigen=='FiveE94'] = '5E94'
    return(df)
}

plot_affinity_nsamples_crosstransfer_ks = function(){
    infile = 'eleven_outfiles/affinity_nsamples_crosstransfer_ks.csv'
    df = read_csv(infile)
    df$sample[df$sample == 'opulent'] = 'datarich'
    df$crossstatus[df$crossstatus=='cross'] = '+T'
    df$crossstatus[df$crossstatus=='nocross'] = '-T'
    get_xtag = function(xtag){
        if(grepl('sample', xtag)){
            parts = strsplit(xtag, '_')[[1]]
            newtag = sprintf('%s_ntrain%s', parts[2], parts[3])
        }
        else{
            parts = strsplit(xtag, '_')[[1]]
            print(parts)
            newtag = sprintf('%s_%s', parts[2], parts[1])
        }
        return(newtag)
    }
    df$xtag = sprintf('%s_%s_%s', df$sample,df$crossstatus, df$ntrainseq)
    df$xtag = sapply(df$xtag, get_xtag)
    print(unique(df$xtag))
    print(df)
    df = fixfive(df)
    df$xtag= factor(df$xtag, levels=c('-T_ntrain700', '+T_ntrain700',
                                        '-T_ntrain7000', '+T_ntrain7000',
                                        '-T_datarich', '-T_native'))
    meddf = df %>% group_by(antigen, ntrainseq, crossstatus, xtag) %>% summarize(medval=median(KSD))
    meddf$medval = round(meddf$medval,2)
    print(meddf)
    ggplot(data=df, mapping=aes(x=xtag, y=KSD)) +
        geom_boxplot(mapping=aes(color=crossstatus)) +
        geom_point(aes(color = crossstatus), size = 5, shape = 21, position = position_jitterdodge()) +
        geom_text(data=meddf, mapping=aes(x=xtag, y=1., label=medval, color=crossstatus)) +
        facet_wrap(~antigen, nrow=2) +
        theme(axis.text.x=element_text(angle=90, size=15),
            # legend.title=element_blank(),
            legend.position='NA',
            axis.title=element_text(size=25)) +
        scale_color_manual(values=c('orange', 'black')) +
        labs(x='', y='KSD generated against native (Affinity)')
    outpng('native_generated', 'crosstransfer', width=15)
}

plot_affinity_nsamples_transfer_ks = function(){
    infile = 'eleven_outfiles/affinity_nsamples_transfer_ks.csv'
    df = read_csv(infile)
    print(df)
    df$model_short[df$model_short == '-T_opulent'] = '-T_datarich'
    df = fixfive(df)
    df$xtag= factor(df$model_short, levels=c('-T_ntrain700', '+T_ntrain700',
                                '-T_ntrain7000', '+T_ntrain7000',
                                '-T_datarich', '-T_native'))
    meddf = df %>% group_by(antigen, model_short,xtag, crossstatus) %>% summarize(medval=median(KSD))
    meddf$medval = round(meddf$medval,2)
    ggplot(data=df, mapping=aes(x=xtag, y=KSD)) +
        geom_boxplot(mapping=aes(color=crossstatus)) +
        geom_point(aes(color = crossstatus), size = 5, shape = 21, position = position_jitterdodge()) +
        geom_text(data=meddf, mapping=aes(x=xtag, y=1., label=medval, color=crossstatus)) +
        facet_wrap(~antigen, nrow=2) +
        theme(axis.text.x=element_text(angle=90, size=15),
        # legend.title=element_blank(),
        legend.position='NA',
        axis.title=element_text(size=25)) +
        scale_color_manual(values=c(my_base_color[2], 'black')) +
        labs(x='', y='KSD against native')
    outpng('native_generated', 'transfer', width=15)
}



affinity_nsamples2 = function(){
    generated_infile = 'eleven_outfiles/antigens_affinity_ntrainingsamples_ntrainingsamples2_generated_sampled.csv'
    # transfer_infile = 'eleven_outfiles/antigens_affinity_ntrainingsamples_transfer_generated.csv'
    native_infile = 'eleven_outfiles/antigens_affinity_native_sampled.csv'
    gdf = read_csv(generated_infile)
    # tdf = read_csv(transfer_infile)
    # print(gdf)
    # change_model_name = function(x){
    #     if (grepl('sample', x)){
    #         return(sprintf('impoverished_%s', x))
    #     }
    #     else{return(sprintf('opulent'))}
    # }
    # model2 = sapply(gdf$model, function(x) change_model_name(x))
    # print(gdf)
    # gdf$model = model2
    ndf = read_csv(native_infile)
    ndf$model = rep('native', dim(ndf)[1])
    df=rbind(gdf, ndf)
    print(unique(df$model))
    # stop()
    df$model = factor(df$model, levels=c(sprintf('impoverished_sample%s_700', c(0:4)),
            sprintf('impoverished_sample%s_7000',c(0:4)),
            sprintf('impoverished_sample%s_10000',c(0:4)),
            sprintf('impoverished_sample%s_20000',c(0:4)),
            sprintf('impoverished_sample%s_30000',c(0:4)),
            sprintf('impoverished_sample%s_40000',c(0:4)),
            sprintf('impoverished_sample%s_50000',c(0:4)),
            sprintf('impoverished_sample%s_60000',c(0:4)),
            'opulent',
            'native'))
    print(unique(df$model))
    print(df)
    meddf = df %>% group_by(antigen,model,origin) %>% summarize(medval=median(Energy), maxval=max(Energy))
    print(meddf)
    meddf$medval = round(meddf$medval,2)
    get_ntrain = function(x){
        print(x)
        if(grepl('impoverished', x)){
            parts = strsplit(x, '_')[[1]]
            model_name=sprintf('ntrain%s', parts[3])
            return(model_name)
        }
        else{return(x)}
    }
    get_ntrain2 = function(x){
        print(x)
        if(grepl('impoverished', x)){
            parts = strsplit(x, '_')[[1]]
            model_name=sprintf('%s_ntrain%s', parts[2],parts[3])
            return(model_name)
        }
        else{return(x)}
    }
    meddf$model_name = sapply(as.character(meddf$model), get_ntrain)
    print(meddf)
    meddf = fixfive(meddf)
    meddf$model_name[meddf$model_name == 'opulent'] = 'ntrain70000'
    meddf$model_name = factor(meddf$model_name, levels=c('ntrain700',
                                                        'ntrain7000',
                                                        'ntrain10000',
                                                        'ntrain20000',
                                                        'ntrain30000',
                                                        'ntrain40000',
                                                        'ntrain50000',
                                                        'ntrain60000',
                                                        'ntrain70000',
                                                        'native'))
    meddf$origin = factor(meddf$origin, levels=c('native', 'generated'))
    labdf = meddf%>% group_by(antigen, origin, model_name) %>% summarize(medval = median(medval), maxval=min(maxval))
    labdf = labdf%>% group_by(antigen) %>% mutate(maxval = max(maxval))
    labdf$medval = round(labdf$medval,0)
    print(labdf)
    # ggplot(data=meddf, mapping=aes(x=model_name, y=medval)) +
    #     geom_boxplot(mapping=aes(color=origin)) +
    #     geom_point(aes(color = origin), size = 5, shape = 21, position = position_jitterdodge()) +
    #     geom_text(data=labdf, mapping= aes(x=model_name, y=maxval*1.13, label=medval, color=origin), size=3) +
    #     facet_wrap(~antigen, nrow=2, scale='free_y') +
    #     labs(y='Affinity (Median Energy)', x='') +
    #     theme(axis.title = element_text(size=25),
    #             axis.text.x = element_text(angle=90),
    #             legend.title = element_blank(),
    #             legend.text = element_text(size=20)) +
    #     scale_color_manual(values=c(nativecol,igemcol)) +
    #     scale_fill_manual(values = c(nativecol, igemcol))
    # outpng('affinity','ntrainingsamples_medval_boxplot', width=16, height=7):
    # stop()
    df$origin = factor(df$origin, levels=c('native', 'generated'))
    df$model_name = sapply(as.character(df$model), get_ntrain2)
    df$model_name[df$model_name == 'opulent'] = 'ntrain70000'
    df$model_name = factor(df$model_name, levels=c(sprintf('sample%s_ntrain700', c(0:4)),
    sprintf('sample%s_ntrain7000',c(0:4)),
    sprintf('sample%s_ntrain10000',c(0:4)),
    sprintf('sample%s_ntrain20000',c(0:4)),
    sprintf('sample%s_ntrain30000',c(0:4)),
    sprintf('sample%s_ntrain40000',c(0:4)),
    sprintf('sample%s_ntrain50000',c(0:4)),
    sprintf('sample%s_ntrain60000',c(0:4)),
    'ntrain70000',
    'native'))
    print(unique(df$model_name))
    df = fixfive(df)
    print(df)
    meddf = df %>% group_by(antigen, model_name) %>% summarize(medval=median(Energy))
    print(meddf)
    meddf$medval = round(meddf$medval, 3)
    ggplot(data=df) +
        geom_density_ridges(mapping=aes(x=Energy, y=model_name, color=origin, fill=origin)) +
        geom_text(data=meddf, mapping=aes(x=medval, y=model_name, label=medval)) +
        facet_wrap(~antigen, nrow=2) +
        labs(y='Model', x='Affinity (Energy)') +
        theme(axis.title = element_text(size=25),
                legend.title = element_blank(),
                legend.text = element_text(size=25),
                axis.text.x=element_text(size=20),
                strip.text = element_text(size=20)) +
        scale_color_manual(values=c(nativecol,igemcol)) +
        scale_fill_manual(values = c(nativecol, igemcol))
    outpng('fig','S2', width=25, height=15)
}


# run stuff
affinity_nsamples2()
