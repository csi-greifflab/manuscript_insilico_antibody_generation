# preprocess data for mainfig2

# import stuff
import pandas as pd
import sys


def get_binding_fold_counts():
    '''
    gets binding fold counts from absolut's output, native and generated
    :return:
    '''
    ginfile = 'eleven_outfiles/antigens_affinity_generated_native_sized.csv'
    ninfile = 'eleven_outfiles/antigens_affinity_native.csv'
    gdf = pd.read_csv(ginfile)
    ndf = pd.read_csv(ninfile)
    outcontents = []
    for antigen in gdf.antigen.unique():
        agdf = gdf[gdf.antigen==antigen]
        andf = ndf[ndf.antigen==antigen]
        folds = set(agdf.Structure.tolist()  + andf.Structure.tolist())
        # print(folds, len(folds))
        print(len(agdf.Structure.unique()))
        print(len(andf.Structure.unique()))
        print(antigen)
        for fold in folds:
            fagdf = agdf[agdf.Structure == fold]
            fandf = andf[andf.Structure == fold]
            # fagcount='NA'
            fagcount=0
            if fagdf.shape[0] != 0:
                fagcount= fagdf.shape[0]
            # fancount = 'NA'
            fancount = 0
            if fandf.shape[0] != 0:
                fancount=fandf.shape[0]
            # fagcontent = [antigen, fold, 'generated', fagcount]
            # fancontent = [antigen, fold, 'native', fancount]
            fangcontent = [antigen, fold, fancount, fagcount]
            # outcontents.append(fagcontent)
            # outcontents.append(fancontent)
            outcontents.append(fangcontent)
    # colnames = ['antigen', 'Structure', 'origin', 'structure_count']
    colnames = ['antigen', 'Structure', 'native_count', 'generated_count']
    outdf = pd.DataFrame(outcontents, columns=colnames)
    outname = 'eleven_outfiles/antigens_native_generated_structure_count.csv'
    outdf.to_csv(outname, index=False)

def check_ld_data():
    '''
    qcs the ldcount data
    :return:
    '''
    infile = 'eleven_outfiles/antigens_native_generated_ldcount_merged.csv'
    df = pd.read_csv(infile)
    print(df)
    adf = df[(df.antigen=='3RAJ') & (df.antigen2=='3RAJ')]
    agdf = adf[df.origin =='generated_generated']
    andf = adf[df.origin=='native_native']
    outname1 = 'eleven_outfiles/tmp1.csv'
    outname2 = 'eleven_outfiles/tmp2.csv'
    agdf.to_csv(outname1, index=False)
    andf.to_csv(outname2, index=False)


# runs tuff

get_binding_fold_counts()
# check_ld_data()



