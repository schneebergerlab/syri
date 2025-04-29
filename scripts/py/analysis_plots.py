def get_snp_indel_plots():
    '''
    Read the VCFdist summary file and generate summary plot for the SNPs/Indels identified using Syri, deepvariant, and
    gatk. Statistics considered:
        1) True Positive counts
        2) Precision score
        3) Recall score
        4) F1_QScore
    :return:
    '''
    import pandas as pd
    from matplotlib import pyplot as plt
    from hometools.plot import cleanax
    def getdata(p, t, v):
        df = pd.read_table(p)
        df['tool'] = t
        df['variant'] = v
        return df



    indir = "/home/ra98jam/d16/projects/sv_benchmark/results/variant_calls/"
    tools = ["syri", "deepvariant", "gatk"]

    syrisnp = getdata(f'{indir}/syri/vcfdist/precision-recall-summary.tsv', 'syri', 'snps')
    syrinofiltsnp = getdata(f'{indir}/syri/vcfdist_nofilt/precision-recall-summary.tsv', 'syri_nofilt', 'snps')
    dvsnp = getdata(f'{indir}/deepvariant/vcfdist/precision-recall-summary.tsv', 'deepvariant', 'snps')
    gatksnp = getdata(f'{indir}/gatk/vcfdist/precision-recall-summary.tsv', 'gatk', 'snps')
    snpdf = pd.concat([syrisnp, syrinofiltsnp, dvsnp, gatksnp])
    snpdf = snpdf.loc[snpdf['VAR_TYPE'] == 'SNP']
    snpdf = snpdf.loc[snpdf['THRESHOLD'] == 'NONE']


    syriindels = getdata(f'{indir}/syri/vcfdist_indels/precision-recall-summary.tsv', 'syri', 'indels')
    syrinofiltindels = getdata(f'{indir}/syri/vcfdist_nofilt_indels/precision-recall-summary.tsv', 'syri_nofilt', 'indels')
    dvindels = getdata(f'{indir}/deepvariant/vcfdist_indels/precision-recall-summary.tsv', 'deepvariant', 'indels')
    gatkindels = getdata(f'{indir}/gatk/vcfdist_indels/precision-recall-summary.tsv', 'gatk', 'indels')
    indelsdf = pd.concat([syriindels, syrinofiltindels, dvindels, gatkindels])
    indelsdf = indelsdf.loc[indelsdf['VAR_TYPE'] == 'INDEL']
    indelsdf = indelsdf.loc[indelsdf['THRESHOLD'] == 'NONE']
    bar_colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red']
    bar_labels = ['syri', 'syri_nofilt', 'deepvariant', 'gatk']
    bar_ylim = [1.05, 1.05, 30]
    x_stats = 'TRUTH_TP  TRUTH_FN QUERY_TP  QUERY_FP  PREC     RECALL    F1_QSCORE'.split()
    
    
    def get_stats_fig(fig, df, title):
        for j, stat in enumerate(x_stats):
            ax = fig.add_subplot(2, 4, j + 1)
            temp = ax.bar(df['tool'], df[stat].round(4), color=bar_colors, label=bar_labels, alpha=0.75,
                          edgecolor='black',
                          width=0.5)
            ax.set_ylabel(stat)
            ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation=45, ha='right')
            ax = cleanax(ax)
            plt.suptitle(title)
            plt.tight_layout()
        return fig

    fig = plt.figure(figsize=(12, 8))
    fig = get_stats_fig(fig, snpdf, 'SNVs')
    # for j, stat in enumerate('PREC    RECALL   F1_QSCORE'.split()):
    
        ax.bar_label(temp)
        ax.set_ylim([0, bar_ylim[j]])
        if stat == 'RECALL':
            ax.set_title('SNVs' if i == 0 else 'Indels', fontsize=14)
    for i, df in enumerate([snpdf, indelsdf]):

    # plt.legend(bbox_to_anchor=[1.01, 0.75], frameon=False)
    plt.tight_layout()
    plt.savefig(f'{indir}/SNP_indels.png', dpi=300)
    plt.savefig(f'{indir}/SNP_indels.pdf')

    return
# END



def get_SVs_plots():
    '''
    Read the VCFdist summary file and generate summary plot for SVs identified using Syri, svim-asm, sniffles2, pbsv, cutesv.
    Statistics considered:
        1) True Positive counts
        2) Precision score
        3) Recall score
        4) F1_QScore
    :return:
    '''
    import pandas as pd
    from matplotlib import pyplot as plt
    from hometools.plot import cleanax
    def getdata(p, t, v):
        df = pd.read_table(p)
        df['tool'] = t
        df['variant'] = v
        return df



    indir = "/dss/dsslegfs01/pn29fi/pn29fi-dss-0016/projects/sv_benchmark/results/variant_calls/"
    tools = ["syri", "svim_asm", "sniffles2", "pbsv", "cutesv"]

    syri = getdata(f'{indir}/syri/vcfdist_unphased_svs/precision-recall-summary.tsv', 'syri', 'SVs')
    syrinofilt = getdata(f'{indir}/syri/vcfdist_nofilt_unphased_svs/precision-recall-summary.tsv', 'syri_nofilt', 'SVs')
    svimasm = getdata(f'{indir}/svim_asm/vcfdist_svs/precision-recall-summary.tsv', 'svim_asm', 'SVs')
    sniffles = getdata(f'{indir}/sniffle2/vcfdist_svs/precision-recall-summary.tsv', 'sniffles2', 'SVs')
    pbsv = getdata(f'{indir}/pbsv//vcfdist_svs/precision-recall-summary.tsv', 'pbsv', 'SVs')
    cutesv = getdata(f'{indir}/cutesv/vcfdist_svs/precision-recall-summary.tsv', 'cutesv', 'SVs')
    cutesvasm = getdata(f'{indir}/cutesv/vcfdist_diploid_v2_svs/precision-recall-summary.tsv', 'cutesv_assembly', 'SVs')


    svdf = pd.concat([syri, syrinofilt, svimasm, sniffles, pbsv, cutesv, cutesvasm])
    svdf = svdf.loc[svdf['VAR_TYPE'] == 'SV']
    svdf = svdf.loc[svdf['THRESHOLD'] == 'NONE']

    """
    syriindels = getdata(f'{indir}/syri/vcfdist_indels/precision-recall-summary.tsv', 'syri', 'indels')
    syrinofiltindels = getdata(f'{indir}/syri/vcfdist_nofilt_indels/precision-recall-summary.tsv', 'syri_nofilt', 'indels')
    dvindels = getdata(f'{indir}/deepvariant/vcfdist_indels/precision-recall-summary.tsv', 'deepvariant', 'indels')
    gatkindels = getdata(f'{indir}/gatk/vcfdist_indels/precision-recall-summary.tsv', 'gatk', 'indels')
    indelsdf = pd.concat([syriindels, syrinofiltindels, dvindels, gatkindels])
    indelsdf = indelsdf.loc[indelsdf['VAR_TYPE'] == 'INDEL']
    indelsdf = indelsdf.loc[indelsdf['THRESHOLD'] == 'NONE']
    """
    bar_colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink']
    bar_labels = ['syri', 'syri_nofilt', 'svim_asm', 'sniffles2', 'pbsv', 'cutesv', 'cutesv_assembly'] 
    bar_ylim = [0.9, 0.9, 10]
    x_stats = 'PREC     RECALL    F1_QSCORE'.split()
    
    
    def get_stats_fig(fig, df, title):
        for j, stat in enumerate(x_stats):
            ax = fig.add_subplot(1, 3, j + 1)
            temp = ax.bar(df['tool'], df[stat].round(4), color=bar_colors, label=bar_labels, alpha=0.75,
                          edgecolor='black',
                          width=0.5)
            ax.set_ylabel(stat)
            ax.set_xticks(ax.get_xticks(), ax.get_xticklabels(), rotation=45, ha='right')
            ax = cleanax(ax)
            plt.suptitle(title)
            plt.tight_layout()
            ax.bar_label(temp)
#            ax.set_ylim([0, bar_ylim[j]])
#            if stat == 'RECALL':
#                ax.set_title('SVs', fontsize=14)

        return fig

    fig = plt.figure(figsize=(14, 4))
    fig = get_stats_fig(fig, svdf, 'SVs')
    #for j, stat in enumerate('PREC    RECALL   F1_QSCORE'.split()):
    
#        ax.bar_label(temp)
#        ax.set_ylim([0, bar_ylim[j]])
#        if stat == 'RECALL':
#            ax.set_title('SNVs' if i == 0 else 'Indels', fontsize=14)
#    for i, df in enumerate([snpdf, indelsdf]):

    # plt.legend(bbox_to_anchor=[1.01, 0.75], frameon=False)
    plt.tight_layout()
    plt.savefig(f'{indir}/SV_v2.png', dpi=300)
    plt.savefig(f'{indir}/SV_v2.pdf')

    return

# END
