import pandas as pd
import seaborn as sns; sns.set(style='white', context='paper')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from statannot import add_stat_annotation
from scipy import stats
from statannot import add_stat_annotation
import itertools


### Graph Swarm Plots of Biomarker Expressions in Each Biopsy ###


ys = ['IL-1β', 'IL-6', 'IL-8', 'IL-10']
location_map = {'SLICE 3': 'Slice 3', 'SLICE 6': 'Slice 6', 'SLICE 9': 'Slice 9', 'Injury': 'Injury', 'Lingula': 'Lingula'}

# Annotate significance levels from an imported table #

def annotate_anova(ax, data, y, anova_path, anova_sheet):
    df = pd.read_excel(anova_path, sheet_name=anova_sheet, index_col=0)
    df = df[y.split(' ')[0]]
    pvalues = []
    box_pairs = []
    for x in df.index:
        p = df[x]
        if p < 0.05:
            pvalues.append(p)
            box_pairs.append(((x, 'RU'), (x, 'LL')))
    add_stat_annotation(ax, data=data, x='EVLP ID', y=y, hue='Location',
                        box_pairs=box_pairs, pvalues=pvalues, perform_stat_test=False,
                        loc='outside', verbose=0)

# Plot horizontal lines on existing graphs #

def horizontal_lines(data, ax, x, y, **kwargs):
    data = data.copy()
    xticks = ax.get_xticklabels()
    xticks = {tick.get_text(): i for i, tick in enumerate(xticks)}
    data['xmin'] = data[x].apply(lambda xval: xticks[str(xval)] - 0.4)
    data['xmax'] = data[x].apply(lambda xval: xticks[str(xval)] + 0.4)
    ax.hlines(y=data[y], xmin=data['xmin'], xmax=data['xmax'], **kwargs)
    return ax

# Make swarm plots with horizontal lines showing mean vluaes and significance annotations #

def make_graph(data_path, data_sheet, anova_path=None, anova_sheet=None, show_fig=False):
    raw_df = pd.read_excel(data_path, sheet_name=data_sheet)
    raw_df = raw_df.drop(columns=['Batch #'])
    df = pd.DataFrame(raw_df)

    fig = plt.figure(figsize=(15, 3))
    gs = GridSpec(1, 4)
    axs = []
    for i in range(len(ys)):
        axs.append(fig.add_subplot(gs[0, i]))

    mypalette = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "orange", "#2ecc71"]
    sns.set_palette(mypalette)

    for ax, y in zip(axs, ys):
        sns.swarmplot(x='Donor ID', y=y, hue='Slice #', data=df, ax=ax, s=3.5)
        ax.axhline(df[y].mean(), color='green', alpha=0.5, linewidth=2, linestyle='dashed')
        df_subj_mean = df.groupby('Donor ID', as_index=False)[y].mean()
        horizontal_lines(data=df_subj_mean, ax=ax, x='Donor ID', y=y, color='black', linewidth=2)
        ax.xaxis.labelpad = 15
        ax.yaxis.labelpad = 5
        ax.set_xlabel(ax.get_xlabel(), fontsize=12.5)
        ax.set_ylabel(ax.get_ylabel(), fontsize=12.5)
        ax.tick_params(axis='both', which='major', labelsize=10.5)
        ax.get_legend().remove()
        if anova_path is not None:
            annotate_anova(ax, df, y, anova_path, anova_sheet)

    handles, labels = axs[0].get_legend_handles_labels()
    labels = [location_map[l] for l in labels]
    fig.legend(handles, labels, loc=(0.05, 0.6))

    fig.tight_layout()

    if show_fig:
        fig.show()
    else:
        fig.savefig(f'Swarm Plots/SwarmPlot in a row.png', dpi=200)

make_graph(r'C:\Users\chaob\Documents\Biopsy Heterogeneity Data Sheet.xlsx', 'Sheet1')


### Make Violin and Box Plots Using Normalized Data ###


def cat_plot(type, x, y, file_name, hue=None):
    df_cat = pd.read_excel(r'C:\Users\chaob\Documents\Biopsy Heterogeneity Data Sheet.xlsx',
                              sheet_name='Violin Plot')

    fig, ax = plt.subplots(figsize=(5, 6))

    if type == 'violin':
        sns.violinplot(x=x, y=y, hue=hue, palette='pastel', data=df_cat, ax=ax)
        fig.tight_layout()
        fig.savefig(file_name, dpi=200)

    if type == 'box':
        sns.boxplot(x=x, y=y, hue=hue, saturation=0.5, showfliers=False, palette='pastel', data=df_cat, ax=ax)
        ax.xaxis.labelpad = 15
        ax.yaxis.labelpad = 5
        ax.set_xlabel(ax.get_xlabel(), fontsize=13)
        ax.set_ylabel(ax.get_ylabel(), fontsize=13)
        ax.tick_params(axis='both', which='major', labelsize=10.5)
        xlabels = [l.get_text() for l in ax.get_xticklabels()]
        if hue is None:
            box_pairs = list(itertools.combinations(xlabels, 2))
            sns.stripplot(x=x, y=y, hue=hue, s=3, data=df_cat, alpha=0.6, palette='tab10', ax=ax)
        else:
            huelabels = df_cat[hue].unique().tolist()
            hue_pairs = list(itertools.combinations(huelabels, 2))
            box_pairs = []
            for xlabel in xlabels:
                for hue1, hue2 in hue_pairs:
                    pair = ((xlabel, hue1), (xlabel, hue2))
                    box_pairs.append(pair)
        add_stat_annotation(ax, data=df_cat, x=x, y=y, hue=hue, box_pairs=box_pairs,
                                perform_stat_test=True, test='t-test_welch',
                                loc='inside', verbose=0, no_ns=True, fontsize='large')
        fig.tight_layout()
        fig.savefig(file_name, dpi=200)

cat_plot('violin', 'Biopsy', 'Cytokine Expression', f'Swarm Plots/Violin plot by injury.png')
cat_plot('box', 'Biopsy', 'Cytokine Expression', f'Swarm Plots/Overall box plot by injury.png')
cat_plot('violin', 'Reason for Rejection', 'Cytokine Expression', f'Swarm Plots/Overall violin plot by injury.png', hue='Biopsy')
cat_plot('box', 'Reason for Rejection', 'Cytokine Expression', f'Swarm Plots/Box plot by injury.png', hue='Biopsy')

# By Cytokine #

# fig_violin_ByCytokine = plt.figure(figsize=(12, 12))
# sns.violinplot(x='Cytokine', y='Cytokine Expression', hue='Biopsy', data=df_violin)
# fig_violin_ByCytokine.tight_layout()
# fig_violin_ByCytokine.savefig(f'Swarm Plots/Violin plot by cytokine.png', dpi=200)
#
# fig_box_ByCytokine = plt.figure(figsize=(12, 12))
# sns.boxplot(x='Cytokine', y='Cytokine Expression', hue='Biopsy', data=df_violin)
# fig_box_ByCytokine.tight_layout()
# fig_box_ByCytokine.savefig(f'Swarm Plots/Box plot by cytokine.png', dpi=200)


### Intra-class Correlations ###


import pingouin as pg

file_paths = [r'C:\Users\chaob\Documents\Data by Donor (3 biopsies).xlsx',
              r'C:\Users\chaob\Documents\Data by Donor (4 biopsies).xlsx',
              r'C:\Users\chaob\Documents\Data by Donor (5 biopsies).xlsx',]
sheet_names = ['#155', '#211', '#250', '#288', '#308', '#318', '#289', '#302']

for p in file_paths:
    for s in sheet_names:
        df = pd.read_excel(p, s)
        icc = pg.intraclass_corr(data=df, targets='Target Cytokine', raters='Slice #',
                             ratings='Expression').round(3)
        print(p, s, '\n', icc)


### Plot Distributions of Coefficients of Variance (CV) in Different Biopsy Groups Compared to the Injury Status ###


df_distribution = pd.read_excel(r'C:\Users\chaob\Documents\Biopsy Heterogeneity Data Sheet.xlsx',
                               sheet_name='Distribution')
df_cat = pd.read_excel(r'C:\Users\chaob\Documents\Biopsy Heterogeneity Data Sheet.xlsx',
                              sheet_name='Violin Plot')

fig_hist_369uninjured = plt.figure(figsize=(12, 12))
plt.hist(df_distribution['Log of intra-Biopsy %CV'])
plt.hist(df_distribution['Log of 3, 6, 9 %CV (Uninjured)'])
plt.legend(['Log of intra-Biopsy %CV', 'Log of 3, 6, 9 %CV (Uninjured)'])
fig_hist_369uninjured.savefig(f'Swarm Plots/Histogram of 369 Uninjured.png', dpi=200)

fig_hist_369injured = plt.figure(figsize=(12, 12))
plt.hist(df_distribution['Log of intra-Biopsy %CV'])
plt.hist(df_distribution['Log of 3, 6, 9 %CV (Injured)'])
plt.legend(['Log of intra-Biopsy %CV', 'Log of 3, 6, 9 %CV (Injured)'])
fig_hist_369injured.savefig(f'Swarm Plots/Histogram of 369 Injured.png', dpi=200)

fig_hist_369LingulaUninjured = plt.figure(figsize=(12, 12))
plt.hist(df_distribution['Log of intra-Biopsy %CV'])
plt.hist(df_distribution['Log of 3, 6, 9 + Lingula %CV (Uninjured)'])
plt.legend(['Log of intra-Biopsy %CV', 'Log of 3, 6, 9 + Lingula %CV (Uninjured)'])
fig_hist_369LingulaUninjured.savefig(f'Swarm Plots/Histogram of 369+Lingula Uninjured.png', dpi=200)

fig_hist_369LingulaInjured = plt.figure(figsize=(12, 12))
plt.hist(df_distribution['Log of intra-Biopsy %CV'])
plt.hist(df_distribution['Log of 3, 6, 9 + Lingula %CV (Injured)'])
plt.legend(['Log of intra-Biopsy %CV', 'Log of 3, 6, 9 + Lingula %CV (Injured)'])
fig_hist_369LingulaInjured.savefig(f'Swarm Plots/Histogram of 369+Lingula Injured.png', dpi=200)


### Various Statistical Tests ###


# Calculate p-values using t-tests and mann-whitney tests #

p_IntraBiopsy_369Uninjured = stats.ttest_ind(df_distribution['Intra-Biopsy %CV'], df_distribution['3, 6, 9 %CV (Uninjured)'], nan_policy='omit', equal_var=True)
p_IntraBiopsy_369Injured = stats.ttest_ind(df_distribution['Intra-Biopsy %CV'], df_distribution['3, 6, 9 %CV (Injured)'], nan_policy='omit', equal_var=False)
p_IntraBiopsy_369LingulaUninjured = stats.ttest_ind(df_distribution['Intra-Biopsy %CV'], df_distribution['3, 6, 9 + Lingula %CV (Uninjured)'], nan_policy='omit', equal_var=False)
p_IntraBiopsy_369LingulaInjured = stats.ttest_ind(df_distribution['Intra-Biopsy %CV'], df_distribution['3, 6, 9 + Lingula %CV (Injured)'], nan_policy='omit', equal_var=False)
p_369Uninjured_369Injured = stats.ttest_ind(df_distribution['3, 6, 9 %CV (Uninjured)'], df_distribution['3, 6, 9 %CV (Injured)'], nan_policy='omit', equal_var=False)
p_369LingulaUninjured_369LingulaInjured = stats.ttest_ind(df_distribution['3, 6, 9 + Lingula %CV (Uninjured)'], df_distribution['3, 6, 9 + Lingula %CV (Injured)'], nan_policy='omit', equal_var=False)

p_IntraBiopsy_369Uninjured = stats.mannwhitneyu(df_distribution['Intra-Biopsy %CV'].dropna(), df_distribution['3, 6, 9 %CV (Uninjured)'].dropna())
p_IntraBiopsy_369Injured = stats.mannwhitneyu(df_distribution['Intra-Biopsy %CV'].dropna(), df_distribution['3, 6, 9 %CV (Injured)'].dropna())
p_IntraBiopsy_369LingulaUninjured = stats.mannwhitneyu(df_distribution['Intra-Biopsy %CV'].dropna(), df_distribution['3, 6, 9 + Lingula %CV (Uninjured)'].dropna())
p_IntraBiopsy_369LingulaInjured = stats.mannwhitneyu(df_distribution['Intra-Biopsy %CV'].dropna(), df_distribution['3, 6, 9 + Lingula %CV (Injured)'].dropna())
p_369Uninjured_369Injured = stats.mannwhitneyu(df_distribution['3, 6, 9 %CV (Uninjured)'].dropna(), df_distribution['3, 6, 9 %CV (Injured)'].dropna())
p_369LingulaUninjured_369LingulaInjured = stats.mannwhitneyu(df_distribution['3, 6, 9 + Lingula %CV (Uninjured)'].dropna(), df_distribution['3, 6, 9 + Lingula %CV (Injured)'].dropna())

print(p_IntraBiopsy_369Uninjured, p_IntraBiopsy_369Injured, p_IntraBiopsy_369LingulaUninjured, p_IntraBiopsy_369LingulaInjured, p_369Uninjured_369Injured, p_369LingulaUninjured_369LingulaInjured)

# Test for skewness #

skew_IntraBiopsy = stats.skew(df_distribution['Intra-Biopsy %CV'], nan_policy='omit')
skew_369Uninjured = stats.skew(df_distribution['3, 6, 9 %CV (Uninjured)'], nan_policy='omit')
skew_369Injured = stats.skew(df_distribution['3, 6, 9 %CV (Injured)'], nan_policy='omit')
skew_369LingulaUninjured = stats.skew(df_distribution['3, 6, 9 + Lingula %CV (Uninjured)'], nan_policy='omit')
skew_369LingulaInjured = stats.skew(df_distribution['3, 6, 9 + Lingula %CV (Injured)'], nan_policy='omit')

skew_Normalized369 = stats.skew(df_cat['Cytokine Expression'][df_cat['Biopsy'] == '3 Slices'], nan_policy='omit')
skew_NormalizedInjury = stats.skew(df_cat['Cytokine Expression'][df_cat['Biopsy'] == 'Injury'], nan_policy='omit')
skew_NormalizedLingula = stats.skew(df_cat['Cytokine Expression'][df_cat['Biopsy'] == 'Lingula'], nan_policy='omit')

print(skew_IntraBiopsy, skew_369Uninjured, skew_369Injured, skew_369LingulaUninjured, skew_369LingulaInjured)
print(skew_Normalized369, skew_NormalizedInjury, skew_NormalizedLingula)

# P-value of Normalized Biopsy Data #

fig_hist_normalized = plt.figure(figsize=(12, 12))
plt.hist(df_cat['Cytokine Expression'][df_cat['Biopsy'] == '3 Slices'])
plt.hist(df_cat['Cytokine Expression'][df_cat['Biopsy'] == 'Injury'])
plt.hist(df_cat['Cytokine Expression'][df_cat['Biopsy'] == 'Lingula'])
plt.legend(['Normalized three slices', 'Normalized injury', 'Normalized lingula'])
fig_hist_normalized.savefig(f'Swarm Plots/Histogram of Normalized Biopsy Data.png', dpi=200)

p_Normalized_369Injury = stats.ttest_ind(df_cat['Cytokine Expression'][df_cat['Biopsy'] == '3 Slices'],
                                            df_cat['Cytokine Expression'][df_cat['Biopsy'] == 'Injury'], nan_policy='omit', equal_var=False)
p_Normalized_369Lingula = stats.ttest_ind(df_cat['Cytokine Expression'][df_cat['Biopsy'] == '3 Slices'],
                                            df_cat['Cytokine Expression'][df_cat['Biopsy'] == 'Lingula'], nan_policy='omit', equal_var=False)
p_Normalized_InjuryLingula = stats.ttest_ind(df_cat['Cytokine Expression'][df_cat['Biopsy'] == 'Injury'],
                                            df_cat['Cytokine Expression'][df_cat['Biopsy'] == 'Lingula'], nan_policy='omit', equal_var=False)

print(p_Normalized_369Injury, p_Normalized_369Lingula, p_Normalized_InjuryLingula)

test_1 = stats.levene(df_distribution['Intra-Biopsy %CV'].dropna(), df_distribution['3, 6, 9 + Lingula %CV (Injured)'].dropna())
print(test_1)

df_cat = pd.read_excel(r'C:\Users\chaob\Documents\Biopsy Heterogeneity Data Sheet.xlsx',
                               sheet_name='Violin Plot')

test_2 = stats.levene(df_cat['Cytokine Expression'][df_cat['Biopsy'] == '3 Slices'].dropna(), df_cat['Cytokine Expression'][df_cat['Biopsy'] == 'Injury'].dropna())
test_3 = stats.levene(df_cat['Cytokine Expression'][df_cat['Biopsy'] == '3 Slices'].dropna(), df_cat['Cytokine Expression'][df_cat['Biopsy'] == 'Lingula'].dropna())
test_4 = stats.levene(df_cat['Cytokine Expression'][df_cat['Biopsy'] == 'Lingula'].dropna(), df_cat['Cytokine Expression'][df_cat['Biopsy'] == 'Injury'].dropna())
print(test_2, test_3, test_4)

p_369_Injury = stats.ttest_ind(df_cat['Cytokine Expression'][df_cat['Biopsy'] == '3 Slices'], df_cat['Cytokine Expression'][df_cat['Biopsy'] == 'Injury'], nan_policy='omit', equal_var=False)
p_369_Lingula = stats.ttest_ind(df_cat['Cytokine Expression'][df_cat['Biopsy'] == '3 Slices'], df_cat['Cytokine Expression'][df_cat['Biopsy'] == 'Lingula'], nan_policy='omit', equal_var=False)
p_Lingula_Injury = stats.ttest_ind(df_cat['Cytokine Expression'][df_cat['Biopsy'] == 'Lingula'], df_cat['Cytokine Expression'][df_cat['Biopsy'] == 'Injury'], nan_policy='omit', equal_var=True)
print(p_369_Injury, p_369_Lingula, p_Lingula_Injury)