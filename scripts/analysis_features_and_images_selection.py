# analysis


from jimg_ncd.nuclei import GroupAnalysis
import pandas as pd

data = pd.read_csv('21q_71q_77q_109q_IS.csv', sep=',', header=0)

# initiate class
ga = GroupAnalysis.load_data(data, ids_col = 'id_name', set_col = 'set')

# check available groups for selection of differential features
ga.groups

# run DFA analysis on example sets

ga.DFA(meta_group_by = 'sets',
   sets = {'disease':['109q', '71q', '77q'],
           'ctrl':['21q']}, 
   n_proc = 5)
 
group_diff_features = ga.get_DFA()



ga.heatmap_DFA(top_n = 10)

fig = ga.get_DFA_plot()

fig.savefig('DFA_groups.svg', dpi=300, bbox_inches='tight')


ga.heatmap_DFA(top_n = 10, scale = True)

fig = ga.get_DFA_plot()

fig.savefig('scaled_DFA_groups.svg', dpi=300, bbox_inches='tight')

group_diff_features.to_excel('DFA_groups.xlsx')

###############################################################################



ga.data_scale()


# run PCA dimensionality reduction
ga.PCA()


# get PCA data, if required
pca_data = ga.get_PCA()


# run PC variance analysis
ga.var_plot()


# get var_data, if required


# get knee_plot, if required
knee_plot = ga.get_knee_plot(show = True)

knee_plot.savefig('knee.svg', dpi=300)

ga.UMAP(PC_num = 10,
     factorize_with_metadata = False, 
     harmonize_sets = True,
     n_neighbors = 5,
     min_dist = 0.02,
     n_components = 2)


# get UMAP_data, if required
UMAP_data = ga.get_UMAP_data()

pd.DataFrame(UMAP_data).to_excel('UMAP.xlsx')


# get UMAP_plots, if required
UMAP_plots = ga.get_UMAP_plots(show = True)

UMAP_plots.keys()

UMAP_plots['PrimaryUMAP'].savefig('UMAP.svg', dpi=300)
                                  
                                  
                                  
ga.db_scan(eps = 0.46,
        min_samples = 50)


# run UMAP_on_clusters
ga.UMAP_on_clusters(min_entities = 50)


# get UMAP_plots, if required
UMAP_plots = ga.get_UMAP_plots(show = True)

UMAP_plots.keys()

UMAP_plots['ClusterUMAP'].savefig('UMAP_clusters.png', dpi=300)
UMAP_plots['ClusterXSetsUMAP'].savefig('ClusterXSetsUMAP.png', dpi=300)
UMAP_plots['ClusterUMAP'].savefig('UMAP_clusters.svg', dpi=300)
UMAP_plots['ClusterXSetsUMAP'].savefig('ClusterXSetsUMAP.svg', dpi=300)




UMAP_plots = ga.get_UMAP_plots(show = True, plot_type = 'html')

UMAP_plots.keys()

import plotly.io as pio


pio.write_html(UMAP_plots['ClusterUMAP'], 'UMAP_clusters.html', auto_open=False)
pio.write_html(UMAP_plots['ClusterXSetsUMAP'], 'ClusterXSetsUMAP.html', auto_open=False)

# run DFA analysis on finl clusters
ga.DFA(meta_group_by = 'clusters',
    sets = {}, 
    n_proc = 5)

dfa_clusters = ga.get_DFA()

dfa_clusters.to_excel('markers_clusters.xlsx')


ga.heatmap_DFA(top_n = 5, scale = False,  figsize = (10,6))

fig = ga.get_DFA_plot()

fig.savefig('markers_heatmap.svg', dpi=300)



ga.heatmap_DFA(top_n = 5, scale = True, p_value = 0.05, figsize = (10,6))

fig = ga.get_DFA_plot()

fig.savefig('scaled_markers_heatmap.svg', dpi=300)

ga.print_avaiable_features()

          
ga.proportion_analysis(grouping_col = 'sets', 
                       val_col = 'nuclei_per_img', 
                       grouping_dict = None, 
                       omit = None)
     
pl = ga.get_proportion_plot(show = True)

fig.savefig('nuclei_stats.svg', dpi=300)

stats = ga.proportion_stats



ga.proportion_analysis(grouping_col = 'sets', 
                       val_col = 'spot_n', 
                       grouping_dict = None, 
                       omit = None)
     
pl = ga.get_proportion_plot(show = True)

fig.savefig('chromatinization_stats.svg', dpi=300)

ga.proportion_stats



# save data

to_save = ga.full_info()

to_save.to_excel('full_data.xlsx')

###############################################################################

import os
import shutil

for c in set(dfa_clusters['valid_group']):
    tmp = dfa_clusters[dfa_clusters['valid_group'] == c]
    tmp = tmp[tmp['log(FC)'] > 0]
    tmp = tmp[tmp['adj_pval'] < 0.05].sort_values('esm',  ascending = False)
    
    tmp2 = to_save.loc[:,list(tmp['feature']) + ['clusters', 'set', 'id_y']]
    tmp2 = tmp2[tmp2['clusters'] == c]
    
    path = f'results/clusters/{c}'
    os.makedirs(path, exist_ok=True)
    
    for s in set(tmp2['set']):
        print(s)
        tmp3 = tmp2[tmp2['set'] == s]
        
        feature_cols = [col for col in tmp3.columns if col not in ['clusters', 'set', 'id_y']]

        medians = tmp3[feature_cols].median()
        
        tmp3['median_distance'] = (tmp3[feature_cols] - medians).abs().sum(axis=1)
        
        top3_idx = tmp3['median_distance'].nsmallest(5).index
        top3_ids = tmp3.loc[top3_idx, 'id_y'].tolist()
                
        folder_path = f'results/{s}_stitched_bright'
        
        files = [f for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f))]
        
        top3_files = [f for f in files if any(str(id_) in f for id_ in top3_ids)]
        
        for f in top3_files:
            src_path = os.path.join(folder_path, f)
            
            original_name = os.path.splitext(f)[0]  
            new_name = f'{original_name}_cluster{c}_{s}.png'
            dst_path = os.path.join(path, new_name)
            shutil.copy2(src_path, dst_path)



import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

vis_data = to_save[to_save['clusters'] != '11']

df = (
    vis_data
        .groupby(['clusters', 'set'])
        .size()
        .reset_index(name='count')
)

set_content = vis_data['set'].value_counts()

df['normalized'] = df['count'] / df['set'].map(set_content)

df['pct'] = (
    df['normalized'] /
    df.groupby('clusters')['normalized'].transform('sum')
)


cmap_name = "Accent"
width = 2
height = 6
font_size = 20
legend_split_col = 1
legend_bbox = (1.3, 0.5)

clusters = sorted(df["clusters"].unique(), key=lambda x: int(x))

sets_unique = ["21q", "71q", "77q", "109q"]

fig, axes = plt.subplots(1, len(clusters), figsize=(width * len(clusters), height))

if len(clusters) == 1:
    axes = [axes]

cmap = plt.get_cmap(cmap_name)
colors = [cmap(i / len(sets_unique)) for i in range(len(sets_unique))]
cmap_dict = dict(zip(sets_unique, colors))


for idx, cluster in enumerate(clusters):
    ax = axes[idx]

    tmp_df = df[df["clusters"] == cluster].reset_index(drop=True)

    values = tmp_df["pct"].values * 100  # wysokość segmentów
    values = np.round(values, 2)

    names = tmp_df["set"].values
    perc = tmp_df["pct"].values * 100
    nums = tmp_df["count"].values

    bottom = 0
    centers = []

    for name, val in zip(names, values):
        ax.bar(
            0,
            val,
            bottom=bottom,
            color=cmap_dict[name],
            width=0.8
        )
        centers.append(bottom + val / 2)
        bottom += val

    for y_center, pct in zip(centers, perc):
        ax.annotate(
            f"{pct:.1f}%",
            xy=(0, y_center),
            xytext=(-0.9, y_center),
            ha="right",
            va="center",
            fontsize=font_size,
            arrowprops=dict(
                arrowstyle="->",
                lw=1,
                color="black"
            ),
        )

    ax.set_ylim(0, 100)
    ax.set_xlim(-1, 1)
    ax.set_xlabel(f"Cluster {cluster}", fontsize=font_size)
    ax.set_xticks([])
    ax.set_yticks([])

    for spine in ax.spines.values():
        spine.set_visible(False)


legend_handles = [
    Patch(facecolor=cmap_dict[label], edgecolor="black", label=label)
    for label in sets_unique
]

fig.legend(
    handles=legend_handles,
    loc="center right",
    bbox_to_anchor=legend_bbox,
    ncol=legend_split_col
)

plt.tight_layout()
plt.show()
       

fig.savefig('proportion.svg', dpi=300, bbox_inches='tight')



###############################################################################
        


from jimg_ncd.nuclei import GroupAnalysis
import pandas as pd

data = pd.read_csv('21q_71q_77q_109q_IS.csv', sep=',', header=0)
data['full_id'] = data['id_name'].astype(str) + ' # ' + data['set']

data = data[data['full_id'].isin(list(set(vis_data['full_id'])))]


# initiate class
ga = GroupAnalysis.load_data(data, ids_col = 'id_name', set_col = 'set')

# check available groups for selection of differential features
ga.groups

# run DFA analysis on example sets

ga.DFA(meta_group_by = 'sets',
   sets = {'disease':['109q', '71q', '77q'],
           'ctrl':['21q']}, 
   n_proc = 5)
 
group_diff_features = ga.get_DFA()



ga.heatmap_DFA(top_n = 10)

fig = ga.get_DFA_plot()

fig.savefig('DFA_groups_reduced.svg', dpi=300, bbox_inches='tight')


ga.heatmap_DFA(top_n = 10, scale = True)

fig = ga.get_DFA_plot()

fig.savefig('scaled_DFA_groups_reduced.svg', dpi=300, bbox_inches='tight')

group_diff_features.to_excel('DFA_groups_reduced.xlsx')




###############################################################################
        
import pandas as pd

to_save = pd.read_excel('results/full_data.xlsx')

to_save = to_save[to_save['full_id'].isin(list(set(vis_data['full_id'])))]

to_save.to_excel('full_data_reduced.xlsx')

import os
import shutil

features = group_diff_features['feature'][group_diff_features['adj_pval'] < 0.05]

for c in features:
    tmp2 = to_save.loc[:, [c] + ['clusters', 'set', 'id_y']]
    
    for s in set(tmp2['set']):
        print(s)
        tmp3 = tmp2[tmp2['set'] == s]
        
        feature_cols = [col for col in tmp3.columns if col not in ['clusters', 'set', 'id_y', 'set2']]

        medians = tmp3[feature_cols].median()
        q25 = tmp3[feature_cols].quantile(0.10)
        q75 = tmp3[feature_cols].quantile(0.90)
        
        quantiles = {'median': medians, 'q25': q25, 'q75': q75}
        
        for q_name, q_values in quantiles.items():
            path = f'results/stats/{q_name}/{c}'
            os.makedirs(path, exist_ok=True)
            
            tmp3[f'{q_name}_distance'] = (tmp3[feature_cols] - q_values).abs().sum(axis=1)
            
            top5_idx = tmp3[f'{q_name}_distance'].nsmallest(5).index
            top5_ids = tmp3.loc[top5_idx, 'id_y'].tolist()
            
            folder_path = f'results/{s}_stitched_bright'
            files = [f for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f))]
            
            top5_files = [f for f in files if any(str(id_) in f for id_ in top5_ids)]
            
            for f in top5_files:
                src_path = os.path.join(folder_path, f)
                original_name = os.path.splitext(f)[0]
                new_name = f'{original_name}_stat_{c}_{s}_{q_name}.png'
                dst_path = os.path.join(path, new_name)
                shutil.copy2(src_path, dst_path)
            
###############################################################################



from jimg_ncd.nuclei import ImagesManagement

img_mg = ImagesManagement.load_from_dict('21q.inuc.npz', '21q')

img_mg.adjust_images(acronyme='bright', color='gray', path_to_images='FC-IS_DATA_separated/21q/21q', fille_name_part='Ch9.ome', eq = False, clahe = False, gamma = 0.7, contrast = 1.2)

img_mg.get_included_acronyms()

img_mg.image_stitching(acronyms = ['bright'], to_results_image = True)  

img_mg.save_prepared_images(acronyme = 'stitched_bright')




img_mg = ImagesManagement.load_from_dict('71q.inuc.npz', '71q')

img_mg.adjust_images(acronyme='bright', color='gray', path_to_images='FC-IS_DATA_separated/71q/71q', fille_name_part='Ch9.ome', eq = False, clahe = False, gamma = 0.7, contrast = 1.2)

img_mg.get_included_acronyms()

img_mg.image_stitching(acronyms = ['bright'], to_results_image = True)  

img_mg.save_prepared_images(acronyme = 'stitched_bright')




img_mg = ImagesManagement.load_from_dict('77q.inuc.npz', '77q')

img_mg.adjust_images(acronyme='bright', color='gray', path_to_images='FC-IS_DATA_separated/77q/77q', fille_name_part='Ch9.ome', eq = False, clahe = False, gamma = 0.7, contrast = 1.2)

img_mg.get_included_acronyms()

img_mg.image_stitching(acronyms = ['bright'], to_results_image = True)  

img_mg.save_prepared_images(acronyme = 'stitched_bright')



img_mg = ImagesManagement.load_from_dict('109q.inuc.npz', '109q')

img_mg.adjust_images(acronyme='bright', color='gray', path_to_images='FC-IS_DATA_separated/109q/109q', fille_name_part='Ch9.ome', eq = False, clahe = False, gamma = 0.7, contrast = 1.2)

img_mg.get_included_acronyms()

img_mg.image_stitching(acronyms = ['bright'], to_results_image = True)  

img_mg.save_prepared_images(acronyme = 'stitched_bright')





img_mg = ImagesManagement.load_from_dict('results/109q.inuc.npz', '109q')

img_mg.adjust_images(acronyme='bright', color='gray', path_to_images='FC-IS_DATA_separated/109q/109q', fille_name_part='Ch9.ome', eq = False, clahe = False, gamma = 0.7, contrast = 1.2)


img_mg.adjust_images(acronyme='DAPI', color='blue', path_to_images='FC-IS_DATA_separated/109q/109q', fille_name_part='Ch7.ome', eq = True, clahe = True, gamma =2, contrast = 1.2)

img_mg.image_merging(acronyms = ['bright', 'DAPI'], ratio_list = [0.9,0.6])  

img_mg.save_prepared_images(acronyme = 'merged_bright_DAPI')

img_mg.get_included_acronymes()

img_mg.image_stitching(acronymes = ['merged_bright_DAPI'], to_results_image = True)  

img_mg.save_prepared_images(acronyme = 'merged_bright_DAPI')




###############################################################################



import pandas as pd

to_save = pd.read_excel('results/full_data.xlsx')
to_save['nuc_cell_ratio'] = to_save['nuclei_area']/to_save['Area_M09']


fet_dict = {'5':['Area_M09', 'Diameter_M09', 'Major Axis_M09'],
            '7':['avg_spot_perimeter', 'sum_spot_perimeter'],
            '14':['nuc_cell_ratio', 'Area_M09', 'Diameter_M09']}

import os
import shutil

for c in fet_dict.keys():
    print(c)
   
    
    tmp2 = to_save.loc[:,fet_dict[c] + ['clusters', 'set', 'id_y']]
    tmp2 = tmp2[tmp2['clusters'] == int(c)]
    
    path = f'results/stats_clusters/{c}'
    os.makedirs(path, exist_ok=True)
    
    for s in set(tmp2['set']):
        print(s)
        tmp3 = tmp2[tmp2['set'] == s]
        
        feature_cols = [col for col in tmp3.columns if col not in ['clusters', 'set', 'id_y']]

        medians = tmp3[feature_cols].median()
        
        tmp3['median_distance'] = (tmp3[feature_cols] - medians).abs().sum(axis=1)
        
        top3_idx = tmp3['median_distance'].nsmallest(10).index
        top3_ids = tmp3.loc[top3_idx, 'id_y'].tolist()
                
        folder_path = f'results/{s}_stitched_bright'
        
        files = [f for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f))]
        
        top3_files = [f for f in files if any(str(id_) in f for id_ in top3_ids)]
        
        for f in top3_files:
            src_path = os.path.join(folder_path, f)
            
            # tworzymy nową nazwę: oryginalna_nazwa + c + s
            original_name = os.path.splitext(f)[0]  # usuwamy rozszerzenie np. ".png"
            new_name = f'{original_name}_cluster{c}_{s}.png'
            dst_path = os.path.join(path, new_name)
            shutil.copy2(src_path, dst_path)





################################################################################


# prepare images for www




import pandas as pd

to_save = pd.read_excel('results/full_data.xlsx')


import os
import re
import shutil
import json


saved_files = {}

for c in set(to_save['clusters']):
    print(c)

    tmp = to_save[to_save['clusters'] == int(c)]

    path = f'docs/images/{c}'
    os.makedirs(path, exist_ok=True)

    saved_files[str(c)] = {}

    for s in set(tmp['set']):
        print(s)

        tmp2 = tmp[tmp['set'] == s]

        folder_path = f'results/{s}_stitched_bright'
        files = [f for f in os.listdir(folder_path)
                 if os.path.isfile(os.path.join(folder_path, f))]

        to_select = list(tmp2['id_y'])
        files_selected = [f for f in files if any(str(id_) in f for id_ in to_select)]

        saved_files[str(c)][str(s)] = []

        for f in files_selected:
            src_path = os.path.join(folder_path, f)

            original_name = os.path.splitext(f)[0]
            original_name = re.sub('_.*', '', original_name)

            new_name = f'{original_name}_set_{s}.png'
            dst_path = os.path.join(path, new_name)

            shutil.copy2(src_path, dst_path)

            # zapisujemy tylko nazwę pliku
            saved_files[str(c)][str(s)].append(new_name)
            
with open("saved_images.json", "w", encoding="utf-8") as f:
    json.dump(saved_files, f, indent=4, ensure_ascii=False)


