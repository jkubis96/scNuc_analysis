
# analysis


# adjustment parameters

from jimg_ncd.nuclei import NucleiFinder

# initiate class
nf = NucleiFinder()

nf.set_adj_image_gamma(gamma = 1.2)

nf.set_adj_image_contrast(contrast = 2)

nf.set_nms(nms = 0.6)

nf.set_prob(prob = 0.3)

nf.set_nuclei_circularity(circ = 0.5)

nf.set_nuclei_size(size = (100,2000))

nf.set_nuclei_min_mean_intensity(intensity = 2000)

nf.set_chromatinization_size(size = (2,1000))

nf.set_chromatinization_ratio(ratio = 0.005)

nf.set_chromatinization_cut_point(cut_point = 1.05)

nf.set_adj_chrom_gamma(gamma = 0.25)

nf.set_adj_chrom_contrast(contrast = 2)

nf.set_adj_chrom_brightness(brightness = 950)


from pickshot import PickShot

import pandas as pd

model = PickShot.download(url = 'https://github.com/jkubis96/PickShot/raw/refs/heads/DAPI/saved_models/DAPI_CNN_(50,50).h5', path_to_save='')


results = model.predict(path_to_images = 'FC-IS_DATA_separated/21q/DAPI', ident_part = 'Ch7', pred_value = 0.7)

results = pd.DataFrame(results)

# series analysis for PicShot images selection

series_results_chromatinization = nf.series_analysis_chromatinization(path_to_images = 'FC-IS_DATA_separated/21q/DAPI', 
                                                  file_extension = 'tif', 
                                                  selected_id = list(results['images_ids'][results['prediciton'] == 'pass']), 
                                                  selection_opt = True, 
                                                  include_img = True, 
                                                  test_series = 0)


from jimg_ncd.nuclei import NucleiDataManagement

# initiate class with NucleiFinder data
ndm = NucleiDataManagement(series_results_chromatinization, '21q')

ndm.save_results_df(path='')

ndm.save_nuc_project(path='')




from jimg_ncd.nuclei import ImagesManagement

indm = ImagesManagement.load_experimental_images(series_results_chromatinization, '21q')

indm.save_raw(path_to_save = '')



# from tqdm import tqdm
# import cv2

# for i in tqdm(series_results_chromatinization.keys()):
#     cv2.imwrite(f'results/21q/{i}.png', img = series_results_chromatinization[i]['img'])


 
del series_results_chromatinization
del ndm


results = model.predict(path_to_images = 'FC-IS_DATA_separated/71q/DAPI', ident_part = 'Ch7', pred_value = 0.7)

results = pd.DataFrame(results)

series_results_chromatinization = nf.series_analysis_chromatinization(path_to_images = 'FC-IS_DATA_separated/71q/DAPI', 
                                                  file_extension = 'tif', 
                                                  selected_id = list(results['images_ids'][results['prediciton'] == 'pass']), 
                                                  selection_opt = True, 
                                                  include_img = True, 
                                                  test_series = 0)


from jimg_ncd.nuclei import NucleiDataManagement

# initiate class with NucleiFinder data
ndm = NucleiDataManagement(series_results_chromatinization, '71q')

ndm.save_results_df(path='')

ndm.save_nuc_project(path='')


from jimg_ncd.nuclei import ImagesManagement

indm = ImagesManagement.load_experimental_images(series_results_chromatinization, '71q')

indm.save_raw(path_to_save = '')



# from tqdm import tqdm
# import cv2

# for i in tqdm(series_results_chromatinization.keys()):
#     cv2.imwrite(f'results/71q/{i}.png', img = series_results_chromatinization[i]['img'])







del series_results_chromatinization
del ndm


results = model.predict(path_to_images = 'FC-IS_DATA_separated/77q/DAPI', ident_part = 'Ch7', pred_value = 0.7)

results = pd.DataFrame(results)


series_results_chromatinization = nf.series_analysis_chromatinization(path_to_images = 'FC-IS_DATA_separated/77q/DAPI', 
                                                  file_extension = 'tif', 
                                                  selected_id = list(results['images_ids'][results['prediciton'] == 'pass']), 
                                                  selection_opt = True, 
                                                  include_img = True, 
                                                  test_series = 0)


from jimg_ncd.nuclei import NucleiDataManagement

# initiate class with NucleiFinder data
ndm = NucleiDataManagement(series_results_chromatinization, '77q')

ndm.save_results_df(path='')

ndm.save_nuc_project(path='')


from jimg_ncd.nuclei import ImagesManagement

indm = ImagesManagement.load_experimental_images(series_results_chromatinization, '77q')

indm.save_raw(path_to_save = '')


# from tqdm import tqdm
# import cv2

# for i in tqdm(series_results_chromatinization.keys()):
#     cv2.imwrite(f'results/77q/{i}.png', img = series_results_chromatinization[i]['img'])





del series_results_chromatinization
del ndm


results = model.predict(path_to_images = 'FC-IS_DATA_separated/109q/DAPI', ident_part = 'Ch7', pred_value = 0.7)

results = pd.DataFrame(results)

series_results_chromatinization = nf.series_analysis_chromatinization(path_to_images = 'FC-IS_DATA_separated/109q/DAPI', 
                                                  file_extension = 'tif', 
                                                  selected_id = list(results['images_ids'][results['prediciton'] == 'pass']), 
                                                  selection_opt = True, 
                                                  include_img = True, 
                                                  test_series = 0)




from jimg_ncd.nuclei import NucleiDataManagement

# initiate class with NucleiFinder data
ndm = NucleiDataManagement(series_results_chromatinization, '109q')

ndm.save_results_df(path='')

ndm.save_nuc_project(path='')



from jimg_ncd.nuclei import ImagesManagement

indm = ImagesManagement.load_experimental_images(series_results_chromatinization, '109q')

indm.save_raw(path_to_save = '')


# from tqdm import tqdm
# import cv2

# for i in tqdm(series_results_chromatinization.keys()):
#     cv2.imwrite(f'results/109q/{i}.png', img = series_results_chromatinization[i]['img'])


###############################################################################


# IS concatenate and projects merging

from jimg_ncd.nuclei import NucleiDataManagement

# merge projects

ndm1 = NucleiDataManagement.load_nuc_dict('21q.nuc')
ndm2 = NucleiDataManagement.load_nuc_dict('71q.nuc')
ndm3 = NucleiDataManagement.load_nuc_dict('77q.nuc')
ndm4 = NucleiDataManagement.load_nuc_dict('109q.nuc')


# load IS data


import pandas as pd

ndm1_is = pd.read_csv('FC-IS_DATA_separated/21q/21q.txt', sep='\t', header=1)
ndm2_is = pd.read_csv('FC-IS_DATA_separated/71q/71q.txt', sep='\t', header=1)
ndm3_is = pd.read_csv('FC-IS_DATA_separated/77q/77q.txt', sep='\t', header=1)
ndm4_is = pd.read_csv('FC-IS_DATA_separated/109q/109q.txt', sep='\t', header=1)


selectes_columns = [ 
   
    'Area_M09',
    'Major Axis_M09',
    'Minor Axis_M09',
    'Aspect Ratio_M09',
    'Diameter_M09',
 ]

ndm1.add_IS_data(ndm1_is, IS_features=selectes_columns)
ndm2.add_IS_data(ndm2_is, IS_features=selectes_columns)
ndm3.add_IS_data(ndm3_is, IS_features=selectes_columns)
ndm4.add_IS_data(ndm4_is, IS_features=selectes_columns)



ndm1.add_experiment([ndm2, ndm3, ndm4])


ndm1.save_mutual_experiments(path='', inc_is=True)


###############################################################################
