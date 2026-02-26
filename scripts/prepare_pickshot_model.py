# adjustment parameters
from jimg_ncd.nuclei import NucleiFinder

# initiate class
nf = NucleiFinder()

image = nf.load_image("FC-IS_DATA_separated/21q/DAPI/3076_Ch7.ome.tif")

nf.input_image(image)

nf.set_adj_image_gamma(gamma = 1.2)

nf.set_adj_image_contrast(contrast = 2)

# Check the basic parameters
nf.current_parameters_nuclei


# Test nms & prob parmeters for nuclei segmentation
nf.nuclei_finder_test()

nf.browser_test()

nf.set_nms(nms = 0.6)

nf.set_prob(prob = 0.3)


nf.find_nuclei()



nf.set_nuclei_circularity(circ = 0.5)


nf.set_nuclei_size(size = (100,800))

nf.set_nuclei_min_mean_intensity(intensity = 2000)


# Check if parameters has changed
nf.current_parameters_nuclei


# Second execution with adjusted parameters of second step of analysis (selection)
nf.select_nuclei()




# Chromatinization parameters
nf.current_parameters_chromatinization

nf.set_chromatinization_size(size = (2,1000))

nf.set_chromatinization_ratio(ratio = 0.005)

nf.set_chromatinization_cut_point(cut_point = 1.05)

nf.current_parameters_chromatinization


# Chromatinization image parameters
nf.current_parameters_img_adj_chro

nf.set_adj_chrom_gamma(gamma = 0.1)

nf.set_adj_chrom_contrast(contrast = 2)

nf.set_adj_chrom_brightness(brightness = 950)

nf.current_parameters_img_adj_chro


# Second execution of the third step (chromatinization analysis)
nf.nuclei_chromatinization()



# series analysis for PicShot images selection

series_results_chromatinization = nf.series_analysis_chromatinization(path_to_images = 'FC-IS_DATA_separated/21q/DAPI', 
                                                  file_extension = 'tif', 
                                                  selected_id = [], 
                                                  selection_opt = True, 
                                                  include_img = True, 
                                                  test_series = 500)

from tqdm import tqdm
import cv2

for i in tqdm(series_results_chromatinization.keys()):
    cv2.imwrite(f'picskshot_data/21q/{i}.png', img = series_results_chromatinization[i]['img'])


 

series_results_chromatinization = nf.series_analysis_chromatinization(path_to_images = 'FC-IS_DATA_separated/71q/DAPI', 
                                                  file_extension = 'tif', 
                                                  selected_id = [], 
                                                  selection_opt = True, 
                                                  include_img = True, 
                                                  test_series = 500)

from tqdm import tqdm
import cv2

for i in tqdm(series_results_chromatinization.keys()):
    cv2.imwrite(f'picskshot_data/71q/{i}.png', img = series_results_chromatinization[i]['img'])


series_results_chromatinization = nf.series_analysis_chromatinization(path_to_images = 'FC-IS_DATA_separated/77q/DAPI', 
                                                  file_extension = 'tif', 
                                                  selected_id = [], 
                                                  selection_opt = True, 
                                                  include_img = True, 
                                                  test_series = 500)

from tqdm import tqdm
import cv2

for i in tqdm(series_results_chromatinization.keys()):
    cv2.imwrite(f'picskshot_data/77q/{i}.png', img = series_results_chromatinization[i]['img'])




series_results_chromatinization = nf.series_analysis_chromatinization(path_to_images = 'FC-IS_DATA_separated/109q/DAPI', 
                                                  file_extension = 'tif', 
                                                  selected_id = [], 
                                                  selection_opt = True, 
                                                  include_img = True, 
                                                  test_series = 500)

from tqdm import tqdm
import cv2

for i in tqdm(series_results_chromatinization.keys()):
    cv2.imwrite(f'picskshot_data/109q/{i}.png', img = series_results_chromatinization[i]['img'])




# model learning


# collect drop and save images paths

samples = ['21q', '71q', '77q', '109q']

import os
import re

full_drop = []
full_good = []

for o in samples:
    print(o)    
    
    
    file_list_dropout = os.listdir(f'picskshot_data/{o}/bad')
    file_list_dropout = [os.path.join(f'FC-IS_DATA_separated/{o}/DAPI', re.sub(r"\.png$", "", file)) for file in file_list_dropout]
    
    file_list_good = os.listdir(f'picskshot_data/{o}/good')
    file_list_good = [os.path.join(f'FC-IS_DATA_separated/{o}/DAPI', re.sub(r"\.png$", "", file)) for file in file_list_good]
    
    
    full_drop = full_drop + file_list_dropout 
    full_good = full_good + file_list_good 
    



from pickshot import TrainingModel

# create instance of the TrainingModel class
model = TrainingModel()




model.images_paths(images_to_drop = full_drop, images_to_save = full_good)

model.train()

model.save_model(name = 'DAPI', path=os.getcwd())


#######################################################################################################################