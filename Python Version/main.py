
import numpy as np
import matplotlib.pyplot as plt
import os, sys
from scipy.io import loadmat
from sklearn.metrics import cohen_kappa_score
import matplotlib.patches as mpatches

import generate_nystrom_graph as model


#%%--------------------Available datasets-----------------------
#   sardania_dataset
#   omodeo_dataset
#   alaska_dataset
#   Madeirinha_dataset
#   katios_dataset
#   dique_dataset
#   SF_dataset
#   Wenchuan_dataset
#   canada_dataset
#   california_flood
#   contest_dataset
#-------------------------------------------------------------

#%% Read the data

data_path = os.path.join(os.path.dirname(os.getcwd()), 'Data') #this apply if you have the same structure as in the github

mat2 = loadmat(os.path.join(data_path, "sardania_dataset.mat"))
t1 = np.array(mat2["before"], dtype=np.float64)
t2 = np.array(mat2["after"], dtype=np.float64)
gt = np.array(mat2["gt"], dtype=bool)

#%% Process the data

n_samples = 50
g_nys = model.GBF_CD(t1, t2, n_samples)

change_map = g_nys.gbf_cd()

kappa = cohen_kappa_score(change_map.ravel(), gt.ravel())

labeled_map, cmap = g_nys.get_rgb_label_map(gt) 

plt.imshow(labeled_map, cmap = cmap, interpolation = 'none', vmin = 0, vmax = 3, aspect='auto')

# Creating legend with color box
pop_a = mpatches.Patch(color = '#5cff1c', label = 'Correct')
pop_b = mpatches.Patch(color = 'blue', label = 'MA')
pop_c = mpatches.Patch(color = 'red', label = 'FA')
plt.legend(handles = [pop_a, pop_b, pop_c], loc = 'upper center', ncol = 3, bbox_to_anchor=(0.5,1.1))

plt.axis('off')
plt.show()