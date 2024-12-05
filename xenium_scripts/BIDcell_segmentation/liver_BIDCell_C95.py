## run a full sample attempt
from bidcell import BIDCellModel
model = BIDCellModel("/cluster/home/t117652uhn/xenium_liver/scripts/liver_BIDCell_C95.yaml")
model.preprocess()


# import os
# import numpy as np
# import h5py

# folder_path = "/cluster/projects/macparland/RE/xenium_liver/C94_2/expr_maps/expr_maps_input_patches_48x48_shift_0/"

# # List all files in the directory
# files = os.listdir(folder_path)

# # Filter out only HDF5 files
# hdf5_files = [file for file in files if file.endswith('.hdf5')]

# for file_name in hdf5_files:
#     patch_fp = os.path.join(folder_path, file_name)
#     print(file_name)
#     h5f = h5py.File(patch_fp, "r")
#     expr = h5f["data"][:].astype(np.float64)
#     h5f.close()
    
#     if expr.shape[0] != 48:
#         # Delete the file
#         os.remove(patch_fp)
#         print("Deleted file:", file_name)


model.train()
model.predict()

