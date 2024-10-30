import pycisTopic
import os
import pickle

out_dir = "/gpfs/Home/haa5704/scenicplus/mESC_new_scenicplus/outs/"

# Path to your .pkl file
file_path = '/gpfs/Home/haa5704/scenicplus/mESC_new_scenicplus/outs/cistopic_obj.pkl'

# Load the pickle file
with open(file_path, 'rb') as file:
    cistopic_obj = pickle.load(file)

os.environ['MALLET_MEMORY'] = '200G'
from pycisTopic.lda_models import run_cgs_models_mallet
# Configure path Mallet
mallet_path="/gpfs/Home/haa5704/scenicplus/Mallet-202108/bin/mallet"
# Run models
models=run_cgs_models_mallet(
    cistopic_obj,
    n_topics=[2, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50],
    n_cpu=12,
    n_iter=500,
    random_state=555,
    alpha=50,
    alpha_by_topic=True,
    eta=0.1,
    eta_by_topic=False,
    tmp_path="/gpfs/Home/haa5704/scenicplus/tmp_mESC/",
    save_path="/gpfs/Home/haa5704/scenicplus/tmp_mESC/",
    mallet_path=mallet_path,
)

pickle.dump(
    models,
    open(os.path.join(out_dir, "models.pkl"), "wb")
)
