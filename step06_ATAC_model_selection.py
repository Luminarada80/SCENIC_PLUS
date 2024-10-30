import pycisTopic 
import os
import pickle

out_dir = "/gpfs/Home/haa5704/scenicplus/mESC_new_scenicplus/outs/"

models_filename = "/gpfs/Home/haa5704/scenicplus/mESC_new_scenicplus/outs/models.pkl"
cistopic_obj_filename = "/gpfs/Home/haa5704/scenicplus/mESC_new_scenicplus/outs/cistopic_obj.pkl"

with open(models_filename, 'rb') as file:
    models = pickle.load(file)

with open(cistopic_obj_filename, 'rb') as file:
    cistopic_obj = pickle.load(file)

from pycisTopic.lda_models import evaluate_models
model = evaluate_models(
    models,
    select_model = 40,
    return_model = True
)

cistopic_obj.add_LDA_model(model)

pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_obj.pkl"), "wb")
)


