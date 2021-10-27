from functools import cache
from math import floor
from json import dump, loads
from pathlib import Path
from fastai.metrics import accuracy
from fastai.vision import ImageDataBunch, cnn_learner, models
import streamlit as st

path = 'images_custom/'

# Load the data and predict Operon-pairs
_data = ImageDataBunch.from_folder(path, test='test_operons/')

# Load the model
learn = cnn_learner(_data, models.resnet18, metrics=accuracy)
learn = learn.load('best511')

def predictor(idx: int, data: ImageDataBunch):
    p = learn.predict(data.test_ds[idx][0])
    is_operon = str(p[0]) == 'Operons'
    return is_operon


def main(genes: set, progress_bar):
    # No need to update learn's data since input image passed manually
    data = ImageDataBunch.from_folder(path, test='test_operons/')
    
    predict_json = Path('.json_files/predict.json')
    predict_dict = loads(predict_json.read_bytes()) if predict_json.exists() else {}

    output = open('operon_pegs.txt', 'w')
    c = 0
    print("Enter predict")
    c = 0
    for i in range(len(data.test_ds)):
        filename = str(data.test_ds.items[i]).split('/')[-1]
        fsplit = filename.removesuffix('.jpg').split('_')
        if genes.issuperset(fsplit):
            c += 1
            progress_bar.progress(min(c/len(genes)*0.80 + .19, 0.99))
            if filename not in predict_dict:
                predict_dict[filename] = predictor(i, data)
            is_operon = predict_dict[filename]
            if is_operon:
                fname = fsplit[0]
                output.write(fname + '\n')
    print("End predict")
    output.close()
    
    with open(predict_json, 'w') as f:
        dump(predict_dict, f)
        
    progress_bar.empty()
        
    

    # Load the data and predict non-Operon-pairs
    if False:
        data = ImageDataBunch.from_folder(path, test='test_noperons/')
        output = open('noperon_pegs.txt', 'w')
        c = 0
        for i in data.test_ds:
         p = learn.predict(i[0])
         filename = str(data.test_ds.items[c]).split('/')[-1]
         fname = filename.split('_')[0]
         if str(p[0]) == 'Noperons':
          output.write(fname + '\n')
         c += 1
        output.close()