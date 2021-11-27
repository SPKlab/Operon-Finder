from functools import cache
from math import floor
from json import dump, loads
from os import makedirs
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


def main(genome_id: str, progress_bar) -> list[int]:
    # No need to update learn's data since input image passed manually
    data = ImageDataBunch.from_folder(path, test=f'test_operons/{genome_id}')
    total = sum(1 for _ in Path(path).joinpath(f'test_operons/{genome_id}').iterdir())
    
    makedirs('.json_files/operons/', exist_ok=True)
    predict_json = Path(f'.json_files/operons/{genome_id}.json')
    if predict_json.exists():
        predict_list = loads(predict_json.read_bytes())
    else:
        predict_list = []

        print("Enter predict")
        c = 0
        for i in range(len(data.test_ds)):
            filename = str(data.test_ds.items[i]).split('/')[-1]
            fsplit = filename.removesuffix('.jpg').split('_')
            c += 1
            progress_bar.progress(min(c/total*0.80 + .19, 0.99))
            is_operon = predictor(i, data)
            if is_operon:
                fname = fsplit[0]
                predict_list.append(int(fname.split(".")[-1]))
        print("End predict")
    
        with open(predict_json, 'w') as f:
            dump(predict_list, f)
        
    progress_bar.empty()
    return predict_list
        
    

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