from functools import cache
from math import floor
from json import dump, loads
from os import makedirs
from pathlib import Path
from fastai.metrics import accuracy
from fastai.vision import ImageDataBunch, cnn_learner, models

path = 'images_custom'

# Load the data and predict Operon-pairs

# Load the model
def get_learn(path):
    makedirs(f"{path}/train", exist_ok=True)
    Path(f"{path}/train/dummy.jpg").touch(exist_ok=True)
    _data = ImageDataBunch.from_folder('images_ecoli', test='test_operons/')
    return cnn_learner(_data, models.resnet18, metrics=accuracy).load('best511')

learn = get_learn(path)

def predictor(idx: int, data: ImageDataBunch) -> tuple[float]:
    p = learn.predict(data.test_ds[idx][0])
    confidence = p[-1][1].item()
    return confidence


def main(genome_id: str, progress_clb) -> dict[str, float]:
    # No need to update learn's data since input image passed manually
    data = ImageDataBunch.from_folder(path, test=f'test_operons/{genome_id}')
    total = sum(1 for _ in Path(path).joinpath(f'test_operons/{genome_id}').iterdir())
    
    predict_list = {}

    print("Enter predict")
    c = 0
    for i in range(len(data.test_ds)):
        filename = str(data.test_ds.items[i]).split('/')[-1]
        fsplit = filename.removesuffix('.jpg').split('_')
        c += 1
        progress_clb(min(0.65 + c/total*0.35, 0.99))
        confidence = predictor(i, data)
        fname = fsplit[0]
        predict_list[fname.split(".")[-1]] = round(confidence, 2)
    print("End predict")
    progress_clb(1.0)
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
