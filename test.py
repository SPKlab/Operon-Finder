from fastai.metrics import accuracy
from fastai.vision import ImageDataBunch, cnn_learner, models

path = 'images_custom/'

def main():
    # Load the data and predict Operon-pairs
    data = ImageDataBunch.from_folder(path, test='test_operons/')

    # Load the model
    learn = cnn_learner(data, models.resnet18, metrics=accuracy)
    learn = learn.load('best511')

    output = open('operon_pegs.txt', 'w')
    c = 0
    for i in data.test_ds:
        p = learn.predict(i[0])
        filename = str(data.test_ds.items[c]).split('/')[-1]
        fname = filename.split('_')[0]
        if str(p[0]) == 'Operons':
         output.write(fname + '\n')
        c += 1
    output.close()

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