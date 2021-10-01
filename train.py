from fastai.callbacks.tracker import SaveModelCallback
from fastai.metrics import accuracy
from fastai.vision import ImageDataBunch, cnn_learner, models
from fastai.vision.transform import get_transforms

path = 'images_ecoli/'
tfms = get_transforms(do_flip=True,flip_vert=False,max_rotate=0,max_zoom=1,max_lighting=None,max_warp=None,p_affine=0,p_lighting=0)
data = ImageDataBunch.from_folder(path,ds_tfms=tfms, bs=32)
learn = cnn_learner(data, models.resnet18, metrics=accuracy)

learn.fit_one_cycle(50, max_lr=0.0001,  callbacks=[SaveModelCallback(learn, every='improvement', monitor='accuracy', name='best_ecoli_model')])
