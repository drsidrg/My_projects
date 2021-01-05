import keras
import matplotlib.pyplot as plt
from keras.datasets import cifar10
from keras.models import Sequential
from keras.layers import Dense, Dropout, Flatten, Conv2D, MaxPooling2D
from pathlib import Path
from keras.datasets import cifar10
from keras.models import model_from_json
from keras.preprocessing import image
import matplotlib.pyplot as plt
import numpy as np
img = image.load_img("cat.png", target_size=(32, 32))

class_names = {
    0: "Plane",
    1: "Car",
    2: "Bird",
    3: "Cat",
    4: "Deer",
    5: "Dog",
    6: "Frog",
    7: "Horse",
    8: "Boat",
    9: "Truck"
}

(x_train, y_train), (x_test, y_test) = cifar10.load_data()

# for i in range(1000):
# image = x_train[i]
# image_class_number = y_train[i][0]
# image_class_number = class_names[image_class_number]


#    plt.imshow(image)
#   plt.title(image_class_number)
#  plt.show()

# Normalising data between 0 to 1
x_train = x_train.astype("float32")
x_train = x_train / 255
x_test = x_test.astype("float32")
x_test = x_test / 255

y_train = keras.utils.to_categorical(y_train, 10)
y_test = keras.utils.to_categorical(y_test, 10)

# create model and add layer
model = Sequential()

model.add(Conv2D(32, (3, 3), padding="same", activation="relu", input_shape=(32, 32, 3)))
model.add(Conv2D(32, (3, 3), activation="relu"))

model.add(MaxPooling2D(pool_size=(2, 2)))
model.add(Dropout(0.25))


model.add(Conv2D(64, (3, 3), padding="same", activation="relu"))
model.add(Conv2D(64, (3, 3), activation="relu"))

model.add(MaxPooling2D(pool_size=(2, 2)))
model.add(Dropout(0.25))

model.add(Flatten())
model.add(Dense(512, activation="relu"))
model.add(Dropout(0.50))

model.add(Dense(10, activation="softmax"))

model.compile(
    loss="categorical_crossentropy",
    optimizer="adam",
    metrics=["accuracy"]
)
model.summary()
##################################################
## now training will start

model.fit(
    x_train,
    y_train,
    batch_size=64,
    epochs=30,
    validation_data=(x_test, y_test),
    shuffle=True
)

# saving data and weights
model_structure = model.to_json()
f = Path("model_structure.json")
f.write_text(model_structure)

model.save_weights("model_weights.h5")

## making predictions -

class_label = [
    "Plane",
    "Car",
    "Bird",
    "Cat",
    "Deer",
    "Dog",
    "Frog",
    "Horse",
    "Boat",
    "Truck"
]

img = image.load_img("cat.png", target_size=(32, 32))
f = Path("model_structure.json")
model_structure = f.read_text()
model = model_from_json(model_structure)
model.load_weights("model_weights.h5")

img = image.load_img("cat.png", target_size=(32, 32))
image_to_test = image.img_to_array(img) / 255
# testing multiple images
list_of_images = np.expand_dims(image_to_test, axis=0)
results = model.predict(list_of_images)

single_result = results[0]
most_likely_class_index = int(np.argmax(single_result))
class_likelihood = single_result([most_likely_class_index])
print("This image is a {} - Likelihood: {2f}".format(class_label, class_likelihood))
