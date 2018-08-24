import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

import shutil
import os


os.chdir(os.path.expanduser('out/manual/'))

if os.path.exists("lattice_animation_frames") and os.path.isdir("lattice_animation_frames"):
    shutil.rmtree("lattice_animation_frames")
os.makedirs("lattice_animation_frames")
if os.path.exists("lattice_animation.html"):
    os.remove("lattice_animation.html")

lattice_file = open("environment.txt")
locations_file = open("locations.txt")

lattice_data = lattice_file.read()
lattice_data = lattice_data.split('\n')
lattice_data_list = []
for row in lattice_data[:-1]:
    rows = row.replace('[', '').split('],')
    list = [map(float, s.replace(']', '').split(',')) for s in rows]
    lattice_data_list.append(list)

locations_data = locations_file.read()
locations_data = locations_data.split('\n')
locations_data_list = []
for row in locations_data[:-1]:
    rows = row.replace('[', '').split('],')
    list = [map(float, s.replace(']', '').split(',')) for s in rows]
    locations_data_list.append(list)

n_bins = len(lattice_data_list[0])

n_time_steps = len(lattice_data_list)


fig = plt.figure()

# ims = []
for i in range(n_time_steps):
    lattice = lattice_data_list[i]
    locations_x = [location[1]*n_bins - 0.5 for location in locations_data_list[i]]
    locations_y = [location[0]*n_bins - 0.5 for location in locations_data_list[i]]

    # clear figure
    plt.clf()

    # plot lattice
    im = plt.imshow(np.asarray(lattice)) #, animated=True)
    plt.axis('off')

    # add locations
    plt.scatter(locations_x, locations_y) #, animated=True)

    plt.savefig('lattice_animation_frames/' + str(i) + '.png')



#
# images = []
# for image in os.listdir('lattice_animation_frames/*'):
#     images.append(Image.open(image))
#
#
# ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True, repeat_delay=1000)


#     ims.append([im])
#
#
# ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True, repeat_delay=1000)
# ani.save('lattice_animation.html')
# plt.show()
