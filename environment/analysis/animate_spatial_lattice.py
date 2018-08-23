import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

lattice_file = open("out/manual/environment.txt")
locations_file = open("out/manual/locations.txt")

lattice_data = lattice_file.read()
lattice_data = lattice_data.split('\n')
lattice_data_list = []
for row in lattice_data[:-1]:
    rows = row.replace('[', '').split('],')
    list = [map(float, s.replace(']', '').split(',')) for s in rows]
    lattice_data_list.append(list)

n_time_steps = len(lattice_data_list)

fig = plt.figure()

ims = []
for i in range(n_time_steps):
    im = plt.imshow(np.asarray(lattice_data_list[i]), animated=True)
    ims.append([im])

ani = animation.ArtistAnimation(fig, ims, interval=50, blit=True, repeat_delay=1000)
ani.save('lattice_animation.html')
plt.show()
