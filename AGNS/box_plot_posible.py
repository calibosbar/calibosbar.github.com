from matplotlib.pylab import *
import matplotlib.cm as cm
import matplotlib.patches as patches

min_val = 0
max_val = 50

my_cmap = cm.get_cmap('jet') # or any other one
norm = matplotlib.colors.Normalize(min_val, max_val) # the color maps work for [0, 1]

x_i = 15
color_i = my_cmap(norm(x_i)) # returns an rgba value

rect = patches.Rectangle((.5, .5), .25, .25, color=color_i) # make your rectangle

cmmapable = cm.ScalarMappable(norm, my_cmap)
cmmapable.set_array(range(min_val, max_val))

figure()
ax = gca()
ax.add_patch(rect)
colorbar(cmmapable)
show()
