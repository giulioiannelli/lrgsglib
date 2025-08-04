import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from mpl_toolkits.axes_grid1 import make_axes_locatable

__all__ = ["make_animation_fromFrames"]

def make_animation_fromFrames(frames, savename="output.mp4", fps=10, dpi=200):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # I like to position my colorbars this way, but you don't have to
    div = make_axes_locatable(ax)
    cax = div.append_axes("right", "5%", "5%")
    cv0 = frames[0]
    im = ax.imshow(cv0)  # Here make an AxesImage rather than contour
    cb = fig.colorbar(im, cax=cax)

    # tx = ax.set_title('Frame 0')
    def animate(i):
        arr = frames[i]
        vmax = np.max(arr)
        vmin = np.min(arr)
        im.set_data(arr)
        im.set_clim(vmin, vmax)
        # tx.set_text('Frame {0}'.format(i))
        # In this version you don't have to do anything to the colorbar,
        # it updates itself when the mappable it watches (im) changes

    print("# of frames: ", len(frames))
    ani = animation.FuncAnimation(fig, animate, frames=len(frames))
    fig.tight_layout()
    writervideo = animation.FFMpegWriter(fps=fps)
    ani.save(savename, writer=writervideo, dpi=dpi)