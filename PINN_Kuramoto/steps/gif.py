import glob

steps = sorted(glob.glob(r"steps/*.png"), key=lambda x: int(x.split(".png")[0].split("\\")[1]))

from PIL import Image

frames = []
for step in steps:
    frame = Image.open(step)
    frames.append(frame)

frames[0].save("output.gif", save_all=True, append_images=frames[1:], loop=0, duration=50)