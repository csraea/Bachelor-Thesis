import sys
import random
import hashlib
from PIL import Image
from PIL import ImageDraw
import os
import sys

READ_MAX = 65536

def save(outFileName, image):
    if not os.path.exists("out"):
        os.mkdir("out")
    os.chdir("out")
    image.save(os.path.basename(outFileName) + "-Hmatrix", 'PNG')
    os.chdir("..")

def drawLayer(imgSize, depth, mode):
    # Draw every block for a particular layer (blocks of a particular pixel size).   Returns an image of the finished layer.
    (imgW, imgH) = imgSize
    (w, h) = getBlockSize(imgSize, depth)
    if w < 1 or h < 1:
        return False
    if mode == 'v':    
        print("w, h: " + str(w) + ", " + str(h))
    layer = Image.new("RGB", imgSize, (0, 0, 0))
    draw = ImageDraw.Draw(layer)
    (x, y) = (0, 0)
    while y < imgH:
        while x < imgW:
            draw.rectangle([(x, y), (x + w - 1, y + h - 1)], fill=getRandomColor())
            x += w
        y += h
        x = 0
    return layer


def geHash(filename):
    # Compute hash of the file
    hash = hashlib.sha256()
    
    with open(filename, 'rb') as file:
        for lump in iter(lambda: file.read(READ_MAX), b''):
            hash.update(lump)
    return hash.hexdigest()


def getRandomColor():
    # Return a tuple of random color values
    rgbColor = []
    for i in range(3):
        rgbColor.append(random.randrange(0, 255))
    return tuple(rgbColor)


def getBlockSize(imgSize, depth):
    # Compute the block w and h for this layer
    (w, h) = imgSize
    while depth > 0:
        h = h / 4
        w = w / 2
        depth -= 1
        if depth < 1:
            break
        h = h / 2
        w = w / 4
        depth -= 1
    return (w, h)


def visualize(filename, mode, ssize):

    name, ext = os.path.splitext(filename)
    if ext != ".fasta" and ext != ".fa" and ext != ".fna":
        print("Error: specified file is not FASTA format")
        sys.exit()

    seed = geHash(filename)
    
    if mode == 'v':
        print("seed: " + seed)
    random.seed(seed)
    depth = 1
    layers = []

    size = (ssize, ssize)
    image = drawLayer(size, depth, mode)

    while 1:
        depth += 1
        layers.append(drawLayer(size, depth, mode))
        if layers[-1] == False:
            break

        opacity = 256 * (depth / 2) / 2 ** (depth - 1)  # can be changed

        if mode == 'v':
            print("opacity: " + str(int(opacity / 2.56)) + "%")
        maska = Image.new('L', size, int(opacity))
        image.paste(layers[-1], (0, 0), maska)

    save(name, image)