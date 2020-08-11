## custom cuboid mesh generation script for Blender 2.83

I have been toying with the idea of writing a python script to generate custom size/resolution cuboids for some time and last weekend, 
I finally got down to it. The challenge was to do it without using operators.

You can try it out by downloading the 'p_cuboid.py' file to your hard drive and opening it in Blender's Text Editor. There are two mesh 
generation functions available. The first one uses the BMesh module (but without operators) and the second one uses the Mesh module function 
from_pydata. There are also a few satelite functions which I have kept separate so they can be used in other mesh development scripts.

I hope you find it useful

Pan
