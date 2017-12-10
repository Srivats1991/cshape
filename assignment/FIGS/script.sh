#!/bin/bash

cd alpha
for i in *.bmp; do convert $i ${i/.bmp}.png; done
cd ..

cd convex_hull
for i in *.bmp; do convert $i ${i/.bmp}.png; done
cd ..

cd crust
for i in *.bmp; do convert $i ${i/.bmp}.png; done
cd ..

cd delaunay
for i in *.bmp; do convert $i ${i/.bmp}.png; done
cd ..

cd points
for i in *.bmp; do convert $i ${i/.bmp}.png; done
cd ..
