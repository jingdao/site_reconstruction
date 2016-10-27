CXX = g++
CC = gcc

all: mark_image match_image site_viewer bundler_viewer im_synth solve_pnp

mark_image: mark_image.c
	$(CC) -ggdb3 -std=gnu99 -o $@ $< -lSDL

match_image: match_image.c
	$(CC) -ggdb3 -std=gnu99 -o $@ $< -lSDL

site_viewer: site_viewer.cpp
	$(CXX) -ggdb3 -std=c++11 -o $@ $< -lSDL -lGL -lGLU

site_viewer_2d: site_viewer_2d.cpp
	$(CXX) -ggdb3 -std=c++11 -o $@ $< -lSDL -lGL -lGLU -lfreetype -I/usr/include/freetype2

bundler_viewer: bundler_viewer.cpp
	$(CXX) -ggdb3 -std=c++11 -o $@ $< -lSDL -lGL -lGLU

im_synth: im_synth.cpp
	$(CXX) -ggdb3 -std=c++11 -o $@ $< -lSDL

solve_pnp: solve_pnp.cpp
	$(CXX) -ggdb3 -std=c++11 -o $@ $< -lopencv_calib3d -lopencv_core

solve_pnp_2d: solve_pnp_2d.cpp
	$(CXX) -ggdb3 -std=c++11 -o $@ $<

imsegment: imsegment.cpp
	$(CXX) -ggdb3 -std=c++11 -o $@ $< -lSDL

imlabel: imlabel.cpp
	$(CXX) -ggdb3 -std=c++11 -o $@ $< -lSDL

solve_accuracy: solve_accuracy.cpp
	$(CXX) -ggdb3 -std=c++11 -o $@ $<

vocab_tree: vocab_tree.cpp
	$(CXX) -ggdb3 -std=c++11 -o $@ $< -llapack -lblas 

