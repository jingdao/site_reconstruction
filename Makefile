CXX = g++
CC = gcc

all: mark_image match_image site_viewer bundler_viewer im_synth solve_pnp

mark_image: mark_image.c
	$(CC) -ggdb3 -std=gnu99 -o $@ $< -lSDL

match_image: match_image.c
	$(CC) -ggdb3 -std=gnu99 -o $@ $< -lSDL

site_viewer: site_viewer.cpp
	$(CXX) -ggdb3 -std=c++11 -o $@ $< -lSDL -lGL -lGLU

bundler_viewer: bundler_viewer.cpp
	$(CXX) -ggdb3 -std=c++11 -o $@ $< -lSDL -lGL -lGLU

im_synth: im_synth.cpp
	$(CXX) -ggdb3 -std=c++11 -o $@ $< -lSDL

solve_pnp: solve_pnp.cpp
	$(CXX) -ggdb3 -std=c++11 -o $@ $< -lopencv_calib3d -lopencv_core
