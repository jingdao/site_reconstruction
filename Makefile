CXX = g++
CC = gcc

all: mark_image match_image

mark_image: mark_image.c
	$(CC) -ggdb3 -std=gnu99 -o $@ $< -lSDL

match_image: match_image.c
	$(CC) -ggdb3 -std=gnu99 -o $@ $< -lSDL
