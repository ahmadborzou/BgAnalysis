DIR = $(shell pwd)
all:
	g++ main.cpp -LDelphes-3.0.10 -lDelphes -o main `root-config --cflags` `root-config --libs` 
