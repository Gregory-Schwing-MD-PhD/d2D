.PHONY: all clean install

all: d2D 

d2D: predictor.c
	gcc  predictor.c -lm -O3 -o d2D-v2.x -I. -DLORENTZ

clean:
	rm *.x
	rm ../d2D.x

install: all
	cp *.x ../bin
