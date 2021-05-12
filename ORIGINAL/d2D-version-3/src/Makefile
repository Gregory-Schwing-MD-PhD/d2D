.PHONY: all clean install

all: d2D
 
d2D: d2D-predictor.c
	gcc  d2D-predictor.c -lm -O3  -mtune=native -m64 -Wall -Wextra -ansi -pedantic -o d2D+.x -I. 

clean:
	rm d2D+.x

install: all
	cp d2D+.x ../bin

uninstall: clean 
	rm ../bin/d2D+.x
