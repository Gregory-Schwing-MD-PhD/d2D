.PHONY: all clean install

all: d2D camcoil cambeta camhelix camppii

d2D: d2D-predictor.c
	gcc  d2D-predictor.c -lm -O3  -mtune=native -m64 -Wall -Wextra -ansi -pedantic -o d2D+.x -I. 

camcoil: coil_pred.c
	gcc  coil_pred.c -lm -O3 -Wall -Wextra -ansi -pedantic -o cam-coil.x

cambeta: beta_pred.c
	gcc  beta_pred.c -lm -O3 -Wall -Wextra -ansi -pedantic -o cam-beta.x

camhelix: helix_pred.c
	gcc  helix_pred.c -lm -O3 -Wall -Wextra -ansi -pedantic -o cam-helix.x

camppii: ppii_pred.c
	gcc  ppii_pred.c -lm -O3 -Wall -Wextra -ansi -pedantic -o cam-ppii.x

clean:
	rm *.x
	rm ../d2D+.x

install: all
	cp *.x ../bin
