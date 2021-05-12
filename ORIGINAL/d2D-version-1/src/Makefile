.PHONY: all clean install

all: d2D camcoil cambeta camhelix camppii tamiola

d2D: predictor.c
	gcc  predictor.c -lm -O3 -g -Wall -Wextra -ansi -pedantic -o d2D-v1.x -I. -DLORENTZ

camcoil: coil_pred.c
	gcc  coil_pred.c -lm -O3 -g -Wall -Wextra -ansi -pedantic -o cam-coil.x

cambeta: beta_pred.c
	gcc  beta_pred.c -lm -O3 -g -Wall -Wextra -ansi -pedantic -o cam-beta.x

camhelix: helix_pred.c
	gcc  helix_pred.c -lm -O3 -g -Wall -Wextra -ansi -pedantic -o cam-helix.x

camppii: ppii_pred.c
	gcc  ppii_pred.c -lm -O3 -g -Wall -Wextra -ansi -pedantic -o cam-ppii.x

tamiola: tamiola.c
	gcc  tamiola.c -lm -O3 -g -Wall -Wextra -ansi -pedantic -o tamiola_coil.x

clean:
	rm *.x
	rm ../d2D.x

install: all
	cp *.x ../bin
