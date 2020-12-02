# librairies de PRIMME
LIBP = -L./primme/ -lprimme

# librairies de JADAMILU_2.0
LIBJ = -L./JADAMILU_2.0/lib/INT32GNU/ -ljadamilu

# includes de PRIMME
INCP = -I./primme/PRIMMESRC/COMMONSRC/ 
# toutes les librairies
LIB = $(LIBP) $(LIBJ) -lm -lblas -llapack -lgfortran

COPT = -O3

default: main

clean: 
	rm *.o 
	rm main data.txt cmd.txt dataTemp.txt cmdTemp.txt

main: main.c main.h prob.o time.o primme.o residu.o plot.o euler_prog.o plot_temp.o
	cc $(COPT) $^ -o $@ $(LIB)

%.o: %.c
	cc $(COPT) -c $< -o $@ $(INCP)


