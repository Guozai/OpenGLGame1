PROG = a3
CC = g++
CPPFLAGS = -g -Wall
LDFLAGS = -lSDL2 -lGLU -lGL -lm
OBJS = assignment3.o shaders.o

$(PROG) : $(OBJS)
		$(CC) $(OBJS) $(LDFLAGS) -o $(PROG)
assignment3.o :
		$(CC) $(CPPFLAGS) -c assignment3.c
shaders.o :
		$(CC) $(CPPFLAGS) -c shaders.c
clean :
		rm -f core $(PROG) $(OBJS)
