# Compiler name.
CC = gcc

# Compiler flags.
CFLAGS = -std=gnu99 -Wall

# Optimalisation flags
#CFLAGS += -O2 -march=native -mtune=native # default
#CFLAGS += -O2 -march=core-avx2	 # intel cpu family 6 gcc > 4.8.2
CFLAGS += -O2 -march=sandybridge # intel cpu family 6 gcc > 4.9.2, boneyard
#CFLAGS += -O2 -march=haswell     # intel cpu family 6 gcc > 4.9.2, inula 

# Linker flags (i.e. external libraries used by the project).
LDFLAGS = -lm

# Output binary.
binary = fccdcgen

# All necessary object files, equivilent of `ls *.c | sed 's/\.c/\.o/'`.
objfiles = $(patsubst %.c,%.o,$(wildcard *.c))

# Master (default) rule for building the binary, $@ is equiv. to $(binary).
$(binary): $(objfiles)
	$(CC) $(CFLAGS) $(objfiles) -o $@ $(LDFLAGS)

# Rules for compiling obj. files, with implicit recepies `$(CC) -o $@` omitted.
algebra.o   : algebra.h
initials.o  : config.h initials.h
utils.o     : config.h data.h mt19937.h mt19937_64.h algebra.h utils.h 
io.o        : data.h io.h 
structure.o : data.h utils.h algebra.h structure.h
fccdcgen.o  : config.h data.h globals.h initials.h io.h utils.h structure.h

# PHONY - targets that are only actions.
.PHONY: all clean debug mark
all: mark $(binary)
	mv $(binary) ..

# Clean the source directory of old object, backup files abd the binary.
clean:
	rm -rf $(objfiles) $(binary) *~

# Recompile the project without optimization but with debugging symbols.
debug: CFLAGS = -std=gnu99 -ggdb3 -O0 -Wall
debug: clean $(binary)

# Update checksum.h with sha1 sum's of all the source files.
mark:
	./mark.sh
# Modify source according to config.ini
configure:
	./config.sh
