# Compiler name.
CC = gcc

# Compiler flags.
CFLAGS = -std=gnu99 -Wall

# Optimalisation flags
CFLAGS += -O2 -march=native -mtune=native

# enable profiling
# CFLAGS += -pg

# Linker flags (i.e. external libraries used by the project).
LDFLAGS = -lm

# Output binary.
binary = xcsg

# All necessary object files, equivilent of `ls *.c | sed 's/\.c/\.o/'`.
objfiles = $(patsubst %.c,%.o,$(wildcard *.c))

# Master (default) rule for building the binary, $@ is equiv. to $(binary).
$(binary): $(objfiles)
	$(CC) $(CFLAGS) $(objfiles) -o $@ $(LDFLAGS)

# Rules for compiling obj. files, with implicit recepies `$(CC) -o $@` omitted.
algebra.o   : algebra.h
io.o        : config.h data.h io.h 
utils.o     : config.h data.h mt19937.h mt19937_64.h algebra.h utils.h 
initials.o  : config.h checksum.h data.h io.h utils.h initials.h
structure.o : data.h utils.h algebra.h structure.h
terminators.o : config.h data.h
xcsg.o  : config.h data.h globals.h initials.h io.h utils.h structure.h terminators.h

# PHONY - targets that are only actions.
.PHONY: all clean debug mark
all: configure mark $(binary)
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
