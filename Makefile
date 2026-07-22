# Compilation target.
TARGET = massey

# Directories
SRCDIR = src
BUILDDIR = build
HEADERDIR = headers


# Compiler.
CC = clang

# Compiler flags. PARI should point to the installation prefix.
CPPFLAGS = -I. -I$(HEADERDIR) -I$(PARI)/include
CFLAGS = -O3 -Wall -Wunused-function -fno-strict-aliasing -fomit-frame-pointer -pipe -march=native -pthread -g
LDFLAGS = -Wl,-O3 -Wl,-rpath,$(PARI)/lib -pthread
LIBS = -L$(PARI)/lib -lpari -lm
#STATIC = -static

# Source files
SOURCES = $(wildcard $(SRCDIR)/*.c)
OBJECTS = $(SOURCES:$(SRCDIR)/%.c=$(BUILDDIR)/%.o)

# Default target
$(BUILDDIR)/$(TARGET): $(OBJECTS) | $(BUILDDIR)
	$(CC) -o $@ $(OBJECTS) $(LDFLAGS) $(LIBS)

# Compile source files to object files
$(BUILDDIR)/%.o: $(SRCDIR)/%.c | $(BUILDDIR)
	$(CC) -c -o $@ $< $(CFLAGS) $(CPPFLAGS)

# Create build directory if it doesn't exist
$(BUILDDIR):
	mkdir -p $(BUILDDIR)

# Phony targets
.PHONY: clean all

all: $(BUILDDIR)/$(TARGET)

clean:
	-$(RM) -r $(BUILDDIR)
	-$(RM) *.o $(TARGET)
