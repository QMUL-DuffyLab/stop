SHELL = /bin/sh
FC    = mpifort
CURDIR = $(shell pwd)
SRCDIR = ${CURDIR}/src
LIB = ${CURDIR}
FLAGS = -std=f2018 -ffree-form -Wall -Werror -pedantic -fcheck=all
SOURCES := $(wildcard ${SRCDIR}/*.f90)
SOURCES = ${SRCDIR}/io.f90 ${SRCDIR}/lattice.f90 ${SRCDIR}/mc.f90 ${SRCDIR}/main.f90
TARGET = stop

OBJECTS = $(patsubst %.f90, %.o, $(SOURCES))

DEBUG = 0
ifeq (${DEBUG}, 1)
	FLAGS += -gdwarf-3
else
	FLAGS += -O2
endif

$(OBJECTS): %.o : %.f90
	$(FC) $(FLAGS) $(DFLAGS) -c -o $@ $<

$(TARGET): $(OBJECTS)
	$(FC) $(FLAGS) $(DFLAGS) -o $(TARGET) $(OBJECTS)

.PHONY: all clean check

all: $(TARGET)

clean:
	rm -f $(OBJECTS) $(TARGET) $(patsubst %.o, %.mod, $(OBJECTS))
