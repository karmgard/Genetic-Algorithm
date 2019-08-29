#
SHELL=/bin/bash

# Make sure the .dependencies file exists, otherwise the include at the bottom will choke
$(shell touch .dependencies)

SRC=main.cpp parameters.cpp population.cpp individual.cpp gnuplot.cpp
HDR=global.h individual.h parameters.h population.h utilities.h gnuplot.h
OBJ=$(subst .cpp,.o,${SRC})

TESTSRC=test.cpp
TESTOBJ=$(subst .cpp,.o,${TESTSRC})

LIBSRC=fitness.cpp utilities.cpp threadpool.cpp
LIBHDR=fitness.h utilities.h threadpool.h
LIBOBJ=$(subst .cpp,.o,${LIBSRC})
LIBBIN=libfitness.so

EXTRA=Makefile .dependencies README.md GPL-3.0 .gitignore

GCC_VERSION=`g++ -dumpversion`
ARCH=x86_64
OS=linux

INCLUDEPATHS= -I. -I ${HOME}/include \
	-I/usr/include/c++/${GCC_VERSION} \
	-I/usr/include/c++/${GCC_VERSION}/${ARCH}-linux-gnu \
	-I/usr/lib/gcc/${ARCH}-linux-gnu/${GCC_VERSION}/include \
	-I/usr/include/${ARCH}-linux-gnu \
	-I/usr/include/${ARCH}-linux-gnu/c++/${GCC_VERSION} \
	`pkg-config libgtop-2.0 --cflags`

LIBSEARCH=-L./ -L${HOME}/lib
LIBRARIES=-lm -lfitness `pkg-config libgtop-2.0 --libs`
DEBUG=0

ifeq (${DEBUG},1)
  CPUOPT=-g3 -Wall -Wunused -pg -fno-strict-aliasing -finline-functions -std=c++11
  LIBS=${LIBSEARCH} ${LIBRARIES} -pg
else
  CPUOPT=-Wall -Wunused -mhard-float -fno-strict-aliasing -finline-functions -std=c++11 -fno-stack-protector
  LIBS=${LIBSEARCH} ${LIBRARIES}
endif

ifneq (${OS},darwin)
  override CPUOPT+= $(patsubst %,%,-Wno-unused-but-set-variable)
  override CPUOPT+= $(patsubst %,%,-std=gnu++0x)
else
  override CPUOPT += $(patsubst %,%,-D__DARWIN__)
endif

CPUOPT+=-D__GXX_EXPERIMENTAL_CXX0X__

BIN=ga

CC=g++ $(CPUOPT) $(INCLUDEPATHS) 
LINK=g++ -o $(BIN) $(OBJ) $(LIBS)
LINKTEST=g++ -o test $(TESTOBJ)

all:	lib $(BIN)

clean:
	rm -f $(BIN) $(OBJ) $(LIBOBJ) $(LIBBIN) *~ *.bak .*.bak gmon.out test libfitness.so

tidy:
	rm -f $(BIN) $(OBJ) $(LIBOBJ) $(LIBBIN)

force:	tidy all

${BIN}:	dep $(OBJ)
	@echo ">>>>>>>>>>>> Linking <<<<<<<<<<<<<"
	$(LINK)
        ifeq (${DEBUG},0)
	  @strip ${BIN} ${LIBBIN}
        endif

test:	dep $(TESTOBJ)
	@echo ">>>>>>>>>>>> Linking <<<<<<<<<<<<<"
	$(LINKTEST)
        ifeq (${DEBUG},0)
	  @strip test
        endif
	@rm test.o

.cpp.o:
	@echo ">>>>>>>>>>>> Compiling $< -> $@ <<<<<<<<<<<<<"
	$(CC) -c $< -o $@

lib:
	@echo ">>>>>>>>>>>> Making Library <<<<<<<<<<<<<"
	$(CC) -fPIC $(CPUOPT) -c ${LIBSRC}
	$(CC) -shared -o ${LIBBIN} ${LIBOBJ} -pthread 

backup:
	@tar -zcf network.tar.gz $(SRC) $(HDR) $(LIBSRC) $(LIBHDR) $(EXTRA)
depend dep:
ifneq (${OS},darwin)
	makedepend  $(INCLUDEPATHS) $(SRC) -f .dependencies;
else
	makedepend -D__DARWIN__ $(INCLUDEPATHS) $(SRC) -f .dependencies;
endif

include .dependencies
# DO NOT DELETE
