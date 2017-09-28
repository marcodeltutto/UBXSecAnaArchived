SRC_PATH=source
FILES=source/Spectrum.cpp source/Spectrum2D.cpp source/AnaTree.C
SOURCES = $(FILES:%.cpp=$(SRC_PATH)/%.cpp)

CFLAGS=-Wall

all:
	`root-config --cxx --cflags --glibs` ${CFLAGS} -I./include -o EventSelection EventSelection.cpp ${FILES} `root-config --glibs`
	`root-config --cxx --cflags --glibs` ${CFLAGS} -I./include -o DrawDataMC DrawDataMC.cpp ${FILES} `root-config --glibs`

draw:
	`root-config --cxx --cflags --glibs` ${CFLAGS} -I./include -o DrawDataMC DrawDataMC.cpp ${FILES} `root-config --glibs`
