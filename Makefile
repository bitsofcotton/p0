CXX=	clang++
LD=	${CXX}

# compiler flags.
#CXXFLAGS+=	-fopenmp
CXXFLAGS+=	-std=c++11 -Ofast -g0 -mtune=native
#CXXFLAGS+=	-std=c++11 -Ofast -gfull -mtune=native
CXXFLAGS+=	-D_FLOAT_BITS_=32
LDFLAGS+=	-lc++
#LDFLAGS+=	-lomp

CLEANFILES= *.o p0

all:		p0
clean:
	@rm -rf ${CLEANFILES}

