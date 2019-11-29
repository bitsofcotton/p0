CXX=	xcrun clang++
LD=	${CXX}

# compiler flags.
#CXXFLAGS+=	-fopenmp
CXXFLAGS+=	-std=c++11 -Ofast -g0 -mtune=native
LDFLAGS+=	-lc++

CLEANFILES= *.o *.dSYM p0

all:		p0
clean:
	@rm -rf ${CLEANFILES}

