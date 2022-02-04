CXX=	clang++

# compiler flags.
CXXFLAGS+=	-std=c++11
#CXXFLAGS+=	-fopenmp -L/usr/local/lib -lomp
#CXXFLAGS+=	-Ofast -mtune=native -gfull
CXXFLAGS+=	-Oz -mtune=native -gfull
#CXXFLAGS+=	-Oz -mno-sse2 -mno-sse -mno-3dnow -mno-mmx -msoft-float -gfull
LDFLAGS+=	-lc++

#CXXFLAGS+=	-D_FLOAT_BITS_=32

clean:
	@rm -rf p0 p0-32 pc pc-32 p0-r p0-r-32
all:	p0 p0-32 pc pc-32 p0-r p0-r-32
p0:
	${CXX} ${CXXFLAGS} -static -o p0 p0.cc
p0-32:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=32 -o p0-32 p0.cc
p0-64:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=64 -o p0-64 p0.cc
p0-128:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=128 -o p0-128 p0.cc
pc:
	${CXX} ${CXXFLAGS} -static -o pc pc.cc
pc-32:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=32 -o pc-32 pc.cc
p0-r:
	${CXX} ${CXXFLAGS} -static -o p0-r p0-r.cc
p0-r-32:
	${CXX} ${CXXFLAGS} -static -D_FLOAT_BITS_=32 -o p0-r-32 p0-r.cc
