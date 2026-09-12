CXX=	clang++
#CXX=	g++
#CXX=	eg++
#CXX=	c++
#CXX=	icpx
#CXX=	amdclang++
#CXX=	nvc++

# compiler flags.
CXXFLAGS+=	-Ofast -mtune=native -gfull
#CXXFLAGS+=	-O3 -mtune=native -g3
#CXXFLAGS+=	-O2 -g3
#CXXFLAGS+=	-O1 -gfull
# This doesn't work, we need operator >>, operator << with ongoing stdlibc++.
#CXXFLAGS+=	-I/usr/local/include -mlong-double-128
#CXXFLAGS+=	-Oz -mtune=native -gfull
#CXXFLAGS+=	-O1 -mtune=native -gfull
#CXXFLAGS+=	-O0 -mtune=native -gfull
#CXXFLAGS+=	-mno-sse2 -mno-sse -mno-3dnow -mno-mmx -msoft-float -fno-omit-frame-pointer
# XXX: this worse decreases multi thread performance.
#CXXFLAGS+=	-pg
#CXXFLAGS+=	--analyze
#CXXFLAGS+=      -D_LIBCPP_HARDENING_MODE_DEBUG
MPFLAGS=	-I/usr/local/include -L/usr/local/lib -lomp -fopenmp
#MPFLAGS=	-I/usr/local/include -L/usr/local/lib -lgomp -fopenmp
CXXFLAGS+=	-std=c++11
#CXXFLAGS+=	-std=gnu++98
LDFLAGS+=	-lc++ -L/usr/local/lib
#LDFLAGS+=	-lestdc++ -L/usr/local/lib
# Same as -mlong-double-128
#LDFLAGS+=	-lquadmath -lm

# lieonn.hh compile options
CXXFLAGS+=	-D_ARCFOUR_
# N.B. this specify after to sum up results
CXXFLAGS+=	-D_P_NOWALK_
# N.B. orthogonalize meaning
#CXXFLAGS+=	-D_P_SEPARATE_
# N.B. on disk datastream with cache = _P_ONDISK_ elements.
#CXXFLAGS+=	-D_P_ONDISK_=4194304
# N.B. GPGPU offloading
#CXXFLAGS+=	-I/usr/local/include -D_P_VULKAN_
#LDFLAGS+=	-L/usr/local/lib -lvulkan
# N.B. _SIMPLEALLOC_=align needs env VM_LIEONN=(mem usage MB).
# XXX: pred function needs huge memory with this because resize(0) doesn't
#      release the memory.
#CXXFLAGS+=	-D_SIMPLEALLOC_=64
# cf. src_alloc.c.diff
#CXXFLAGS+=	-I/usr/local/include -D_MIMALLOC_
#LDFLAGS+=	-L/usr/local/lib -lmimalloc
# N.B. _FLOAT_BITS_=bits for internal integer only calculation.
# N.B. this is NOT compatible with _P_VULKAN_.
#CXXFLAGS+=	-D_FLOAT_BITS_=32
#CXXFLAGS+=	-D_FLOAT_BITS_=64
#CXXFLAGS+=	-D_FLOAT_BITS_=128
# N.B. omit assertion, may have buggy but vasty speed up.
#CXXFLAGS+=	-D_OMIT_ASSERT_
# N.B. only use size_t and ssize_t for calculation on first order logic and
#      first order arithmetic operator + however memory reference needs
#      multiply sum operation on operator [].
#CXXFLAGS+=	-D_PERSISTENT_
# N.B. sed -e s/static\ inline//g | sed -e s/inline//g
#CXXFLAGS+=	-D_OLDCPP_ -ftemplate-depth-99

CLEANFILES= *.o p0 p0-32

clean:
	@rm -rf ${CLEANFILES}

all:	p0 p0-32

p0:
	${CXX} ${CXXFLAGS} -o p0.o -c p0.cc
	${CXX} p0.o ${LDFLAGS} -o p0
p0-32:
	${CXX} ${CXXFLAGS} -D_FLOAT_BITS_=32 -o p0-32.o -c p0.cc
	${CXX} p0-32.o ${LDFLAGS} -o p0-32

