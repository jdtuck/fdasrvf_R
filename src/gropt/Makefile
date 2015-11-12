.PHONY : clean

ERR = $(shell which icpc>/dev/null; echo $$?)
ifeq "$(ERR)" "0"
    CXX = icpc
	CFLAGS =  -fPIC -O3 -g # C flags
	LDFLAGS = -shared # linking flags
	LIBS = -lmkl_rt -lpthread
else
	CFLAGS =  -fPIC -std=c++11 -O3 -g  # C flags
	LDFLAGS = -shared # linking flags
	LIBS = -lblas -llapack
endif


LIB=gropt
SUFFIX=so
SOURCES = $(shell echo src/*.cpp)
HEADERS = $(shell echo incl/*.h)
OBJECTS = $(SOURCES:.cpp=.o)
INC = -Iincl/
TARGET=$(LIB).$(SUFFIX)

all: $(TARGET)

clean:
	rm -f $(OBJECTS) $(TARGET)

src/%.o: src/%.cpp
	$(CXX) $(CFLAGS) $(INC) -c -o $@ $<

$(TARGET) : $(OBJECTS)
	$(CXX) $(LDFLAGS) -o $@ $^ $(LIBS)
