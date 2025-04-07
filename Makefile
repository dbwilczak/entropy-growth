# CAPD installation settings
CAPD = ../CAPD/build/bin
INCLUDE = `$(CAPD)/capd-config --cflags` -I./include
LIBS =  `$(CAPD)/capd-config --libs`
# compiler 
CXX = g++ -std=c++17 -O2 -s $(INCLUDE)

# programs to compile 
SOURCES = $(wildcard src/*.cpp)

PROGS = $(patsubst src/%.cpp,%,$(SOURCES)) 
DEPS = $(patsubst src/%.cpp,dep/%.d,$(SOURCES))

all: $(PROGS) $(DEPS)

include $(DEPS)

dep/%.d: src/%.cpp
	$(CXX) -MM -MT obj/$*.o $< > $@

%: src/%.cpp dep/%.d
	$(CXX) $(CXXFLAGS) $< -o $@ $(LIBS) -lpthread

clean:
	rm -f dep/*.d $(PROGS)

