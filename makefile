# makefile
# disable all warnings
CXXFLAGS += -w

simGRO: *.cpp *.h
	g++ GROgu.cpp -o GROgu -framework GLUT -framework OpenGL -framework Cocoa
	g++ simGRO.cpp -o simGRO
