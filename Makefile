OBJECTS=parser.tab.o lexer.yy.o main.o opengl.o helper.o RgbImage.o
SOURCES=*.cpp *.h
CFLAGS=-g
LDFLAGS=-lglut -lGL -lGLU

all: render

render: $(OBJECTS)
	g++ $(CFLAGS) -o render $^ $(LDFLAGS)

clean:
	rm *.o lexer.yy.cpp parser.tab.cpp parser.tab.hpp render

parser.tab.cpp parser.tab.hpp: parser.ypp
	bison -d parser.ypp

lexer.yy.cpp: lexer.lex
	flex -+ -olexer.yy.cpp lexer.lex

.cpp.o:
	g++ -std=c++11 -c $(CFLAGS) -o $@ $< $(LDFLAGS)

.PHONY: all clean
