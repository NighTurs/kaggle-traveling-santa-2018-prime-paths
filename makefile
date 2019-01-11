a.out : file_man.o global.o gpx.o main.o util_functions.o  	
	g++ -Wall -O3 -lm -lpthread file_man.o global.o gpx.o main.o util_functions.o -o a.out

file_man.o : file_man.cpp	
	g++ -Wall -O3 -o file_man.o -c file_man.cpp

global.o : global.cpp	
	g++ -Wall -O3 -o global.o -c global.cpp

gpx.o : gpx.cpp	
	g++ -Wall -O3 -o gpx.o -c gpx.cpp

main.o : main.cpp	
	g++ -Wall -O3 -lm -lpthread -o main.o -c main.cpp

util_functions.o : util_functions.cpp	
	g++ -Wall -O3 -o util_functions.o -c util_functions.cpp

clean :
	rm *.o
