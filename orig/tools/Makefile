ROOT = /home/alastair/multiuser
INCLUDE = -I/home/alastair/multiuser/include
CFLAGS = -g -Wall -Wextra $(INCLUDE)
LIBS = -lm
CC = g++

#gamma_objects = main.o $(ROOT)/src/symerr.o

dbm_objects = dbmhz_to_watts.o

watts_objects = watts_to_dbmhz.o

#gamma: $(gamma_objects) main.c $(ROOT)/src/symerr.c
#	cc -g -o gamma $(gamma_objects) $(LIBS)
dbm: $(dbm_objects) dbmhz_to_watts.c
	g++ -g -o dbmhz_to_watts $(dbm_objects) $(LIBS) 
watts: $(watts_objects) watts_to_dbmhz.c
	g++ -g -o watts_to_dbmhz $(watts_objects) $(LIBS)
clean:
