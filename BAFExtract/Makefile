all: BAFExtract 

CC = g++
comp_flags = -c -Wall -O3
exec_name = BAFExtract

# Define pattern rule for building object files.
%.o: %.cpp
	${CC} ${comp_flags} $< -o $@

objs = src/main.o \
src/utils.o \
src/ansi_string.o \
src/genomics_coords.o \

BAFExtract: ${objs}
	${CC} -O3 -o bin/${exec_name} ${objs}

clean:
	rm -f *.o ${objs} bin/${exec_name} 
