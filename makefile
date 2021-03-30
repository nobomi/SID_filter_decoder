LINK_TARGET = sid_filt_decode.exe

OBJS =  \
 pocketfft/pocketfft.o \
 main.o

REBUILDABLES = $(OBJS) $(LINK_TARGET)

clean : 
	rm -f $(REBUILDABLES)

all : $(LINK_TARGET)
	echo All done

$(LINK_TARGET) : $(OBJS)
	gcc -o $@ $^

%.o : %.c
	gcc -o $@ -c $<

main.o : pocketfft/pocketfft.h
pocketfft.o : pocketfft/pocketfft.h


