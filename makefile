CC     = gcc
CFLAGS = 

INCDIRS = -I/usr/lib/x86_64-linux-gnu/glib-2.0/include -I/usr/include/glib-2.0

LIBS = -lm -lglib-2.0 

SRCS = main.c thpool.c 
OBJS = $(SRCS:.c=.o)
EXE  = Optimization 

OPTIFLAG = -O0

DBGDIR = debug
DBGEXE = $(DBGDIR)/$(EXE)
DBGOBJS = $(addprefix $(DBGDIR)/, $(OBJS))
DBGCFLAGS = -g -DDEBUG

RELDIR = release
RELEXE = $(RELDIR)/$(EXE)
RELOBJS = $(addprefix $(RELDIR)/, $(OBJS))
RELCFLAGS = -DNDEBUG

TESTDIR = test
TESTEXE = $(TESTDIR)/$(EXE)
TESTOBJS = $(addprefix $(TESTDIR)/, $(OBJS))

.PHONY: all clean debug prep release test

all: prep release

####  Target defines
release: OPTIFLAG = -O3
relwdbg: RELCFLAGS += -g
relwdbg:  OPTIFLAG = -O3
test: RELCFLAGS += -DTEST 
test: OPTIFLAG = -O3
fast: OPTIFLAG = -Ofast
debug: OPTIFLAG = -O0
####

debug: prep $(DBGEXE)

$(DBGEXE): $(DBGOBJS)
	$(CC) $(CFLAGS) $(OPTIFLAG) $(DBGCFLAGS) -o $(DBGEXE) $^ $(LIBS) 

$(DBGDIR)/%.o: %.c
	$(CC) -c  $(CFLAGS) $(OPTIFLAG) $(DBGCFLAGS) $(INCDIRS) -o $@ $<

test: prep $(TESTEXE)

relwdbg: prep $(RELEXE)

release: prep $(RELEXE)

$(RELEXE): $(RELOBJS)
	$(CC) $(CFLAGS) $(RELCFLAGS) $(OPTIFLAG) -o $(RELEXE) $^ $(LIBS)

$(RELDIR)/%.o: %.c
	$(CC) -c $(CFLAGS) $(RELCFLAGS) $(OPTIFLAG) $(INCDIRS) -o $@ $<

$(TESTEXE): $(TESTOBJS)
	$(CC) $(CFLAGS) $(RELCFLAGS) $(OPTIFLAG) -o $(TESTEXE) $^ $(LIBS)

$(TESTDIR)/%.o: %.c
	$(CC) -c $(CFLAGS) $(RELCFLAGS) $(OPTIFLAG) $(INCDIRS) -o $@ $<

prep:
	@mkdir -p $(DBGDIR) $(RELDIR) $(TESTDIR)

remake: clean all

clean:
	rm -f $(RELEXE) $(RELOBJS) $(DBGEXE) $(DBGOBJS) $(TESTEXE) $(TESTOBJS)
