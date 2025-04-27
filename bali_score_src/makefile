# Makefile — MSF/GCG‐only build of bali_score (no XML, no Expat)

CC       = cc
CFLAGS   = -std=c99 -O2 -Wall -DGCG -I.
LDFLAGS  = -lm

# Explicitly list source files (exclude readxml.c to drop XML support)
SRCS     = init.c util.c bali_score.c xmlstub.c
OBJS     = $(SRCS:.c=.o)

.PHONY: all clean

all: bali_score

bali_score: $(OBJS)
	$(CC) -o $@ $(OBJS) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f bali_score $(OBJS)
