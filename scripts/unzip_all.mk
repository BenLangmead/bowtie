# Unzips all of the .gz files in the current directory
# Use with -jN to do N in parallel

gzs=$(wildcard *.gz)
gunzs=$(patsubst %.gz,%,$(gzs))

all: $(gunzs)

%: %.gz
	zcat < $< > $@
