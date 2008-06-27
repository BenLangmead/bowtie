# Exhaustively calculate difference-cover samples
#
# Note: difference-cover algorithm in diff_sample isn't good enough
# yet to calculate dc-128 even in a small number of days.  I don't
# recommend running this until the exhaustive algorithm can be made
# more clever.

# Simultaneously calculate difference covers for these periodicities
all: dc-128 dc-256 dc-512 dc-1024

dc-%:
	./diff_sample -e `echo $@ | sed s/.*-//` > $@
