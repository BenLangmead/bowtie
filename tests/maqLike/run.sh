#!/bin/sh
# Don't give it -m; we can't really use vmatch here
perl scripts/regression_test.pl \
        tests/maqLike/ref.fa \
        tests/maqLike/queries.mfa
