#!/usr/bin/env bash

(zcat $1 | head -100 | grep ^#; zcat $1| grep -v ^# | sort -k1,1d -k2,2n;) | bgzip -c > $2 