#!/bin/sh

# Copyright 2018 GRAIL, Inc. All rights reserved.
# Use of this source code is governed by the Apache 2.0
# license that can be found in the LICENSE file.

out=results.tsv

cat 170501*/*.conta.tsv | grep "^conta_version" | head -n 1 > header.tsv
cat 170501*/*.conta.tsv | grep -v "^conta_version" > all.tsv
cat header.tsv all.tsv > $out
rm header.tsv all.tsv
