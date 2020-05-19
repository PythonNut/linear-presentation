#!/bin/bash

rm *.aux
rm *.pdf
rm *.log

for i in *.tex; do
    [ -f "$i" ] || break

    echo $i

    # run pdflatex in quiet mode. See
    # https://tex.stackexchange.com/a/459470
    : | pdflatex -halt-on-error $i | grep '^!.*' -A200 --color=always
done
