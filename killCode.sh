#!/bin/bash

pids=$(ps -u $USER -o pid=,args= | grep '[x]stella' | awk '{print $1}')

if [ -z "$pids" ]; then
    echo "No xstella process found."
    exit 0
fi

echo "Killing xstella process(es):"
echo "$pids"

kill $pids
