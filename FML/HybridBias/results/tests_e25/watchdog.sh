#!/bin/bash
# Memory watchdog for one run: records the node's available memory and extra swap every second, and kills the user's
# nbody processes if available memory drops below 15 GB or more than 512 MB of swap is used. Usage: watchdog.sh <log>
LOG=$1
swap0=$(awk '/SwapTotal/{t=$2}/SwapFree/{f=$2}END{print int((t-f)/1024)}' /proc/meminfo)
min_avail=999999; max_swap=0
while true; do
    read avail swap < <(awk -v s0=$swap0 '/MemAvailable/{a=$2}/SwapTotal/{t=$2}/SwapFree/{f=$2}END{print int(a/1048576), int((t-f)/1024)-s0}' /proc/meminfo)
    (( avail < min_avail )) && min_avail=$avail
    (( swap > max_swap )) && max_swap=$swap
    echo "min_available_GB $min_avail max_extra_swap_MB $max_swap" > $LOG
    if (( avail < 15 || swap > 512 )); then
        echo "KILLED at available ${avail} GB, extra swap ${swap} MB" >> $LOG
        pkill -u $USER -f "/nbody_float "
        sleep 5
    fi
    sleep 1
done
