#!/bin/bash

# File: ~/headless.sh
# Control Mac display & sleep remotely via SSH
# Usage: ./headless.sh start | wake | stop

PIDFILE="/tmp/caffeinate_headless.pid"

start() {
    echo "[*] Starting headless mode..."
    # Keep system awake in background
    caffeinate -dimsu &
    echo $! > "$PIDFILE"
    sleep 1
    # Turn off the display
    pmset displaysleepnow
    echo "[+] Display off, system awake."
}

wake() {
    echo "[*] Waking display..."
    caffeinate -u -t 1
    echo "[+] Display should be on now."
}

stop() {
    echo "[*] Stopping headless mode..."
    if [[ -f "$PIDFILE" ]]; then
        kill "$(cat "$PIDFILE")" 2>/dev/null
        rm -f "$PIDFILE"
        echo "[+] Caffeinate stopped."
    else
        echo "[!] No caffeinate process found."
    fi
}

case "$1" in
    start) start ;;
    wake)  wake ;;
    stop)  stop ;;
    *) echo "Usage: $0 {start|wake|stop}" ;;
esac

