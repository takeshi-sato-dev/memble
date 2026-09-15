#!/usr/bin/env bash
#
# curve2_run.sh : run the six systems of the second curve over several days.
#
# The machine holds three devices and only two of them are compute cards. The
# third is a display adapter, and GROMACS distributes jobs over all three when
# it is left to choose, which puts a third of the work on a card that returns
# about fifteen nanoseconds a day. This script hides the display adapter with
# CUDA_VISIBLE_DEVICES, gives each job one card and a fixed number of threads,
# and runs the six in pairs, one pair at a time.
#
# Six systems at four hours each on two cards is twelve hours, which is longer
# than a working day. The script is therefore made to be started in the morning
# and stopped in the evening, every day, with the same two commands:
#
#   morning   STOP_AT=18:00 bash curve2_run.sh run
#   evening   nothing, if STOP_AT was given; otherwise bash curve2_run.sh stop
#
# Every job resumes from its own checkpoint and every system that has already
# finished is skipped, so the morning command may be given on any number of
# days until `status` reports all six finished.
#
# Usage:
#   bash curve2_run.sh run            the pairs that are not yet finished
#   STOP_AT=18:00 bash curve2_run.sh run     the same, stopping at 18:00
#   bash curve2_run.sh pair 2         only the second pair
#   bash curve2_run.sh status         how far each system is
#   bash curve2_run.sh stop           stop every job of this curve
#
set -uo pipefail

export CUDA_DEVICE_ORDER=PCI_BUS_ID
ROOT=${ROOT:-$HOME/counts_run/curve2}
GMX=${GMX:-$(command -v gmx || command -v gmx_mpi)}
DEFF=${DEFF:-step7_production}
NT=${NT:-16}
CARDS=${CARDS:-"0 2"}
NSTEPS=${NSTEPS:-12500000}
PAIRS=("-20 -12" "-6 0" "6 12")
FLAG="$ROOT/.stopped"

work_of() { echo "$ROOT/c2d${1}/c2d${1}_work"; }

last_step() {
  local w; w=$(work_of "$1")
  local s=""
  for f in "$w/prod_run.log" "$w/log_step7.txt"; do
    [ -f "$f" ] || continue
    local v
    v=$(tr '\r' '\n' < "$f" 2>/dev/null | grep -a "^step " | tail -1 | awk '{print $2}')
    [ -n "$v" ] && s=$v
  done
  echo "${s:-0}"
}

finished() {
  local w; w=$(work_of "$1")
  grep -qa "Finished mdrun" "$w/prod_run.log" 2>/dev/null && return 0
  grep -qa "Finished mdrun" "$w/log_step7.txt" 2>/dev/null && return 0
  local s; s=$(last_step "$1")
  [ "${s%%.*}" -ge "$NSTEPS" ] 2>/dev/null && return 0
  return 1
}

start_one() {
  local d=$1 card=$2 slot=$3
  local w; w=$(work_of "$d")
  if [ ! -f "$w/$DEFF.tpr" ]; then
    echo "  d$d  no $DEFF.tpr, skipped"
    return 1
  fi
  if finished "$d"; then
    echo "  d$d  already finished, skipped"
    return 1
  fi
  local cpi=""
  [ -f "$w/$DEFF.cpt" ] && cpi="-cpi $DEFF.cpt -append"
  (
    cd "$w" || exit 1
    CUDA_VISIBLE_DEVICES=$card \
    $GMX mdrun -deffnm "$DEFF" $cpi -v \
        -ntmpi 1 -ntomp "$NT" -pin on -pinoffset $((slot * NT)) -pinstride 1 \
        -nb gpu -gpu_id 0 \
        >> prod_run.log 2>&1
  ) &
  echo "  d$d  started on card $card, threads $NT, pid $!"
  return 0
}

run_pair() {
  local i=$1
  local -a two; read -r -a two <<< "${PAIRS[$i]}"
  local -a cards; read -r -a cards <<< "$CARDS"
  if finished "${two[0]}" && finished "${two[1]}"; then
    echo "pair $((i + 1)): d${two[0]} and d${two[1]} are both finished, skipped"
    return 0
  fi
  echo "pair $((i + 1)): d${two[0]} and d${two[1]}   $(date '+%F %H:%M')"
  local slot=0 started=0
  for j in 0 1; do
    start_one "${two[$j]}" "${cards[$j]}" "$slot" && started=$((started + 1))
    slot=$((slot + 1))
    sleep 3
  done
  [ "$started" -eq 0 ] && return 0
  wait
  echo "pair $((i + 1)) ended   $(date '+%F %H:%M')"
}

arm_clock() {
  [ -n "${STOP_AT:-}" ] || return 0
  local target now
  now=$(date +%s)
  target=$(date -d "today $STOP_AT" +%s 2>/dev/null) || {
    echo "STOP_AT is not a time this shell understands: $STOP_AT"; exit 1; }
  [ "$target" -le "$now" ] && target=$(date -d "tomorrow $STOP_AT" +%s)
  echo "will stop at $(date -d "@$target" '+%F %H:%M'), $(( (target - now) / 60 )) minutes from now"
  (
    sleep $((target - now))
    touch "$FLAG"
    pkill -f "mdrun -deffnm $DEFF"
  ) &
  WATCH=$!
}

case "${1:-}" in
  run)
    [ -n "$GMX" ] || { echo "STOP: gmx not found"; exit 1; }
    mkdir -p "$ROOT"
    rm -f "$FLAG"
    WATCH=""
    arm_clock
    for i in "${!PAIRS[@]}"; do
      run_pair "$i"
      if [ -f "$FLAG" ]; then
        echo ""
        echo "stopped by the clock. Give the same command tomorrow; every system"
        echo "resumes from its checkpoint and the finished ones are skipped."
        break
      fi
    done
    [ -n "$WATCH" ] && kill "$WATCH" 2>/dev/null
    rm -f "$FLAG"
    echo ""
    bash "$0" status
    ;;
  pair)
    [ -n "$GMX" ] || { echo "STOP: gmx not found"; exit 1; }
    rm -f "$FLAG"
    WATCH=""
    arm_clock
    run_pair $(( ${2:-1} - 1 ))
    [ -n "$WATCH" ] && kill "$WATCH" 2>/dev/null
    rm -f "$FLAG"
    ;;
  status)
    left=0
    for d in -20 -12 -6 0 6 12; do
      s=$(last_step "$d")
      if finished "$d"; then
        printf "d%-4s %10s / %s   finished\n" "$d" "$s" "$NSTEPS"
      else
        left=$((left + 1))
        printf "d%-4s %10s / %s   %5.1f%%\n" "$d" "$s" "$NSTEPS" \
          "$(awk -v a="$s" -v b="$NSTEPS" 'BEGIN{print 100*a/b}')"
      fi
    done
    echo ""
    echo "$left of 6 still to run"
    nvidia-smi --query-gpu=index,name,utilization.gpu --format=csv,noheader 2>/dev/null
    ;;
  stop)
    touch "$FLAG"
    pkill -f "mdrun -deffnm $DEFF"
    sleep 3
    pgrep -af "mdrun -deffnm $DEFF" || echo "stopped"
    ;;
  *)
    echo "usage: bash curve2_run.sh run | pair N | status | stop"
    echo "       STOP_AT=18:00 bash curve2_run.sh run"
    ;;
esac
