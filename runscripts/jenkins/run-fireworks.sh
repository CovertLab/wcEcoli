set -e

# Commented rapidfire command below produces seg fault after 2 hr and 10 min (see #764)
# Could replace singleshot loop with rapidfire if fixed
# Singleshot might seg fault as well for long single tasks over 2 hr and 10 min

# rlaunch rapidfire --nlaunches 0 2> >(tee -a stderr.log >&2)
while [ $(lpad get_fws -s READY -d count) -ge 1 ]; do
  rlaunch singleshot 2> >(tee -a stderr.log >&2)
done

N_FAILS=$(lpad get_fws -s FIZZLED -d count)

if [ $N_FAILS -gt 0 ]; then
  echo -e "\nSims failed. $(git status | head -1) $(git rev-parse --short HEAD):"
  lpad get_fws -s FIZZLED
  cat stderr.log
  mv out/2* /scratch/PI/mcovert/wc_ecoli/failed/
fi

test $N_FAILS = 0
