#!/bin/bash

today=$(date +"%Y-%m-%d")
# Define the directory path where you want to store the flag file
flags_dir="./CR/output/flags"
echo "Done1"
# Create the directory if it doesn't exist
mkdir -p "$flags_dir"
if [ -f "$flags_dir/completion_flag_${today}" ]; then
  echo "removing"
    rm "$flags_dir/completion_flag_${today}"
fi
echo "over"

# Array to store the PIDs of background jobs
pids=()
echo pids

for l in 1 #{1..2} # cause 1 probabilty for CR
  do for k in 1 #{1..2}  #crit(both surv and endpoint = stick together), 1-2
    do for j in {1..2}  #size = 1,2
      do for i in {1..2} #propensity = 1,2
        do for h in 1 #{1..2} #{1..4} beta = 1,2,3,4
          #ncauses = 2
   	        #endpoint = CR
              do sbatch -t 07-05:00:00  --mail-type=fail --mail-user=cwzhou@email.unc.edu --mem=9000 CR_S0run.sh CR01.Simulation_Run.R 1 1 $h $i $j $k $k $l $1  #$1 = date
      	    #done
      	  #done
        done
      done
    done
  done
done

echo "Jobs submitted"
#names(arg)[1:7] = c("endpoint", "ncauses", "beta", "propensity", "size", "crit_surv", "crit_endpoint")
# do this for arg8 = 1 and then arg8 = 2 (change with CR02 script right now)
echo "${pids[@]}"

# Wait for all background jobs to complete
for pid in "${pids[@]}"; do
    wait "$pid"
done

echo "wait over"

# Create the completion flag file in the specified directory with today's date as part of the filename
touch "$flags_dir/completion_flag_${today}"