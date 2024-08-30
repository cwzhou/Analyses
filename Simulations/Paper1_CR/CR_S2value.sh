for l in {1..2} # cause 1 probabilty for CR # change to just 1 for simple_exp
  do for k in 1 #{1..2}  #crit(both surv and endpoint = stick together), 1-2
    do for j in {1..2}  #size = 1,2
      do for i in {1..2} #propensity = 1,2
        do for h in {1..2} #{1..4} beta = 1,2,3,4
          #ncauses = 2
            do for g in {1..2} # censor is either low (1) or high (2)
   	      #endpoint = CR
                do sbatch -t 07-05:00:00  --mail-type=fail --mail-user=cwzhou@email.unc.edu --mem=9000 CR_S0run.sh CR01.Simulation_Run.R 1 $g 1 $h $i $j $k $k $l $1  #$1 = date
      	      #done
            done
      	  #done
        done
      done
    done
  done
done
#names(arg)[1:9] = c("endpoint", "censor", "ncauses", "beta",
#                    "propensity", "size", "crit_surv",
#                    "crit_endpoint", "cause1prob")
# do this for arg8 = 1 and then arg8 = 2 (change with CR02 script right now)
