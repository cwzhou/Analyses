for k in {1..2}  #crit(both surv and endpoint = stick together), 1-2
  do for j in 1  #size = 1,2
    do for i in 1  #propensity = 1,2
      #beta = one version
        #ncauses = 2
          #endpoint = CR
            do sbatch -t 07-09:00:00  --mail-type=fail --mail-user=cwzhou@email.unc.edu --mem=10000 CR_S0run.sh CR01.Simulation_Parameters.R 1 1 1 $i $j $k $k $1  #$1 = date
          #done
        #done
      #done
    done
  done
done

