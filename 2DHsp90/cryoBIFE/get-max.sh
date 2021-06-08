#### Bash script to get Gmax
### run: ./

nstep=`awk '$1=="Steps"' $1 | sort -nrk4 | head -1 | awk '{print $3}'`

awk -v nstep=$nstep '{if(NR==nstep)for(x=2;x<=NF;x++)print $x}' Gs > Gmax-$1

