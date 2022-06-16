#!/bin/bash

sigma=$1
T=$2
N=$3
q=$4
a=$5
r=$6
i=$7
nRepeat=$8
uid=$9

[[ $# -ne 9 ]] && echo "Error: You need to give 9 arguments." && exit 1


pc0475=0.628244 #FROM OAS





case $q in
	0) pc=$pc0 ;;
	0.1) pc=$pc01 ;;
	0.2) pc=$pc02 ;;
	0.3) pc=$pc03 ;;
	0.4) pc=$pc04 ;;
	0.41) pc=$pc041 ;;
	0.42) pc=$pc042 ;;
	0.43) pc=$pc043 ;;
	0.44) pc=$pc044 ;;
	0.45) pc=$pc045 ;;
	0.46) pc=$pc046 ;;
	0.462) pc=$pc0462 ;;
	0.464) pc=$pc0464;;
	0.466) pc=$pc0466;;
	0.468) pc=$pc0468;;
	0.469) pc=$pc0469;;
	0.47) pc=$pc047 ;;
	0.471) pc=$pc0471 ;;
	0.472) pc=$pc0472 ;;
	0.473) pc=$pc0473 ;;
	0.474) pc=$pc0474;;
	0.475) pc=$pc0475;;
	0.476) pc=$pc0476;;
	0.478) pc=$pc0478;;
	0.48) pc=$pc048 ;;
	0.49) pc=$pc049 ;;
	0.5) pc=$pc05 ;;
	0.51) pc=$pc051 ;;
	0.52) pc=$pc052 ;;
	0.53) pc=$pc053 ;;
	0.54) pc=$pc054 ;;
	0.56) pc=$pc056 ;;
	0.57) pc=$pc057 ;;
	0.58) pc=$pc058 ;;
	0.59) pc=$pc059 ;;
	0.6) pc=$pc06 ;;
	0.62) pc=$pc062 ;;
	0.64) pc=$pc064 ;;
	0.66) pc=$pc066 ;;
	0.67) pc=$pc067 ;;
	0.68) pc=$pc068 ;;
	0.69) pc=$pc069 ;;
	0.7) pc=$pc07 ;;
	0.71) pc=$pc071 ;;
	0.72) pc=$pc072 ;;
	0.73) pc=$pc073 ;;
	0.74) pc=$pc074 ;;
	0.75) pc=$pc075 ;;
	0.76) pc=$pc076 ;;
	0.77) pc=$pc077 ;;
	0.78) pc=$pc078 ;;
	0.79) pc=$pc079 ;;
	0.8) pc=$pc08 ;;
	0.82) pc=$pc082 ;;
	0.84) pc=$pc084 ;;
	0.85) pc=$pc085 ;;
	0.86) pc=$pc086 ;;
	0.88) pc=$pc088 ;;
	0.884) pc=$pc0884 ;;
	0.886) pc=$pc0886 ;;
	0.89) pc=$pc089 ;;
	0.9) pc=$pc09 ;;
	0.91) pc=$pc091 ;;
	0.92) pc=$pc092 ;;
	*)
		echo "Error: q_c is unknown for sigma=$sigma"
		exit 2
		;;
esac

p=$(echo "" | awk '{printf("%lf", (pc + a * r**i))}' pc=$pc a=$a r=$r i=$i)
qq=$(echo "" | awk '{printf("%.3f", (q))}' q=$q)
sigmaa=$(echo "" | awk '{printf("%.2f", (sigma))}' sigma=$sigma)


params=${sigmaa}_${T}_${N}_${qq}_${a}_${r}_${i}_${nRepeat}_${uid}
logFile=Log/LTCP_2D_$params.table
outputFile=Result/LTCP_2D_$params.table
#outputFile2=Result/CP_Dynamics_$params.table


exec 2>> $logFile

echo "$(date +%Y-%m-%d_%H:%M:%S): start $$ $(hostname)" >> $logFile

#time Bin/DP_ContactModelSteady $p $T $N >> $outputFile #2>> $outputFile2

Bin/LTCP_2D $p $q $sigma $T $N $nRepeat > $outputFile #2> $outputFile2
ret=$?;
echo -n "$(date +%Y-%m-%d_%H:%M:%S): " >> $logFile
if [ $ret -eq 0 ]; then
	echo "normal exit"
else
	echo "abnormal exit: $ret"
fi >> $logFile

#for((t=1; t<$nRepeat; t++)); do
#	binn/DP_ContactModelSteady $p $T $N #>> $outputFile #2>> $outputFile2
#	#sleep 1
#done >> $outputFile #&2>> $logFile
#echo "$(date +%Y-%m-%d_%H:%M:%S): end $$ $(hostname)" >> $logFile


















