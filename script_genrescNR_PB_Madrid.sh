#!/bin/bash
#$ -cwd

#script_genrescNR_PB.sh <d> <LAMB> <BETA> <AVEs> <AVEh> <SEED>
#To be used only in a machine with /state/partition1 directory

rm script_genrescNR_PB.sh.*

##########################################################################################
# Outline                                                                                #
#                                                                                        #
#           __ 'NINDTP'R'NINDTP' (line with same size as threatened population and GR)   #
#          |                                                                             #
#          |__ 'NINDTP'R'NINDAL' (line with different size and GR)                       #
# NINDTP---|                                                                             #
#          |__ 'NINDTP'NR'NINDAL' (line with different size and without GR)              #
#          |                                                                             #
#          |                  __ 'NINDTP'RR'NINDTP'                                      #
#          |                 |                       (line with same size as threatened  #
#          |__ 'NINDTP'NR  __|__ 'NINDTP'RR'NINDAL'  population without GR, and          #
#                            |                       subsequent formation of sublines)   #
#                            |__ 'NINDTP'NR'NINDTP'                                      #
#                            |                                                           #
#                            |__ 'NINDTP'NR'NINDAL'                                      #
#                                                                                        #
##########################################################################################


####################################### ARGUMENTS ######################################

#Check number of arguments
if [ $# -ne 6 ]  
then
	echo "Usage: $0 <d> <LAMB> <BETA> <AVEs> <AVEh> <SEED>" 
	exit 1
fi

#Set arguments
d=$1
LAMB=$2
BETA=$3
AVEs=$4
AVEh=$5
SEED=$6

############################### VARIABLES AND DIRECTORIES ##############################

#Parameters
Vs=0
Opt=0
NR=0

NCRO=300
NINDTP=1376
NINDAL=43
GEN=250
tR=83
MIG=1
INT=1

m_TP=0
m_AL=0
dist=0

#Working directory
WDIR=$PWD 
mkdir -p $WDIR/genrescPB_results/L$LAMB.k$NCRO.s$AVEs.h$AVEh.beta$BETA/TP$NINDTP.AL$NINDAL.M$MIG.INT$INT.g$tR
DIR="genrescPB_results/L$LAMB.k$NCRO.s$AVEs.h$AVEh.beta$BETA/TP$NINDTP.AL$NINDAL.M$MIG.INT$INT.g$tR"

#Scratch directory
mkdir -p /state/partition1/noeliaGR$d/$SLURM_JOBID/
SCDIR="/state/partition1/noeliaGR$d" 

############################# TRANSFER OF FILES TO SCRATCH #############################

#Copy all files in scratch directory
cp seedfile $SCDIR/$SLURM_JOBID/
cp genrescNR_PB $SCDIR/$SLURM_JOBID/

#File with information of node and directory
touch $WDIR/$SLURM_JOBID.`hostname`.`date +%HH%MM`

#Move to scratch directory
cd $SCDIR/$SLURM_JOBID

####################################### GENRESC ########################################

START=$(date +%s)
time ./genrescNR_PB>>out<<@
0
$SEED
10000	NIND Base Population
$NINDTP	NIND Threatened Population
$NINDAL	NIND Alternative Lines
$m_TP	Migrants for lines with NINDTP individuals
$m_AL	Migrants for lines with NINDAL individuals
0	Gender of Migrants (males 0, males&females 1)
99	Lenght genome (99=free)
$NCRO	NCRO (max 2000)(Neu=Ncro)
30	NLOCI (2-30)
$LAMB	Lambda_s
0.0	Lambda_L
$BETA	beta_s
$AVEs	ave |s|
1	dom model (0=cnt; 1:variable)
$AVEh	ave h
$Vs	Stabilizing selection (Vs)
1	VE
$Opt	Optimal
0	neutral (0: no, 1:yes)
$GEN	generations
$tR	Formation of R and NR lines (max 'generations')
$tR	Formation of sublines ('tLINES' - 'generations')
$NR	For R lines, initial NR-generations after formation of line
$MIG	Number of migrations (99: periodic)
$INT	Generation intervals of migration (lines NINDTP)
$INT	Generation intervals of migration (lines NINDAL)
$dist	Generations since migration (distribution -w- among replicates)
0.0	Minimum fitness value for extinction (99: without extinction)
100	Replicates
@

END=$(date +%s)
DIFF=$(( $END - $START ))
echo "genresc took 		$DIFF seconds" >> timefile

###################### TRANSFER OF FILES TO MAIN DIRECTORY ############################

cp -r $SCDIR/$SLURM_JOBID/seedfile $WDIR
cp -r $SCDIR/$SLURM_JOBID/timefile $WDIR/$DIR/
cp -r $SCDIR/$SLURM_JOBID/dfilename.dat $WDIR/$DIR/
cp -r $SCDIR/$SLURM_JOBID/natpop.dat $WDIR/$DIR/
cp -r $SCDIR/$SLURM_JOBID/mutation.dat $WDIR/$DIR/
#cp -r $SCDIR/$SLURM_JOBID/genfile_L1.dat $WDIR/$DIR/genfileL1_"$NINDTP"R"$NINDTP".dat
#cp -r $SCDIR/$SLURM_JOBID/genfile_L2.dat $WDIR/$DIR/genfileL2_"$NINDTP"R"$NINDAL".dat
cp -r $SCDIR/$SLURM_JOBID/genfile_L3.dat $WDIR/$DIR/genfileL3_"$NINDTP"NR"$NINDAL".dat
#cp -r $SCDIR/$SLURM_JOBID/genfile_L4subl1.dat $WDIR/$DIR/genfileL4sub1_"$NINDTP"RR"$NINDTP".dat
#cp -r $SCDIR/$SLURM_JOBID/genfile_L4subl2.dat $WDIR/$DIR/genfileL4sub2_"$NINDTP"RR"$NINDAL".dat
cp -r $SCDIR/$SLURM_JOBID/genfile_L4subl3.dat $WDIR/$DIR/genfileL4sub3_"$NINDTP"NR"$NINDTP".dat
#cp -r $SCDIR/$SLURM_JOBID/genfile_L4subl4.dat $WDIR/$DIR/genfileL4sub4_"$NINDTP"NR"$NINDAL".dat
cp -r $SCDIR/$SLURM_JOBID/distribution_qsh.dat $WDIR/$DIR/
cp -r $SCDIR/$SLURM_JOBID/summary_outline.dat $WDIR/$DIR/
cp -r $SCDIR/$SLURM_JOBID/out $WDIR/$DIR/
#cp -r $SCDIR/$SLURM_JOBID/rescfile.dat $WDIR/$DIR/
cp -r $SCDIR/$SLURM_JOBID/extinction.dat $WDIR/$DIR/
#cp -r $SCDIR/$SLURM_JOBID/distribution_w.dat $WDIR/$DIR/
cp -r $SCDIR/$SLURM_JOBID/repfile.dat $WDIR/$DIR/

############################# CLEANING OF SCRATCH ####################################

rm -r $SCDIR/$SLURM_JOBID/
rm $WDIR/$SLURM_JOBID.*


