#!/bin/bash
#SBATCH -N 1            # number of nodes
#SBATCH -n 1            # number of "tasks" (default: 1 core per task)
#SBATCH -t 7-00:00:00   # time in d-hh:mm:ss
#SBATCH --mem=15000mb
#SBATCH -o ../slurms/sc2-full-test.slurm.%j.out # file to save job's STDOUT (%j = JobId)
#SBATCH -e ../slurms/sc2-full-test.slurm.%j.err # file to save job's STDERR (%j = JobId)
#SBATCH --mail-type=ALL # Send an e-mail when a job starts, stops, or fails
#SBATCH --mail-user=----@gmail.com # Mail-to address
#SBATCH --job-name=sc2-full-test_230111


#turn on debugging output
set -x;

echo "Original folder contents:";
ls;
echo "untar slim";
tar -xzf slim.tar.gz;
echo "untar python libraries";
tar -xzf python_env.tar.gz
echo "All folder contents:";
ls * -ar;

echo "slim folder contents:";
ls build/;

export PATH=$_CONDOR_SCRATCH_DIR/build:$PATH;
export PYTHONPATH=$PWD/python_env;



#define parameter space for Model 0
MU_dist=("2.135e-7" "2.135e-6" "2.135e-5");
#REPRO_dist=("1" "3" "5");
K_dist=("5e3" "5e4" "1e5");
INIT_dist=("1" "5" "100");
RUNTIME_dist=("336" "672" "1008");

#define parameter space for Model 1
R_dist=("1e-5" "5.5e-5" "10e-5");

#define parameter space for Model 3
XI_dist=("0.0001" "0.003" "0.1");
BURSTN_dist=("20" "100" "200");

#define parameter space for Model 7
DFE_dist=("c(4, 1)" "c(1,1)" "c(1, 4)");



#Load all arguments into one vector; only the parameter selections into another
arguments=(${@});
paramSel=(${@:2});

echo "My params: ${arguments[@]}";

#prepare to interpret what options being run
#get the model number from commandline
modelNum=${arguments[0]};
modelChew=$modelNum;
maxModel=$((2 ** 4 - 1));
modelOpt=();
moreModel=1;

#convert model number into boolean vector
while (($moreModel))
do
	modelOpt=(${modelOpt[@]} $(expr $modelChew % 2));
	modelChew=$(($modelChew / 2));
	maxModel=$(($maxModel / 2));
	if (($maxModel <= 0))
	then
		moreModel=0;
	fi
done

#begin building slim command
slimCommand='slim';
slimCommand=$slimCommand' -d GENOMESIZE=30000 -d MODEL='${arguments[0]};

#prepare to interpret required parameters being run
#first element of paramSel corresponds to required parameters
reqdNum=${paramSel[0]};
reqdChew=$reqdNum;
maxReqd=$((3 ** 4 - 1));
moreReqd=1;
reqdOpt=();
#convert required parameter number into base-3 vector
while(($moreReqd))
do
	reqdOpt=(${reqdOpt[@]} $(expr $reqdChew % 3));
	reqdChew=$(($reqdChew / 3));
	maxReqd=$(($maxReqd / 3));
	if (($maxReqd <= 0));
	then
		moreReqd=0;
	fi
done

echo ${reqdOpt[@]};

#build slim command using base-3 vector
slimCommand=$slimCommand' -d MU='${MU_dist[${reqdOpt[0]}]};
slimCommand=$slimCommand' -d REPRO=1';
slimCommand=$slimCommand' -d K='${K_dist[${reqdOpt[1]}]};
slimCommand=$slimCommand' -d INIT='${INIT_dist[${reqdOpt[2]}]};
slimCommand=$slimCommand' -d RUNTIME='${RUNTIME_dist[${reqdOpt[3]}]};

#if model number uses recombination, add that to slim command
if ((${modelOpt[0]}))
then
	slimCommand=$slimCommand' -d R='${R_dist[${paramSel[1]}]}
	#slimCommand=$slimCommand' -d R='${R_dist[$3]}
fi

#if model number uses progeny skew, add that to slim command
if ((${modelOpt[1]}))
then
	pSkOpt=${paramSel[2]};
	
	slimCommand=$slimCommand' -d XI='${XI_dist[$((pSkOpt % 3))]};
	slimCommand=$slimCommand' -d BURSTN='${BURSTN_dist[$(((pSkOpt / 3) % 3))]};
fi

#if model number uses DFE, add that to slim command
if ((${modelOpt[2]}))
then
	slimCommand=$slimCommand" -d CLASSN=2 -d 'DFDISTR=c(\"return 0.0;\", \"return -1.0;\")' -d 'DFCLASS=${DFE_dist[${paramSel[3]}]}'"
fi

#add burnin options
slimCommand=$slimCommand" -d 'BURNINPARAMS=c(1, 1e5, 2.5e4, \"sc2-full2_burnin25k.slout\")'";

#add output stem
outputStem="sc2-full2_"$modelNum"."${paramSel[0]}"."${paramSel[1]}"."${paramSel[2]}"."${paramSel[3]}"."${paramSel[4]};
slimCommand=$slimCommand" -d \"OUTPUTSTEM='"${outputStem}"'\"";

#add slim script
slimCommand=$slimCommand" sc2-full2_sc2-sim-v0.2.slim";

#call slim command
echo $slimCommand;
#eval $slimCommand;

echo "simulation run.";
python3 sc2-full2_stats-v0.2.py ${outputStem}"_partial_100 m 30000"
python3 sc2-full2_stats-v0.2.py ${outputStem}"_partial_1000 m 30000"
echo "folder contents now:";
ls;
