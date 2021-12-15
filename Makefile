
world: train_deriv validate

train_deriv: ./derivation/train.R
	nohup Rscript ./derivation/train.R T > ./derivation/logs/train-uti.log 2>&1
	nohup Rscript ./derivation/train.R F > ./derivation/logs/train-le.log 2>&1

validate: ./validation_study.R ./derivation/results/revised_trained_uti.RData ./derivation/results/revised_trained_le.RData
	nohup Rscript ./validation_study.R > ./logs/validate.log 2>&1
