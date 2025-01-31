"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
source ~/.bashrc
echo "Creating conda environment"
micromamba create -f $1/environment.yml --yes
echo "Activating conda environment"
micromamba activate fixienv