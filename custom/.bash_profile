# .bash_profile

# Get the aliases and functions
if [ -f ~/.bashrc ]; then
	. ~/.bashrc
fi

# User specific environment and startup programs

PATH=$PATH:$HOME/.local/bin:$HOME/bin

export PATH

#PS1=' \d \t>\n\u@\H:\w\$>'

#alias ls='ls --color=auto'
PS1='\e[1;92m \d \t >\e[m\n\u@\H:\w \e[1;92m\$>\e[m'

if [ "$(tty) = /dev/tty1" ]; then
	cd Documents/2D/pic2d/
	module load arcanist/03112016
	module load cmake/3.4.3
	module load superlu/4.3
	module load clang/4.0
	module load munge/0.5.11
	module load slurm/15.08.8
	module load gcc/5.2.0
	module load openmpi/1.10.0/gcc/5.2.0
	module load hdf5-ompi1100-gcc520/1.8.15-p1
fi

