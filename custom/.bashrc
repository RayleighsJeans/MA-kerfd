# .bashrc

# Source global definitions
#if [ -f /etc/bashrc ]; then
#	. /etc/bashrc
#fi

# Uncomment the following line if you don't like systemctl's auto-paging feature:
# export SYSTEMD_PAGER=

# User specific aliases and functions

alias ls='ls -lah --color=auto'
alias vi='vim'

PS1='\e[1;92m \d \t >\e[m\n\u@\H:\w \e[1;92m\$>\e[m'
