#
# ~/.bashrc
#

# If not runnig interactiveley, don't do anything
[[ $- != *i* ]] && return

alias ls='ls -lah --color=auto'
alias vi='vim'
alias hpc-access='ssh hpc-access'
alias brain='ssh brain'

PS1='\e[1;92m\d \t>\e[m\n \u@\H: \w \e[1;92m\$>\e[m'

#if [ "$(tty)" = "/dev/tty1" ]; then
#	startx
#fi
