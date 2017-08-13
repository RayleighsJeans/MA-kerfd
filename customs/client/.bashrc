#
# ~/.bashrc
#

# If not runnig interactiveley, don't do anything
[[ $- != *i* ]] && return


alias ls='ls -lah --color=auto --group-directories-first --time-style="+%F, %T "'
alias vi='vim'
alias clean-tex="bash ~/.clean-tex.sh"

alias brain='ssh brain'
alias brain-wX11='ssh brain -Y' 

alias pdf='epdfview'

PS1='\e[1;92m\u@\H >\e[m\n \e[1;34m\W\e[m \e[1;92m\$ >\e[m'

if [ "$(tty)" = "/dev/tty1" ]; then
	startx
fi
