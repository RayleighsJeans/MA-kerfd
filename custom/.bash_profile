#
# ~/.bash_profile
#

[[ -f ~/.bashrc ]] && . ~/.bashrc

alias vim="vim -S ~/.vimrc"
alias vi='vim'
alias ls='ls -lah --color=auto'
alias hpc-access='ssh hpc-access'
alias brain='ssh brain'

#if [ "$(tty)" = "/dev/tty1" ]; then
#	startx
#fi
	
