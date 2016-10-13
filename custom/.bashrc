#
# ~/.bashrc
#

# If not runnig interactiveley, don't do anything
[[ $- != *i* ]] && return


alias ls='ls --color=auto'
PS1='\e[1;92m \d \t >\n\u@\H:\e[m\w \e[1;92m\$>\e[m'
