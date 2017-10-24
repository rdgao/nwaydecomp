# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# User specific aliases
alias matlab='matlab -singleCompThread'
alias matlabnd='matlab -singleCompThread -nodesktop'
alias ll='ls -laFh'
alias lll='ls -lFh'
alias qstatme='qstat -u roevdmei'


# set editor
export EDITOR=vim

# load modules
module load matlab/2015a