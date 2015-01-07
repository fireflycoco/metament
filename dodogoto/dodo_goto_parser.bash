#!/bin/bash 

gotodatafile=~/.dodo_goto.data
touch $gotodatafile

onhelp()
{
	echo "Usage:"
	echo 'alias g="source PATH_TO_DODOGOTO/dodo_goto_parser.bash"'
	echo 'Add a tag'
	echo '    g add TAG'
	echo 'List all tags'
	echo '    g'
	echo 'Jump to TAG'
	echo '    g TAG'
	echo 'Delete a tag'
	echo '    g del TAG'
	echo 'Delete all tags'
	echo '    g delall'
}

search() 
{
	bookmark=$1
	res=`grep "^${bookmark}:" $gotodatafile`
	if [ "$res" != "" ] ; then
		echo $res | cut -d ":" -f 2-
	fi
}

if [ $# == 0 ] ; then
	cat $gotodatafile | sort | sed -e "s/:/:\t/"
elif [ $# == 1 ] ; then
	if [ "$1" == "--help" ] || [ "$1" == "-h" ] ; then
		echo "`onhelp`"
	elif [ "$1" == "delall" ] ; then
		echo "Delete all? yes/no"
		read answer
		if [ "$answer" ==  "yes" ] ; then
			echo "Deleting everything."
			rm -f $gotodatafile
			touch $gotodatafile
		fi
	else
		pathname=`search $1`
		if [ "$pathname" != "" ] && [ "$pathname" != "." ] ; then
			cd "$pathname"
		fi
	fi
elif [ $# == 2 ] ; then
	if [ $1 == "del" ] ; then
		pathname=`search $2`
		if [ "${pathname}" != "" ] ; then
			echo "Deleting ${2}: 	${pathname}"
			sed -i -e "/^${2}:/d" $gotodatafile
		fi
	elif [ $1 == 'add' ] ; then
		sed -i -e "/^${2}:/d" $gotodatafile
		pathname="`pwd`"
		echo Adding "${2}: 	${pathname}"
		echo "${2}:${pathname}" >> $gotodatafile
	fi
fi

