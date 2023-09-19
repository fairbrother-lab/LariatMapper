#!/bin/bash

# String multiplier from https://newbedev.com/multiplying-strings-in-bash-script

# Paramter 1: A string 
# Parameter 2: The number of line breaks to add at the end of the string
# Returns a string with the time since seconds concatenated at the beginning and the current time concatenated at the end
timestamp(){
    time_elapsed="$(($SECONDS / 3600))h:$(($SECONDS % 3600 / 60))m:$(($SECONDS % 3600 % 60))s"     #If the time is less than 10,000 hours, the max length of this string is 28
	time_now="$(date '+%H:%M:%S %m/%d')"
    time="$time_now  -  $time_elapsed"
    time_len=${#time}
    
    space_for_spaces="$((31-$time_len -1))"s     # -1 to account for a space in statement="$time $1"
    spacer=$(printf "%$space_for_spaces")
    spaces=${spacer// /' '}
    time="$time |${spaces}"

    statement="$time $1"

    if [ $# == 2 ] && [ $2 != 0 ];then
        line_breaking="$2s"
        line_breaker=$(printf "%$line_breaking")
        line_breaks=${line_breaker// /'\n'}
		# If $2 is positive, add \n's to the back. If it's negative, add them to the front
		if [[ $2 -gt 0 ]];then
        	statement="$statement$line_breaks"
		elif [[ $2 -lt 0 ]];then
			statement="${line_breaks}${statement}"
		fi
    fi
    
    echo -e "${statement}"
}

