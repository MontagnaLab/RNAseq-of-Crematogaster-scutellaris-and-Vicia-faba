for i in $(cat $1); do if grep $i $2; then echo -n ""; else echo $i ; fi; done
