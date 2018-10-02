for i in *K/; do echo $i; cd $i; sleep 1; for j in grainindex_vs_atomindex_from_train*soft; do aa=`echo $j | cut -c30-99`; ../a.out $aa $j;  done; cd ..; echo ""; done
