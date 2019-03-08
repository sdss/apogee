foreach file ( $* )
set n = `grep lsf $file | wc -l`
echo $file $n
grep lsf $file | awk "NR==int($n/2)"
end

