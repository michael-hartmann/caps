echo "# L/R, L, R, T, ldim, F*(L+R)/(ħc)"

for file in pc_gk/slurm-*.out; do
    echo `grep -v "#" $file`
done
