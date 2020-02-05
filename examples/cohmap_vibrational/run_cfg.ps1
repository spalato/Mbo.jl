$cfg=$args
if (!$cfg) {
    $cfg=$(ls *.yaml)
}
foreach ($c in $cfg) {
    julia vibrational_ia.jl $c
    python ../plot_linear.py $c
    python ../plot_2d.py $c
    python ./cohmap_vib.py $c 
}
