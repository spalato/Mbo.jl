$cfg=$args
if (!$cfg) {
    $cfg=$(ls *.yaml)
}
foreach ($c in $cfg) {
    julia kubo.jl $c
    python ../plot_linear.py $c
    python ../plot_2d.py $c
}