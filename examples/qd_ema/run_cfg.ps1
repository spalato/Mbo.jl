$cfg=$args
if (!$cfg) {
    $cfg=$(ls *.yaml)
}
foreach ($c in $cfg) {
    julia --project ema_BX.jl $c
    python ../plot_linear.py $c
    python ../plot_2d.py $c
    python plot_peak_dyn.py $c
}